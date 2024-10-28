#!/usr/bin/env python3
import os
import subprocess
import threading
import WDL
import json
import re
import pandas as pd
import numpy as np
from pathlib import Path

import gcsfs
import requests
from google.cloud import storage

from copy import deepcopy

import cromwell.initialize as ini
import cromwell.preprocess_run as pre
from cromwell.constants import *
import cromwell.run_monitor as crm

from datetime import datetime
import dateutil
import time
import traceback
from zoneinfo import ZoneInfo


class CromwellManager:
    def __init__(self, run_name, inputs_file, json_template_path, wdl_path, save_specific_outputs=[],
                 batch=None, limit=None, n_parallel_workflows=N_PARALLEL_WORKFLOWS,
                 add_requester_pays_parameter=True, restart=False, batches_precomputed=False,
                 submission_sleep=30, check_freq=120, quiet=False):
        """
        Initialize a Cromwell manager.

        PARAMETERS
        ----------
        run_name: name of run
        inputs_file: pandas flat file with column names to be added to the JSON for submission.
            Each row will be a job unless "batch" is set in which case rows will be grouped.
        json_template_path: local path to json template file. Will be augmneted based on inputs_file.
        wdl_path: local path to WDL to submit
        save_specific_outputs: list. Can be empty. If it contains items, will track only those specific outputs.
        batch: How many rows to add to each batch. If None, no batching is done.
        n_limit: number to submit
        n_parallel_workflows: number of workflows to allow at once
        add_requester_pays_parameter: if true will add req pays as a parameter for WDL
        restart: if true, will ignore prior state found in the inputs file
        batches_precomputed: if true, will rely on "batch" column to identify batches
        addl_sub_interval_sec:
        submission_retries:
        submission_sleep: time to sleep between each submission
        check_freq: how many seconds to wait before checking again for new submissions
        cromwell_timeout:
        """
        # initialize cromwell
        ini.main(quiet)
        self.app_name, self.app_status, self.app_url, self.env = ini.get_cromwell_url()
        self.cloud_fs = gcsfs.GCSFileSystem(project=PROJECT, requester_pays=True)
        
        # run parameters
        self.run_name = run_name
        self.output = os.path.join(BUCKET, run_name)
        self.samples_df_output_path = os.path.join(self.output, 'samples_df.csv')
        self.pipeline_status_output_path = os.path.join(self.output, 'pipeline_status.tsv')
        self.wdl_path = wdl_path
        self.batch = batch
        self.n_parallel_workflows = n_parallel_workflows
        self.submission_sleep = submission_sleep
        self.check_frequency = check_freq
        self.quiet = quiet

        # preliminary run statistics
        self.n_pending = 0
        self.n_submitted = 0
        self.n_running = 0
        self.n_success = 0
        self.n_fail = 0
        self.n_review = 0

        # process json with parameters
        with open(json_template_path) as json_in:
            json_template = json.load(json_in)

        json_template, json_prefix = pre.verify_json_template(json_template, quiet=quiet)
        self.json_template = json_template
        self.workflow_name = json_prefix

        if add_requester_pays_parameter:
            self.update_template('requester_pays_project', PROJECT)

        # process input flat file
        if limit is not None:
            inputs_file = inputs_file.head(limit)
        inputs_file = pre.initialize_with_status(inputs_file, restart)
        if batches_precomputed and (np.sum(inputs_file.columns == 'cromwell_id') == 1):
            print('Batches precomputed, column located.', flush=True)
        else:
            inputs_file = pre.compute_batches(inputs_file, batch)
        self.sample_file = inputs_file

        # create placeholder for workflow status file
        self.initialize_workflow_status_file(save_specific_outputs, restart)

        # create parameters file
        run_params = self.get_initial_parameters_json()
        self.write_parameters_json(run_params)

        # update table based on what is running in Cromwell
        self.running_workflows = {}
        self.rebuild_running_jobs() # since just starting, we need to remake workflow objects
        self.update_all_status() # then check to see if these have updated status


    def run_pipeline(self, submission_retries, cromwell_timeout, skip_waiting=False):
        # Primary function to submit and monitor all jobs from this pipeline.
        # This spawns two processes. The first will submit all jobs.
        # The second will update a df containing run information per workflow.

        # submit the initial set of jobs and then periodically check to submit more as they complete, up to 40k samples
        self.update_run_statistics()
        self.lock = threading.Lock()

        self.exit_signal = threading.Event()

        thread_queue = threading.Thread(target=self.queue_all_jobs, args=(submission_retries, cromwell_timeout), daemon=True)
        thread_monitor = threading.Thread(target=self.monitor_for_job_metrics, daemon=True)

        thread_queue.start()
        thread_monitor.start()
        print(f'Pipeline {self.run_name} launched.', flush=True)

        if not skip_waiting:
            try:
                while not self.exit_signal.is_set():
                    if not thread_queue.is_alive() and not thread_monitor.is_alive():
                        break
                    time.sleep(1)
            except KeyboardInterrupt:
                print('Quitting all threads...', flush=True)
                self.exit_signal.set()
            
            thread_queue.join()
            thread_monitor.join()
            print(f'Pipeline {self.run_name} completed.', flush=True)
        else:
            print('ALERT: pipeline threads have been initialized and may be running in the background.', flush=True)


    def queue_all_jobs(self, submission_retries, cromwell_timeout):
        # This function handles submission of all jobs.
        try:
            with self.lock:
                this_pending = self.n_pending
            
            while (this_pending > 0) and not self.exit_signal.is_set():
                #reload submission parameters each iteration to allow them to be tuned over the course of the run
                with self.lock:
                    self.update_parameters_from_disk()
                
                time_now = datetime.now(tz=ZoneInfo('America/New_York')).isoformat()
                if not self.quiet:
                    print(f'{time_now}: Checking to submit jobs...', flush=True)

                #Now, do the submission
                with self.lock:
                    _ = self.submit_jobs(addl_sub_interval_sec=self.submission_sleep,
                                         submission_retries=submission_retries,
                                         cromwell_timeout=cromwell_timeout)
                    check_frequency = self.check_frequency

                # wait before updating the job status and attempting to submit more jobs
                time.sleep(check_frequency)

                with self.lock:
                    # update the status results
                    self.update_all_status()
                    this_pending = self.n_pending

                    # check whether we are having excessive failures (reload status_counts to get updated numbers)
                    if (self.n_fail + self.n_success) > 0:
                        pct_failed = (self.n_fail/(self.n_fail+self.n_success))*100

                        # if so, stop the loop
                        if (self.n_fail > MIN_FAIL_NUM) and (pct_failed > MAX_FAIL_PCT):
                            msg = (f'{pct_failed:0.1f}% of terminated jobs have ended in failure, which is greater\n'
                                f'than the threshold setting of {MAX_FAIL_PCT}%. Halting the submission loop.\n'
                                'Please check to see if something is wrong.')
                            with open('./PIPELINE_SUBMISSION_ERROR', 'w') as out:
                                out.write(msg + '\n')
                            print(msg, flush=True)
                            break
            
            else:
                print('All samples submitted. Submission loop terminating.', flush=True)

        except Exception:
            excpt_str = traceback.format_exc()
            print(excpt_str, flush=True)
            with open('./PIPELINE_SUBMISSION_ERROR', 'w') as out:
                out.write(excpt_str + '\n')
            raise
        

    def monitor_for_job_metrics(self):

        print('Starting monitoring. Monitoring thread will run until there are no more jobs in a non-terminal state.', flush=True)
        
        with self.lock:
            run_metrics_holder = {k: [] for k in self.workflow_status.columns}
            this_workflow_name = self.workflow_name

        while not self.exit_signal.is_set():

            with self.lock:
                self.update_all_status()
                if not self.quiet:
                    self.print_status(flush=True)
                n_non_terminal = self.n_running + self.n_pending + self.n_submitted
                final_check = n_non_terminal == 0
                
                # Get completed workflow IDs
                completed_workflow_ids = list(self.get_samples_with_status('Succeeded').cromwell_id.unique())
                # If there was an existing run_metrics file, filter for just the workflow_ids that haven't been 
                # recorded there yet.
                completed_workflow_ids = [elt for elt in completed_workflow_ids if elt not in self.workflow_status['cromwell_id'].values]
            
            # now get the results for all of the newly-succeeded workflows
            if not self.quiet and (len(completed_workflow_ids) > 0):
                print('Collating information on newly completed workflows...', flush=True)
            for idx, workflow_id in enumerate(completed_workflow_ids):
                if idx and not idx%10:
                    if not self.quiet:
                        print(f'Processed {idx} newly completed workflows', flush=True)
                attempts = 4
                to_raise = None
                while attempts >= 0:
                    try:
                        run_meta_resp = subprocess.run(['cromshell', '-t', '60', '--no_turtle', '--machine_processable', 'metadata', 
                                                        '--dont-expand-subworkflows', workflow_id], 
                                                        check=True, capture_output=True)
                    except Exception as err:
                        print(f'Error retrieving info about workflow {workflow_id}. Retrying {attempts} more times.', flush=True)
                        to_raise = err
                        attempts -= 1
                        time.sleep(15)
                    else:
                        break
                else:
                    raise to_raise
                
                run_meta = json.loads(run_meta_resp.stdout.decode())
                run_metrics = run_metrics_holder.copy()

                # the usual stuff
                run_metrics['cromwell_id'].append(workflow_id)
                run_metrics['status'].append(run_meta['status'])
                run_metrics['start_time'].append(run_meta['start'])
                run_metrics['end_time'].append(run_meta['end'])
                run_metrics['runtime'].append(str(dateutil.parser.isoparse(run_meta['end'])
                                                - dateutil.parser.isoparse(run_meta['start'])))
                
                outputs_to_find = [k for k in run_metrics.keys() if re.search(f'^{this_workflow_name}.+', k)]
                for output in outputs_to_find:
                    run_metrics[output].append(run_meta['outputs'].get(output, 'Not found'))
                
                with self.lock:
                    self.workflow_status = pd.concat([self.workflow_status, pd.DataFrame(run_metrics)]).reset_index(drop=True)

            if len(completed_workflow_ids) > 0:
                with self.lock:
                    self.save_pipeline_status_file()
            
            if final_check:
                print('All workflows have completed.', flush=True)
                break
            else:
                time.sleep(self.check_frequency)


    def submit_jobs(self, addl_sub_interval_sec, submission_retries, cromwell_timeout):
        """
        Performs one round of submissions.
        """
        self.update_run_statistics()
        if self.n_pending == 0:
            print('No samples with a status of "Submission pending". Nothing new to submit.', flush=True)
        else:
            #group the samples ready for submission into batches of size batch_size, record the inputs, and do the submission
            submission_records_root = os.path.join(self.output, 'cromwell_submissions')
            
            #to submit
            to_submit_df = self.get_samples_with_status('Submission pending')
            ct_to_submit = max(min(self.n_pending, self.n_parallel_workflows - self.n_running), 0)
            print(f'Found {str(self.n_pending)} samples awaiting submission, and {str(self.n_running)} '
                'currently running samples. Given the provided submission and concurrency limits, '
                f'will submit {str(ct_to_submit)} items.', flush=True)
            if self.batch is not None:
                n_batches = ct_to_submit // self.batch
                first_batch = list(to_submit_df.head(1).batch)[0]
                print(f'Since batch mode is enabled, will submit {str(n_batches)} of size {str(self.batch)}, starting from batch {str(first_batch)}.', flush=True)
                ct_to_submit = n_batches * self.batch
            else:
                n_batches = ct_to_submit

            batch_list = to_submit_df.batch.unique()
            batch_to_submit = batch_list[0:min(n_batches, len(batch_list))]
            to_submit_this_round = to_submit_df[to_submit_df['batch'].isin(batch_to_submit)]
            if self.batch is None:
                to_submit_by_batch = to_submit_this_round
            else:
                to_submit_by_batch = to_submit_this_round.groupby('batch').agg(list).reset_index()

            for _, row in to_submit_by_batch.iterrows():
                col_names = [x for x in row.index if x not in STATUS_COLS]
                this_batch = row.batch

                if self.batch is None:
                    batch_root = submission_records_root
                    dct_update = {x: row[x] for x in col_names}
                else:    
                    batch_root = os.path.join(submission_records_root, f'batch_{str(this_batch)}')
                    if not self.quiet:
                        print('Saving file lists for use with batched Cromwell imports.', flush=True)
                        print('Note: this involves transforming input names:', flush=True)
                    dct_update = {}
                    for col in col_names:
                        this_uri = os.path.join(batch_root, f'{col}_list.txt')
                        with self.cloud_fs.open(this_uri, 'w') as out:
                            out.write('\n'.join(row[col]) + '\n')
                        if not self.quiet:
                            print(f'{col} -> {col}_list', flush=True)
                        dct_update.update({f'{col}_list': this_uri})
                
                run_specific_json = self.create_run_specific_parameters(dct_update)

                # create workflow
                workflow = CromwellWorkflow(run_specific_json, batch_root=batch_root, wdl=self.wdl_path,
                                            manager_url=self.get_url(), manager_token=self.get_token(),
                                            cloud_fs=self.cloud_fs, timeout=cromwell_timeout,
                                            n_retries=submission_retries, batch_id=this_batch, quiet=self.quiet)

                # add workflow to running samples
                self.add_running_workflow(workflow)

                #write the info about the submitted jobs back to the full samples_df dataframe
                self.update_status_by_col('batch', this_batch, 'cromwell_id', workflow.get_id())
                self.update_status_by_col('batch', this_batch, 'run_status', workflow.get_status())
                self.update_status_by_col('batch', this_batch, 'inputs_uri', batch_root)

                # update sample file on GCP
                self.save_sample_file()
                
                # wait between batch submissions, if requested
                time.sleep(addl_sub_interval_sec)
            
            print('Submission round complete.', flush=True)
            return self.samples_df_output_path
    

    def resubmit_failed_workflows(self, addl_sub_interval_sec, submission_retries, cromwell_timeout, specific_id=None):
        """
        If specific_id is not None, will resubmit a specific job. Otherwise will resubmit all jobs labeled "Failed".
        """
        if specific_id is None:
            to_resubmit = self.get_samples_with_status('Failed')
        else:
            to_resubmit = self.sample_file.loc[self.sample_file['cromwell_id' == specific_id]].copy()

        if to_resubmit.shape[0] == 0:
            print('Nothing to resubmit.')
        else:
            print(f'Resubmitting {str(to_resubmit.shape[0])} workflows...')

            for _, row in to_resubmit.iterrows():
                this_batch = row.batch
                this_id = row.cromwell_id
                this_batch_root = row.inputs_uri
                with self.cloud_fs.open(os.path.join(this_batch_root, 'input_parameters.json'), 'r') as j:
                    run_specific_json = json.load(j)

                #preserve the old sub_id.txt in preparation for resubmission and making a new one
                prev_sub_id_uri_glob = os.path.join(this_batch_root, 'prev_sub_id*')
                try:
                    res = subprocess.run(['gsutil', '-u', PROJECT, 'ls', prev_sub_id_uri_glob], check=True, capture_output=True)
                except subprocess.CalledProcessError as err:
                    if err.stderr.decode().strip() == 'CommandException: One or more URLs matched no objects.':
                        prev_sub_idx = 0
                    else:
                        raise err
                else:
                    prev_sub_id_list = res.stdout.decode().strip().split('\n')
                    prev_sub_idx = len(prev_sub_id_list)
                prev_sub_id_uri = os.path.join(this_batch_root, f'prev_sub_id{prev_sub_idx}.txt')
                sub_id_uri = os.path.join(this_batch_root, 'sub_id.txt')
                subprocess.run(['gsutil', '-q', '-u', PROJECT, 'mv', sub_id_uri, prev_sub_id_uri], check=True)

                #submit the batch as a cromwell workflow
                workflow = CromwellWorkflow(run_specific_json, batch_root=this_batch_root, wdl=self.wdl_path,
                                            manager_url=self.get_url(), manager_token=self.get_token(),
                                            cloud_fs=self.cloud_fs, timeout=cromwell_timeout,
                                            n_retries=submission_retries, batch_id=this_batch, quiet=self.quiet)
                
                # add workflow to running samples
                self.add_running_workflow(workflow)

                #write the info about the submitted jobs back to the full samples_df dataframe
                self.update_status_by_col('batch', this_batch, 'cromwell_id', workflow.get_id())
                self.update_status_by_col('batch', this_batch, 'run_status', workflow.get_status())
                self.update_status_by_col('batch', this_batch, 'inputs_uri', this_batch_root)
                
                print(f'Batch {this_id} resubmitted with workflow ID {workflow.get_id()}.')

                # update sample file on GCP
                self.save_sample_file()

                # wait between batch submissions, if requested
                time.sleep(addl_sub_interval_sec)


    def rebuild_running_jobs(self):
        samples_to_check = self.get_samples_with_status('Running')
        if len(samples_to_check) > 0:
            to_check = samples_to_check[['batch','cromwell_id','inputs_uri']].drop_duplicates()
            for _, row in to_check.iterrows():
                this_batch = row.batch
                this_id = row.cromwell_id
                this_batch_root = row.inputs_uri

                with self.cloud_fs.open(os.path.join(this_batch_root, 'input_parameters.json'), 'r') as j:
                    run_specific_json = json.load(j)

                workflow = CromwellWorkflow(run_specific_json, batch_root=this_batch_root, wdl=self.wdl_path,
                                            manager_url=self.get_url(), manager_token=self.get_token(),
                                            cloud_fs=self.cloud_fs, timeout=None,
                                            n_retries=None, batch_id=this_batch, id=this_id, status='Running', quiet=self.quiet)
                
                self.add_running_workflow(workflow)
            self.update_run_statistics()


    def update_all_status(self):
        self.update_run_statistics()
        if not self.quiet:
            print(f'Updating status for {str(self.n_running)} workflow ids.', flush=True)
        
        this_running_workflows = self.running_workflows.copy()
        for idx, (id, workflow) in enumerate(self.running_workflows.items()):
            if (idx % 20) == 0 and not self.quiet:
                print(f'{idx} workflows checked.', flush=True)

            if workflow.get_id() != id:
                raise ValueError('ERROR: ID mismatch when updating workflow status.')
            
            this_updated_status = workflow.update_status()
            if this_updated_status is not None:
                if this_updated_status in COMPLETED_STATUS:
                    del this_running_workflows[id]
                self.update_status_by_col('cromwell_id', id, 'run_status', this_updated_status)
            
            self.save_sample_file()
        
        self.running_workflows = this_running_workflows
        self.update_run_statistics()


    def add_running_workflow(self, workflow):
        self.running_workflows.update({workflow.get_id(): workflow})


    def update_template(self, key, value):
        self.json_template.update({f'{self.workflow_name}.{key}': value})


    def create_run_specific_parameters(self, dct):
        json_temp = deepcopy(self.json_template)
        json_temp.update({f'{self.workflow_name}.{k}': v for k, v in dct.items()})
        return(json_temp)


    def initialize_workflow_status_file(self, save_specific_outputs, restart):

        def custom_import_resolver(uri, path, importer):
            if uri.startswith("https://"):
                # Download the file from the remote URL
                path_this = os.path.abspath("./wdl_imports")
                if not os.path.isdir(path_this):
                    os.mkdir(path_this)
                
                local_path = os.path.join(path_this, os.path.basename(uri))
                if not self.quiet:
                    print(f'Moving {os.path.basename(local_path)} to local fs from networked location for parsing...')
                response = requests.get(uri)
                with open(local_path, 'wb') as f:
                    f.write(response.content)
                return WDL.read_source_default(local_path, path, importer)
            else:
                # Use the default importer for non-https URIs
                return WDL.read_source_default(uri, path, importer)

        this_wdl = WDL.load(self.wdl_path, read_source=custom_import_resolver)
        status = {'cromwell_id':[],
                  'status':[],
                  'start_time':[],
                  'end_time':[],
                  'runtime':[]}
        if len(save_specific_outputs) > 0:
            dict_outputs = {f'{self.workflow_name}.{item.name}': [] for item in this_wdl.workflow.outputs if item in save_specific_outputs}
        else:
            dict_outputs = {f'{self.workflow_name}.{item.name}': [] for item in this_wdl.workflow.outputs}
        status.update(dict_outputs)
        
        if self.cloud_fs.exists(self.pipeline_status_output_path) and not restart:
            this_df = pd.read_csv(self.pipeline_status_output_path, sep='\t',
                                  storage_options={'project':PROJECT, 'requester_pays':True})
            tf1 = all([k in this_df.columns for k in status.keys()])
            tf2 = all([col in status.keys() for col in this_df.columns])
            if tf1 and tf2:
                self.workflow_status = this_df  
            else:
                print('Found workflow status file with columns: ' + ', '.join([k for k in this_df.columns]))
                print('Expected workflow status file with columns:' + ', '.join([col for col in dict_outputs.keys()]))
                raise ValueError(f'Found a file at {self.pipeline_status_output_path} which does not contain the exact expected schema. Either remove this file or specify restart=True.')
        else:
            self.workflow_status = pd.DataFrame(status)


    def update_run_statistics(self):
        # adds statistics about what items have each of a variety of statuses
        self.n_pending = self.get_samples_with_status('Submission pending').cromwell_id.unique().shape[0]
        self.n_submitted = self.get_samples_with_status('Submitted').cromwell_id.unique().shape[0]
        self.n_running = self.get_samples_with_status('Running').cromwell_id.unique().shape[0]
        self.n_success = self.get_samples_with_status('Succeeded').cromwell_id.unique().shape[0]
        self.n_fail = self.get_samples_with_status('Failed').cromwell_id.unique().shape[0]
        self.n_review = self.get_samples_with_status('Manual review').cromwell_id.unique().shape[0]

        if (self.n_running + self.n_submitted) != len(self.running_workflows):
            raise ValueError('ERROR: it does not make sense for running workflows in the sample table to be different than recorded running workflows.')


    def print_status(self, flush):
        print('============== STATUS REPORT ===============', flush=flush)
        print(f'{str(self.n_submitted)} jobs are submitted but not yet running.', flush=flush)
        print(f'{str(self.n_running)} jobs are running.', flush=flush)
        print(f'{str(self.n_pending)} jobs are awaiting submission.', flush=flush)
        print(f'{str(self.n_success)} jobs have succeeded.', flush=flush)
        print(f'{str(self.n_fail)} jobs have failed.', flush=flush)
        print('============================================', flush=flush)


    def get_samples_with_status(self, status):
        assert(status in ALLOWED_STATUS)
        return self.sample_file.loc[self.sample_file['run_status'] == status].copy()


    def update_status_by_col(self, filter_col, filter_val, target_col, update_to):
        new_sample_file = self.sample_file.copy()
        if filter_val not in new_sample_file[filter_col].unique():
            raise ValueError('ERROR: value not found in sample table.')
        new_sample_file.loc[new_sample_file[filter_col] == filter_val, target_col] = update_to
        self.sample_file = new_sample_file


    def get_url(self):
        return self.app_url
    

    def get_token(self):
        return self.env['token']
     

    def get_initial_parameters_json(self):
        # set pipeline submission parameter values
        pipe_params = {'RUN_NAME': self.run_name,
                       'BATCH_SIZE': 1 if self.batch is None else self.batch,
                       'NUM_CONCURRENT': self.n_parallel_workflows,
                       'MIN_FAIL_NUM': MIN_FAIL_NUM,
                       'MAX_FAIL_PCT': MAX_FAIL_PCT,
                       'SUBMISSION_WAIT': self.submission_sleep,
                       'SUB_CHECK_WAIT': self.check_frequency,
                       'PIPELINE_STATUS_URI': self.pipeline_status_output_path}

        return pipe_params


    def update_parameters_from_disk(self):
        # Supports updating the number of concurrent workflows, submission wait, check frequency while the pipeline runs.
        params_uri = os.path.join(self.output, 'pipeline_submission_params.json')
        with self.cloud_fs.open(params_uri, 'r') as j:
            data = json.load(j)

        self.n_parallel_workflows = data['NUM_CONCURRENT']
        self.submission_sleep = data['SUBMISSION_WAIT']
        self.check_frequency = data['SUB_CHECK_WAIT']


    def write_parameters_json(self, pipe_params):
        #save to local fs
        params_path = './pipeline_submission_params.json'
        with open(params_path, 'w') as out:
            json.dump(pipe_params, out)
        #save to cloud
        params_uri = os.path.join(self.output, 'pipeline_submission_params.json')
        with self.cloud_fs.open(params_uri, 'w') as out:
            json.dump(pipe_params, out)
        
        print(f'To change parameter values for the running submission process, edit the file here:\n{params_uri}')


    def save_sample_file(self):
        self.sample_file.to_csv(self.samples_df_output_path, sep='\t', index=False,
                                storage_options={'project':PROJECT, 'requester_pays':True})


    def save_pipeline_status_file(self):
        self.workflow_status.to_csv(self.pipeline_status_output_path, sep='\t', index=False,
                                    storage_options={'project':PROJECT, 'requester_pays':True})


class CromwellWorkflow:
    def __init__(self, param_json, batch_root, wdl,
                 manager_url, manager_token, cloud_fs,
                 timeout, n_retries, batch_id,
                 id=None, status=None, quiet=False):                 
        self.json = param_json
        self.batch_root = batch_root
        self.token = manager_token
        self.batch_id = batch_id

        if id is not None or status is not None:
            # if ID is supplied, this implies that this was a running workflow
            # thus we don't need to restart it
            if id is None or status is None:
                raise ValueError('Must provide both ID and status when constructing CromwellWorkflow.')
            self.id = id
            self.status = status
            self.status_url = self.assemble_status_url(manager_url)
        
        else:
            this_n_retries = n_retries
            
            # dump json
            with cloud_fs.open(self.get_json_path(), 'w') as out:
                json.dump(self.json, out)

            local_json = './input_parameters.json'
            with open(local_json, 'w') as out:
                json.dump(self.json, out)

            # submit job
            cmd_arg_list = ['cromshell', '-t', str(timeout), '--no_turtle', '--machine_processable', 'submit', wdl, local_json]
            batch_submission_cmd_uri = os.path.join(self.batch_root, 'batch_submission_cmd.txt')
            with cloud_fs.open(batch_submission_cmd_uri, 'w') as out:
                out.write(' '.join(cmd_arg_list) + '\n')
            while True:
                try:
                    sub_resp = subprocess.run(cmd_arg_list, check=True, capture_output=True)
                except subprocess.CalledProcessError as e:
                    if this_n_retries == 0:
                        print(f'Error output (code {str(e.returncode)}):')
                        print((e.stderr).decode('utf-8'))
                        print('Command: ' + ' '.join(cmd_arg_list))
                        raise e
                    this_n_retries -= 1
                else:
                    break
            sub_info = json.loads(sub_resp.stdout.decode())
            sub_id_txt_uri = os.path.join(batch_root, 'sub_id.txt')

            self.id = sub_info['id']
            self.status = sub_info['status']
            self.status_url = self.assemble_status_url(manager_url)

            with cloud_fs.open(sub_id_txt_uri, 'w') as out:
                out.write(sub_info['id'] + '\n')

            if not quiet:
                print(f'Batch {str(int(self.batch_id))} submitted ({self.get_id()}).')


    def update_status(self):
        status_json = self.check_status()
        if status_json['status'] != self.status:
            self.status = status_json['status']
            return self.status
        else:
            return None


    def get_json_path(self):
        return os.path.join(self.batch_root, 'input_parameters.json')


    def get_id(self):
        return self.id


    def get_status(self):
        return self.status


    def check_status(self, detailed=False, cromwell_run_prefix='cromwell-execution'):
        if detailed:
            bucket_id = os.path.basename(BUCKET)
            storage_client = storage.Client()
            bucket = storage_client.get_bucket(bucket_id)
            return crm.check_success_single(storage_client, bucket, self.id, cromwell_run_prefix, False)
        
        else:
            r = requests.get(
                self.status_url,
                timeout=5,
                verify=True,
                headers={
                    'Referer': 'https://notebooks.firecloud.org',
                    'Authorization': f'Bearer {self.token}'
                }
            )
            r.raise_for_status()
            return r.json()


    def assemble_status_url(self, manager_url):
        return f'{manager_url}/api/workflows/v1/{self.id}/status'


def isolate_failed_shards(samples_df_uri, failed_workflows_list, cromwell_run_prefix='cromwell-execution'):
    #get samples_df
    samples_df = pandas.read_csv(samples_df_uri, index_col=0,
                                 storage_options={'project':os.getenv('GOOGLE_PROJECT'), 'requester_pays':True})
    
    # iterate over workflow IDs and change any samples without a FAIL status in the Cromwell run
    # to "Submission pending" status in the samples_df so that they can be packaged into a new
    # batch and re-submitted
    for crom_id in failed_workflows_list:
        #check that this is a valid workflow id
        if crom_id not in samples_df['cromwell_id'].values:
            print(f'Workflow ID {crom_id} is not found. Skipping.')
            continue

        #get relevant subset of the samples_df
        workflow_samples = samples_df.loc[samples_df['cromwell_id'] == crom_id].copy()

        #check to be sure these are failed samples
        workflow_status = workflow_samples['run_status'].unique()
        if len(workflow_status) != 1 or workflow_status[0] != 'Failed':
            print(f'Workflow {crom_id} has a status other than "Failed" (status: {", ".join(workflow_status)}). Skipping.')
            continue

        #get sample submission order for the requested workflow
        samples_order_uri = os.path.join(workflow_samples['inputs_uri'].unique()[0], 'sample_name_list.txt')
        samples_order = pandas.read_csv(samples_order_uri, header=None, 
                                        storage_options={'project':os.getenv('GOOGLE_PROJECT'),
                                                         'requester_pays':True})
        
        #get the table of metadata about the shards of this workflow
        run_info = get_detailed_workflow_status(crom_id, cromwell_run_prefix=cromwell_run_prefix)[0]
        
        #identify the failed shards and update the status of the non-failed ones so they will 
        # be resubmitted in a new batch
        failed_idx = [int(elt.split('-')[-1]) for elt in run_info.loc[run_info['status'] == 'FAIL', 'shard'].values]
        to_resub_idx = workflow_samples.loc[~workflow_samples['person_id'].isin(samples_order.iloc[failed_idx, 0].to_list())].index
        samples_df.loc[to_resub_idx, 'run_status'] = 'Submission pending'
        to_rev_idx = workflow_samples.loc[workflow_samples['person_id'].isin(samples_order.iloc[failed_idx, 0].to_list())].index
        samples_df.loc[to_rev_idx, 'run_status'] = 'Manual review'

    #save the results back to the cloud
    samples_df.to_csv(samples_df_uri, storage_options={'project':os.getenv('GOOGLE_PROJECT'), 'requester_pays':True})
    return samples_df

