#!/usr/bin/env python3
import os
import subprocess
import json
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

from datetime import datetime, timedelta
import dateutil
import time
import traceback
from zoneinfo import ZoneInfo


class CromwellManager:
    def __init__(self, run_name, inputs_file, json_template_path, wdl_path,
                 batch=None, limit=None, n_parallel_workflows=N_PARALLEL_WORKFLOWS,
                 add_requester_pays_parameter=True, restart=False, batches_precomputed=False,
                 submission_sleep=30, check_freq=120):
        """
        Initialize a Cromwell manager.

        PARAMETERS
        ----------
        run_name: name of run
        inputs_file: pandas flat file with column names to be added to the JSON for submission.
            Each row will be a job unless "batch" is set in which case rows will be grouped.
        json_template_path: local path to json template file. Will be augmneted based on inputs_file
        wdl_path: local path to WDL to submit
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
        ini.main()
        self.app_name, self.app_status, self.app_url, self.env = get_cromwell_url()
        self.cloud_fs = gcsfs.GCSFileSystem(project=PROJECT, requester_pays=True)
        
        # run parameters
        self.run_name = run_name
        self.output = os.path.join(BUCKET, run_name)
        self.samples_df_output_path = os.path.join(self.output, 'samples_df.csv')
        self.wdl_path = wdl_path
        self.batch = batch
        self.n_parallel_workflows = n_parallel_workflows
        self.submission_sleep = submission_sleep
        self.check_frequency = check_freq

        # preliminary run statistics
        self.n_pending = 0
        self.n_running = 0
        self.n_success = 0
        self.n_fail = 0
        self.n_review = 0

        # process json with parameters
        with open(json_template_path) as json_in:
            json_template = json.load(json_in)

        json_template, json_prefix = pre.verify_json_template(json_template)
        self.json_template = json_template
        self.workflow_name = json_prefix

        if add_requester_pays_parameter:
            self.update_template('requester_pays_project', PROJECT)

        # process input flat file
        if limit is not None:
            inputs_file = inputs_file.head(limit)
        inputs_file = pre.initialize_with_status(inputs_file, restart)
        if batches_precomputed and (np.sum(inputs_file.columns == 'cromwell_id') == 1):
            print('Batches precomputed, column located.')
        else:
            inputs_file = pre.compute_batches(inputs_file, batch)
        self.sample_file = inputs_file

        # create parameters file
        run_params = self.get_initial_parameters_json()
        self.write_parameters_json(run_params)

        # update table based on what is running in Cromwell
        self.running_workflows = {}
        self.rebuild_running_jobs() # since just starting, we need to remake workflow objects
        self.update_all_status() # then check to see if these have updated status


    def submit(self, addl_sub_interval_sec, submission_retries, cromwell_timeout):
        """
        Performs one round of submissions.
        """
        self.update_run_statistics()
        if self.n_pending == 0:
            print('No samples with a a status of "Submission pending". Nothing new to submit.')
        else:
            #group the samples ready for submission into batches of size batch_size, record the inputs, and do the submission
            submission_records_root = os.path.join(self.output, 'cromwell_submissions')
            
            #to submit
            to_submit_df = self.get_samples_with_status('Submission pending')
            ct_to_submit = max(min(self.n_pending, self.n_parallel_workflows - self.n_running), 0)
            print(f'Found {str(self.n_pending)} samples awaiting submission, and {str(self.n_running)} '
                'currently running samples. Given the provided submission and concurrency limits, '
                f'will submit {str(ct_to_submit)} items.')
            if self.batch is not None:
                n_batches = ct_to_submit // self.batch
                first_batch = list(to_submit_df.head(1).batch)[0]
                print(f'Since batch mode is enabled, will submit {str(n_batches)} of size {str(self.batch)}, starting from batch {str(first_batch)}.')
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
                col_names = [x for x in row.index if x != 'batch']
                this_batch = row.batch

                if self.batch is None:
                    batch_root = submission_records_root
                    dct_update = {x: row[x] for x in col_names}
                else:    
                    batch_root = os.path.join(submission_records_root, f'batch_{str(this_batch)}')
                    print('Saving file lists for use with batched Cromwell imports.')
                    print('Note: this involves transforming input names:')
                    dct_update = {}
                    for col in col_names:
                        this_uri = os.path.join(batch_root, f'{col}_list.txt')
                        with self.cloud_fs.open(this_uri, 'w') as out:
                            out.write('\n'.join(row[col]) + '\n')
                        print(f'{col} -> {col}_list')
                        dct_update.update({f'{col}_list': this_uri})
                
                run_specific_json = self.create_run_specific_parameters(dct_update)
                print(batch_root)

                # create workflow
                workflow = CromwellWorkflow(run_specific_json, batch_root=batch_root, wdl=self.wdl_path,
                                            manager_url=self.get_url(), manager_token=self.get_token(),
                                            cloud_fs=self.cloud_fs, timeout=cromwell_timeout,
                                            n_retries=submission_retries, batch_idx=this_batch)

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
            
            print('Submission complete.')
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
                                            n_retries=submission_retries, batch_idx=this_batch)
                
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
                                            n_retries=None, batch_idx=this_batch, id=this_id, status='Running')
                
                self.add_running_workflow(workflow)
            self.update_run_statistics()


    def update_all_status(self):
        self.update_run_statistics()
        print(f'Updating status for {str(self.n_running)} workflow ids.')
        
        for idx, (id, workflow) in enumerate(self.running_workflows.items()):
            if (idx % 20) == 0:
                print(f'{idx} workflows checked.')

            if workflow.get_id() != id:
                raise ValueError('ERROR: ID mismatch when updating workflow status.')
            
            this_updated_status = workflow.update_status()
            if this_updated_status is not None:
                self.update_status_by_col('cromwell_id', id, 'run_status', this_updated_status)
            self.save_sample_file()


    def add_running_workflow(self, workflow):
        self.running_workflows.update({workflow.get_id(): workflow})


    def update_template(self, key, value):
        self.json_template.update({f'{self.workflow_name}.{key}': value})


    def create_run_specific_parameters(self, dct):
        json_temp = deepcopy(self.json_template)
        json_temp.update({f'{self.workflow_name}.{k}': v for k, v in dct.items()})
        return(json_temp)


    def update_run_statistics(self):
        # adds statistics about what items have each of a variety of statuses
        self.n_pending = self.get_samples_with_status('Submission pending').shape[0]
        self.n_running = self.get_samples_with_status('Running').shape[0]
        self.n_success = self.get_samples_with_status('Succeeded').shape[0]
        self.n_fail = self.get_samples_with_status('Failed').shape[0]
        self.n_review = self.get_samples_with_status('Manual review').shape[0]

        if self.n_running != len(self.running_workflows):
            raise ValueError('ERROR: it does not make sense for running workflows in the sample table to be different than recorded running workflows.')


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
        pipeline_status_path = f'./{self.run_name}.pipeline_status.tsv'
        pipeline_status_uri = os.path.join(self.output, f'{self.run_name}.pipeline_status.tsv')

        pipe_params = {'RUN_NAME': self.run_name,
                       'BATCH_SIZE': 1 if self.batch is None else self.batch,
                       'NUM_CONCURRENT': self.n_parallel_workflows,
                       'MIN_FAIL_NUM': MIN_FAIL_NUM,
                       'MAX_FAIL_PCT': MAX_FAIL_PCT,
                       'SUBMISSION_WAIT': self.submission_sleep,
                       'SUB_CHECK_WAIT': self.check_frequency,
                       'PIPELINE_STATUS_PATH': pipeline_status_path,
                       'PIPELINE_STATUS_URI': pipeline_status_uri}

        return pipe_params


    def update_parameters_from_disk(self):
        # Supports updating the number of concurrent workflows, submission wait, check frequency while the pipeline runs.
        params_uri = os.path.join(self.out, 'pipeline_submission_params.json')
        with self.cloud_fs.read(params_uri, 'r') as j:
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
        params_uri = os.path.join(self.out, 'pipeline_submission_params.json')
        with self.cloud_fs.open(params_uri, 'w') as out:
            json.dump(pipe_params, out)
        print(f'To change parameter values for the running submission process, edit the file here:\n{params_uri}')


    def save_sample_file(self):
        self.sample_file.to_csv(self.samples_df_output_path, sep='\t', index=False,
                                storage_options={'project':PROJECT, 'requester_pays':True})


class CromwellWorkflow:
    def __init__(self, param_json, batch_root, wdl,
                 manager_url, manager_token, cloud_fs,
                 timeout, n_retries, batch_id,
                 id=None, status=None):
        self.json = param_json
        self.batch_root = batch_root
        self.status_url = self.assemble_status_url(manager_url)
        self.token = manager_token
        self.batch_id = batch_id

        if id is not None or status is not None:
            # if ID is supplied, this implies that this was a running workflow
            # thus we don't need to restart it
            if id is None or status is None:
                raise ValueError('Must provide both ID and status when constructing CromwellWorkflow.')
            self.id = id
            self.status = status
        
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
                except subprocess.CalledProcessError:
                    if this_n_retries == 0:
                        raise
                    this_n_retries -= 1
                else:
                    break
            sub_info = json.loads(sub_resp.stdout.decode())
            sub_id_txt_uri = os.path.join(batch_root, 'sub_id.txt')

            self.id = sub_info['id']
            self.status = sub_info['status']

            with cloud_fs.open(sub_id_txt_uri, 'w') as out:
                out.write(sub_info['id'] + '\n')

            print(f'Batch {str(self.batch_id)} submitted ({self.get_id()}).')


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


def get_succeeded_job_metrics(samples_df_uri, run_name, force_reload=False):
    run_metrics_uri = os.path.join(os.getenv("WORKSPACE_BUCKET"), f'{run_name}/run_metrics.csv')
    run_metrics = {'cromwell_id':[],
                   'status':[],
                   'start_time':[],
                   'end_time':[],
                   'runtime':[],
                   'merging_log':[],
                   'merged_calls':[],
                   'merged_coverage':[],
                   'merged_statistics':[]}
    #get samples_df
    samples_df = pandas.read_csv(samples_df_uri, index_col=0,
                                 storage_options={'project':os.getenv('GOOGLE_PROJECT'), 'requester_pays':True})
    workflow_ids = samples_df.loc[samples_df['run_status'] == 'Succeeded', 'cromwell_id'].unique()
    #get run_metrics (if it already exists), otherwise just make it an empty version of the one we are building
    if not force_reload:
        try:
            existing_run_metrics = pandas.read_csv(run_metrics_uri, index_col=0, sep='\t',
                                                   storage_options={'project':os.getenv('GOOGLE_PROJECT'), 
                                                                    'requester_pays':True})
        except:
            existing_run_metrics = pandas.DataFrame(run_metrics)
        else:
            #if there is an existing run_metrics file, filter for just the workflow_ids that haven't been 
            # recorded there yet.
            workflow_ids = [elt for elt in workflow_ids if elt not in existing_run_metrics['cromwell_id'].values]
    else:
        existing_run_metrics = pandas.DataFrame(run_metrics)

    # now get the results for all of the newly-succeeded workflows
    for idx, workflow_id in enumerate(workflow_ids):
        if idx and not idx%5:
            print(f'Processed {idx} workflows')
        attempts = 4
        to_raise = None
        while attempts >= 0:
            try:
                run_meta_resp = subprocess.run(['cromshell', '-t', '60', '--no_turtle', '--machine_processable', 'metadata', 
                                                '--dont-expand-subworkflows', workflow_id], 
                                                check=True, capture_output=True)
            except Exception as err:
                print(f'Error retrieving info about workflow {workflow_id}. Retrying {attempts} more times.')
                to_raise = err
                attempts -= 1
                time.sleep(15)
            else:
                break
        else:
            raise to_raise
        run_meta = json.loads(run_meta_resp.stdout.decode())
        run_metrics['cromwell_id'].append(workflow_id)
        run_metrics['status'].append(run_meta['status'])
        run_metrics['start_time'].append(run_meta['start'])
        run_metrics['end_time'].append(run_meta['end'])
        run_metrics['runtime'].append(str(dateutil.parser.isoparse(run_meta['end'])
                                          - dateutil.parser.isoparse(run_meta['start'])))
        try:
            run_metrics['merging_log'].append(run_meta['calls']['MitochondriaPipelineWrapper.MergeMitoMultiSampleOutputsInternal'][-1]['backendLogs']['log'])
        except KeyError:
            run_metrics['merging_log'].append('Not found')
        run_metrics['merged_calls'].append(run_meta['outputs'].get('MitochondriaPipelineWrapper.merged_calls', 'Not found'))
        run_metrics['merged_coverage'].append(run_meta['outputs'].get('MitochondriaPipelineWrapper.merged_coverage', 'Not found'))
        run_metrics['merged_statistics'].append(run_meta['outputs'].get('MitochondriaPipelineWrapper.merged_statistics', 'Not found'))

    run_metrics = pandas.concat([existing_run_metrics, pandas.DataFrame(run_metrics)]).reset_index(drop=True)
    run_metrics.to_csv('./run_metrics.csv', sep='\t')
    run_metrics.to_csv(run_metrics_uri, sep='\t',
                       storage_options={'project':os.getenv('GOOGLE_PROJECT'), 'requester_pays':True})
    print(run_metrics_uri)
    return run_metrics


def run_pipeline():
    #Now, submit all workflows until 40k samples have terminated
    pipeline_status = {'timestamp':[],
                       'Submission pending':[],
                       'Submitted':[],
                       'Running':[],
                       'Succeeded':[],
                       'Failed':[]}

    # submit the initial set of jobs and then periodically check to submit more as they complete, up to 40k samples
    try:
        while numpy.sum(test_run['run_status'] != 'Submission pending') < 1000:
            #reload submission parameters each iteration to allow them to be tuned over the course of the run
            with cloud_fs.open(params_uri, 'r') as params_in:
                params_json = json.load(params_in)
            run_name = params_json['RUN_NAME']
            batch_size = params_json['BATCH_SIZE']
            num_concurrent = params_json['NUM_CONCURRENT']
            min_fail_num = params_json['MIN_FAIL_NUM']
            max_fail_pct = params_json['MAX_FAIL_PCT']
            submission_wait_sec = params_json['SUBMISSION_WAIT']
            sub_check_wait_min = params_json['SUB_CHECK_WAIT']
            pipeline_status_path = params_json['PIPELINE_STATUS_PATH']
            pipeline_status_uri = params_json['PIPELINE_STATUS_URI']

            #Now, do the submission
            test_run_uri = pss.submit_cromwell_workflows(test_run, run_name=run_name, batch_size=batch_size, 
                                                         mtSwirl_root='/home/jupyter/mtSwirl_fork/mtSwirl/',
                                                         num_concurrent_crams=num_concurrent,
                                                         addl_sub_interval_sec=submission_wait_sec)

            time_now = datetime.now(tz=ZoneInfo('America/New_York')).isoformat()
            print(time_now)

            # update and record the summary of the pipeline status
            status_counts = test_run['run_status'].value_counts()
            status_counts_dict = status_counts.to_dict()
            pipeline_status['timestamp'].append(time_now)
            for k in sorted(set(pipeline_status.keys()) | set(status_counts_dict.keys())):
                if k == 'timestamp':
                    continue
                try:
                    pipeline_status[k].append(status_counts_dict.get(k, 0))
                except KeyError:
                    #any less-frequent pipeline status types will be added to the dataframe on the fly
                    pipeline_status[k] = ([0]*(len(pipeline_status['timestamp'])-1)) + [status_counts_dict[k]]
            pipeline_status_df = pandas.DataFrame(pipeline_status)
            pipeline_status_df.to_csv(pipeline_status_path, sep='\t', index=False)
            pipeline_status_df.to_csv(pipeline_status_uri, sep='\t', index=False,
                                      storage_options={'project':os.getenv('GOOGLE_PROJECT'), 'requester_pays':True})

            # wait before updating the job status and attempting to submit more jobs
            time.sleep(60*sub_check_wait_min)

            # update the status results
            test_run = pss.update_cromwell_status(test_run_uri, verbose=False)
            test_run['person_id'] = test_run['person_id'].astype(str)

            # check whether we are having excessive failures (reload status_counts to get updated numbers)
            status_counts = test_run['run_status'].value_counts()
            status_counts_dict = status_counts.to_dict()
            success_count = status_counts_dict.get('Succeeded', 0)
            failure_count = status_counts_dict.get('Failed', 0)
            pct_failed = (failure_count/(failure_count+success_count))*100

            # if so, stop the loop
            if (failure_count > min_fail_num) and (pct_failed > max_fail_pct):
                msg = (f'{pct_failed:0.1f}% of terminated jobs have ended in failure, which is greater\n'
                       f'than the threshold setting of {max_fail_pct}%. Halting the submission loop.\n'
                       'Please check to see if something is wrong.')
                with open('./PIPELINE_SUBMISSION_ERROR', 'w') as out:
                    out.write(msg + '\n')
                print(msg)
                break
        else:
            print('All samples submitted. Submission loop ending.')
    except Exception:
        excpt_str = traceback.format_exc()
        print(excpt_str)
        with open('./PIPELINE_SUBMISSION_ERROR', 'w') as out:
            out.write(excpt_str + '\n')
        raise