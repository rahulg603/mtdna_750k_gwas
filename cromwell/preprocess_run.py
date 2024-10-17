#!/usr/bin/env python3
import json
from copy import deepcopy
import pandas as pd
import numpy as np
from cromwell.constants import *


def verify_json_template(j, prefix_delimiter='.'):
    keys_split = [k.split(prefix_delimiter) for k, _ in j.items()]
    
    if len([x for x in keys_split if len(x) < 2]) > 0:
        raise ValueError(f'Elements of the JSON field names do not follow appropriate convention of workflow_name.input_name.')
    
    workflow_names = set([x[0] for x in keys_split])

    if len(workflow_names) != 1:
        raise ValueError('Elements of the template JSON must contain exactly one workflow name. More than one was found.')
    
    print(f'Processing inputs under name {list(workflow_names)[0]}...')
    return j, list(workflow_names)[0]


def initialize_with_status(df, restart):
    if restart or (np.sum(df.columns == 'cromwell_id') == 0):
        df = df.assign(**{'cromwell_id': 'Submission pending'})
    if restart or (np.sum(df.columns == 'inputs_uri') == 0):
        df = df.assign(**{'inputs_uri': 'Submission pending'})
    if restart or (np.sum(df.columns == 'run_status') == 0):
        df = df.assign(**{'run_status': 'Submission pending'})
    
    all_status = set(df['run_status'])
    if not all_status.issubset(ALLOWED_STATUS):
        raise ValueError('ERROR: malformed status was input into CromwellManager.')

    return df


def compute_batches(df_in, batch_size):
    """ Will compute batch IDs within a dataframe. If any samples are completed they will 
    retain their batch (if present) or get a missing value.
    """
    df = deepcopy(df_in)
    if batch_size is None:
        print(f'Batches will not be used.')
        batch_size = 1
    else:
        print(f'Computing new batches of size {batch_size}. Batches will not include any non-pending samples.')
    
    df_not_submit = df.loc[df['run_status'] != 'Submission pending']
    if df_not_submit.shape[0] > 0:
        if np.sum(df_not_submit.columns == 'batch') == 0:
            raise ValueError('ERROR: there must be a batch column initialized if there are already submitted samples.')
        batch_offset = max(df_not_submit.batch)
    else:
        batch_offset = 0
    
    df_submit = df.loc[df['run_status'] == 'Submission pending']
    df_submit['batch'] = 0
    row_idx = 0
    batch_idx = batch_offset + 1
    ele_in_batch = 1
    while row_idx < df_submit.shape[0]:
        df_submit.at[df_submit.index[row_idx],'batch'] = batch_idx
        row_idx += 1
        if ele_in_batch == batch_size:
            ele_in_batch = 1
            batch_idx += 1
        else:
            ele_in_batch += 1

    return pd.concat([df_not_submit, df_submit], axis=0)
    