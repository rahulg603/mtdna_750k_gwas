#!/usr/bin/env python3
import os

BUCKET = os.getenv("WORKSPACE_BUCKET")
PROJECT = os.getenv("GOOGLE_PROJECT")

N_PARALLEL_WORKFLOWS = 1000

MIN_FAIL_NUM = 5 # Allow this many failures before worrying about the percentage failure
MAX_FAIL_PCT = 15.0 # Halt the submission loop if more than this percentage of batches fail

ALLOWED_STATUS = ('Submission pending', 'Submitted', 'Running', 'Succeeded', 'Failed', 'Manual review')
COMPLETED_STATUS = ('Succeeded', 'Failed', 'Manual review')
STATUS_COLS = ['cromwell_id', 'inputs_uri', 'run_status', 'batch']