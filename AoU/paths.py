#!/usr/bin/env python3
import hail as hl
import os

# AoU DRC paths
#ARRAY_PATH = os.getenv('MICROARRAY_HAIL_STORAGE_PATH')
#EXOME_PATH = os.getenv('WGS_EXOME_SPLIT_HAIL_PATH')
#WGS_PATH = os.getenv('WGS_ACAF_THRESHOLD_SPLIT_HAIL_PATH')
ARRAY_PATH = "gs://fc-aou-datasets-controlled/v7/microarray/hail.mt_v7.1"
EXOME_PATH = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/exome_v7.1/splitMT/hail.mt"
WGS_PATH = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/splitMT/hail.mt"
AUX_PATH = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux'
DATASET = os.getenv("WORKSPACE_CDR")

# Workspace paths
#BUCKET = os.getenv('WORKSPACE_BUCKET')
BUCKET = 'gs://fc-secure-34d99fb0-3748-4367-8197-b01069a4a7f9/'
REFERENCE_PATH = os.path.join(BUCKET, 'reference_files')
HAIL_GWAS_PATH = os.path.join(BUCKET, 'hail_gwas')
GWAS_PATH = os.path.join(BUCKET, 'saige_gwas')
GENO_PATH = os.path.join(GWAS_PATH, 'genotypes')
PHENO_PATH = os.path.join(GWAS_PATH, 'phenotypes')
COVAR_PATH = os.path.join(GWAS_PATH, 'covariates')
RESULTS_PATH = os.path.join(GWAS_PATH, 'results')
TEMP_PATH = os.path.join(BUCKET, 'tmp')

# Call rate
CALLRATE_CUTOFF = 0.9


def get_path_vat():
    return os.path.join(AUX_PATH, 'vat/vat_complete_v7.1.bgz.tsv.gz')

def get_ancestry_flat_file():
    return os.path.join(AUX_PATH, 'ancestry/ancestry_preds.tsv')

def get_genomic_metric_flat_file():
    return os.path.join(AUX_PATH, 'qc/genomic_metrics.tsv')

def get_aou_flagged_samples():
    return os.path.join(AUX_PATH, 'qc/flagged_samples.tsv')

def get_aou_hq_sites():
    return os.path.join(AUX_PATH, 'ancestry/merged_sites_only_intersection.vcf.bgz')

def get_aou_related_samples():
    return os.path.join(AUX_PATH, 'relatedness/relatedness_flagged_samples.tsv')


def get_all_by_all_qc_samples():
    return os.path.join(REFERENCE_PATH, '241021_allbyall_qc_sample_ids.tsv')
