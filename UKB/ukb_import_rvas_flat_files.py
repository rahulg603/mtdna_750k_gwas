import hail as hl
import os, re
import secrets
from tqdm import tqdm
from google.cloud import storage
from urllib.parse import urlparse

from utils.SaigeImporters import mwzj_hts_by_tree, GENE_COL_KEY_FIELDS, GENE_ROW_KEY_FIELDS
from AoU.sumstats import unify_saige_gene_ht_schema, split_initial_saige_gene_ht


PHENO_KEY_FIELDS = ['trait_type','phenocode','pheno_sex','coding','modifier']
TEMP_DIR = 'gs://temp-7day/'


def read_log_file(path):
    """Reads log file content from local or GCS."""
    if path.startswith("gs://"):
        # Parse GCS path
        parsed = urlparse(path)
        bucket_name = parsed.netloc
        blob_path = parsed.path.lstrip("/")

        # Download content from GCS
        client = storage.Client()
        bucket = client.bucket(bucket_name)
        blob = bucket.blob(blob_path)
        content = blob.download_as_text().splitlines()
    else:
        # Local file
        with open(path, "r") as f:
            content = f.readlines()
    return content


def dx_get_cases_controls_from_logs(log_list):
    """
    Different from the utils function because we have to strip the first part of the string.
    """
    cases = controls = -1
    for log in log_list:
        try:
            lines = read_log_file(log)
            for line in lines:
                line = re.sub(r'^\[.*?\]\s', '', line)
                line = line.strip()
                if line.startswith('Analyzing'):
                    fields = line.split()
                    if len(fields) == 6:
                        try:
                            cases = int(fields[1])
                            controls = int(fields[4])
                            break
                        except ValueError:
                            print(f'Could not load number of cases or controls from {line}.')
                elif line.endswith('samples were used in fitting the NULL glmm model and are found in sample file') or \
                        line.endswith('samples have been used to fit the glmm null model') or \
                        line.endswith('samples will be used for analysis'):
                    # This is ahead of the case/control count line ("Analyzing ...") above so this should be ok
                    fields = line.split()
                    try:
                        cases = int(fields[0])
                    except ValueError:
                        print(f'Could not load number of cases or controls from {line}.')
            return cases, controls
        except:
            pass
    return cases, controls


def get_this_n(this_flat):
    noext_path = os.path.splitext(os.path.splitext(this_flat)[0])[0]
    null_log = noext_path + '.variant.log'
    if hl.hadoop_exists(null_log):
        cases, _ = dx_get_cases_controls_from_logs([null_log])
    else:
        print('WARNING: did not find log for ' + os.path.basename(noext_path) + '.')
        this_base = os.path.join(os.path.dirname(noext_path), 'phenotypes_munged', os.path.basename(noext_path).replace('_merged','') + '.tsv')
        data = hl.import_table(this_base, types={'value':hl.tfloat64})
        cases = data.filter(hl.is_defined(data.value)).count()
    
    print('For ' + os.path.basename(noext_path) + ', got N=' + str(cases) + '.')
    return cases


def generate_dict(this_name):
    split_name = this_name.split('-')
    if len(split_name) != 5:
        raise ValueError(f'ERROR: phenotype name {this_name} does not conform to expected specification.')
    if not re.search('age_accum', split_name[1]):
        if re.search('_nochip$', split_name[1]):
            split_name[1] = split_name[1].replace('_nochip','') + '_incl_nonaccum_nochip'
        else:
            split_name[1] = split_name[1] + '_incl_nonaccum'
    if re.search('invHL', split_name[1]) or re.search('1minHL', split_name[1]):
        split_name[1] = split_name[1].replace('_count','')
    
    return {k: v for k, v in zip(PHENO_KEY_FIELDS, split_name)}


def saige_merge_raw_ukb_rvas_flat(path, read_previous=False, overwrite=True, n_partitions=200):
    inner_mode = 'overwrite' if overwrite else '_read_if_exists'
    merged_mt_path = os.path.join(path, 'results_EUR.mt')

    if read_previous and hl.hadoop_exists(f'{merged_mt_path}/_SUCCESS'):
        return None

    all_flat_files = [x['path'] for x in hl.hadoop_ls(path) if not x['is_dir'] and re.search('_merged.tsv.gz$', x['path'])]
    all_flat_basename = [os.path.basename(os.path.splitext(os.path.splitext(x)[0])[0]) for x in all_flat_files]
    all_dicts = [generate_dict(x) for x in all_flat_basename]
    all_ht = [os.path.join(path, x) + '.ht' for x in all_flat_basename]

    for this_flat, (this_ht, this_dict) in zip(all_flat_files, zip(all_ht, all_dicts)):
        if overwrite or not hl.hadoop_exists(f'{this_ht}/_SUCCESS'):
            ht = hl.import_table(this_flat, force=True, 
                                 types={'max_MAF':hl.tfloat64, 
                                        'Pvalue': hl.tfloat, 'Pvalue_Burden': hl.tfloat, 'Pvalue_SKAT': hl.tfloat,
                                        'BETA_Burden': hl.tfloat, 'SE_Burden': hl.tfloat, 'MAC': hl.tint,
                                        'Number_rare': hl.tint, 'Number_ultra_rare': hl.tint})
            ht = ht.transmute(gene_symbol = ht.Region,
                              group = ht.Group,
                              **{x: this_dict[x] for x in PHENO_KEY_FIELDS})
            ht = ht.annotate_globals(n_cases = get_this_n(this_flat),
                                     n_controls = hl.missing(hl.tint32),
                                     heritability = hl.missing(hl.tfloat64),
                                     saige_version = hl.missing(hl.tstr),
                                     inv_normalized = "True")
            ht = ht.key_by('gene_symbol', 'group', 'max_MAF', *PHENO_KEY_FIELDS)
            ht.repartition(20).write(this_ht, overwrite=overwrite)

    internal_temp_dir = f'{TEMP_DIR}/EUR/gene/{secrets.token_urlsafe(12)}/'
    if len(all_ht) > 0:
        row_keys = ['gene_symbol']
        initial_col_keys = PHENO_KEY_FIELDS
        initial_hts = [unify_saige_gene_ht_schema(hl.read_table(x), x, internal_temp_dir, row_keys, initial_col_keys) for x in tqdm(all_ht)]
        new_col_keys = GENE_COL_KEY_FIELDS
        row_keys = GENE_ROW_KEY_FIELDS
        col_keys = initial_col_keys + list(new_col_keys)
        all_hts = [y for x in initial_hts for y in split_initial_saige_gene_ht(x, row_keys, new_col_keys)]

        print('Schemas unified. Starting joining...')

        mt = mwzj_hts_by_tree(all_hts, internal_temp_dir, col_keys, debug=True, gene_analysis=True,
                              inner_mode=inner_mode, repartition_final=n_partitions)
        
        print(f'Unioned MTs...')
        print('After merge schema...')
        mt.describe()
        mt = mt.checkpoint(f'{internal_temp_dir}/staging.mt', **{inner_mode: True})

        # patch keys
        key_set = PHENO_KEY_FIELDS + list(GENE_COL_KEY_FIELDS)
        mt = mt.key_cols_by(**{x: hl.case(missing_false=True)
                            .default(mt[x])
                            for x in key_set})

        if mt.inv_normalized.dtype == hl.tstr:
            mt = mt.annotate_cols(inv_normalized=hl.bool(mt.inv_normalized))

        mt = mt.filter_cols(mt.phenocode != "")
        mt = mt.key_rows_by(*row_keys)

        print('Prior to output schema...')
        mt.describe()

        mt.write(merged_mt_path, overwrite=True)
    
    return None