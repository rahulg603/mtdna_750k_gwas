#!/usr/bin/env python3
import hail as hl
import os, re


######### CONSTANTS ##########
PHENO_KEY_FIELDS = ('trait_type', 'phenocode', 'pheno_sex', 'modifier')
MIN_CASES = 50
CHROMOSOMES = list(map(str, range(1, 23))) + ['X', 'XY']
SEXES = ('both_sexes', 'females', 'males')
POPS = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']
BASE_NONPC_COVARS = ['sex','age','age2','age_sex','age2_sex','site_bcm','site_uw']

PHENO_DESCRIPTION_FIELDS = ('description', 'description_more', 'coding_description', 'category')
PHENO_COLUMN_FIELDS = ('n_cases_both_sexes', 'n_cases_females', 'n_cases_males', *PHENO_DESCRIPTION_FIELDS)


######### PATHING ##########
def get_custom_ukb_pheno_mt_path(pheno_folder, suffix):
    return os.path.join(pheno_folder, f'/mt/phenotype_{suffix}.mt')


def get_custom_phenotype_summary_backup_path(pheno_folder, suffix, curdate):
    return os.path.join(pheno_folder, f'/summary/all_pheno_summary_{suffix}_before_{curdate}.txt.bgz')


def get_custom_phenotype_summary_path(pheno_folder, suffix, extension = 'ht'):
    return os.path.join(pheno_folder, f'/summary/phenotype_{suffix}.{extension}')


def get_custom_munged_pheno_path(pheno_folder, suffix):
    return os.path.join(pheno_folder, f'/mt/munged/munged_raw_phenotype_{suffix}.mt')


def get_pheno_export_dir(pheno_folder, suffix, pop):
    return os.path.join(pheno_folder, f'/exported/{suffix}/{pop}')


def format_pheno_dir(pheno):
    return pheno.replace("/", "_")


def pheno_dict_to_str(dct, for_save=False):
    return '-'.join([format_pheno_dir(dct[x])
                     if x == 'phenocode' and for_save
                     else dct[x] for x in PHENO_KEY_FIELDS if x in dct])


def pheno_str_to_dict(str):
    split_str = str.split('-')
    if len(split_str) != len(PHENO_KEY_FIELDS):
        raise ValueError('Input phenotype to pheno_str_to_dict does not have the expected number of fields.')
    return dict(zip(PHENO_KEY_FIELDS, split_str))


def get_pheno_output_path(pheno_export_dir, pheno_coding_trait, extension = '.tsv'):
    extended_suffix = f'{pheno_coding_trait["pheno_sex"]}-{pheno_coding_trait["modifier"]}'
    return os.path.join(pheno_export_dir, f'{pheno_coding_trait["trait_type"]}-{format_pheno_dir(pheno_coding_trait["phenocode"])}-{extended_suffix}{extension}')


def get_base_covariates_path(cov_folder):
    return os.path.join(cov_folder, f'/base/ht/baseline_covariates.ht')


def get_null_model_path(gs_output_path, suffix, pop):
    return os.path.join(gs_output_path, f'/null_glmm/{suffix}/{pop}')


def get_null_model_file_paths(null_model_path, pheno_dict, analysis_type):
    root = get_pheno_output_path(null_model_path, pheno_dict, '')
    return f'{root}.rda', f'{root}.{analysis_type}.varianceRatio.txt'


def get_result_path(gs_output_path, suffix, pop):
    return os.path.join(gs_output_path, f'/result/{suffix}/{pop}')


def get_results_prefix(pheno_results_dir, pheno_key_dict, chromosome):
    prefix = os.path.join(pheno_results_dir, f'/result_') + pheno_dict_to_str(pheno_key_dict, True)
    return f'{prefix}_chr{chromosome}'


def get_results_files(results_pre, analysis_type):
    if analysis_type == 'variant':
        return f'{results_pre}.single_variant.txt', None
    else:
        return f'{results_pre}.txt', f'{results_pre}.txt.singleAssoc.txt'


def get_merged_ht_path(gs_output_path, suffix, pop, pheno_dct):
    result_dir = get_result_path(gs_output_path, suffix, pop)
    return f'{get_pheno_output_path(result_dir, pheno_dct, "")}/variant_results.ht'

######### IMPORT UTILS ##########
def parse_bucket(gs_bucket):
    if re.search('^gs://', gs_bucket):
        return gs_bucket
    else:
        return f'gs://{gs_bucket}'


def get_covariates_with_custom(cov_folder, custom = None):
    new_covariates = hl.read_table(get_base_covariates_path(cov_folder))
    if custom is None:
        return new_covariates, []
    else:
        custom_covariates = hl.import_table(custom, impute=True)
        custom_covariates = custom_covariates.key_by(s = hl.str(custom_covariates.s))
        cust_covar_list = [x for x in custom_covariates.row if x not in custom_covariates.key]
        return new_covariates.annotate(**custom_covariates[new_covariates.key]), cust_covar_list


def pheno_ht_to_mt(pheno_ht: hl.Table, data_type: str, special_fields: str = ('age', 'sex')):
    """
    Input Hail Table with lots of phenotype row fields, distill into
    MatrixTable with either categorical or continuous data types
    as entries

    :param Table pheno_ht: Input hail Table with phenotypes as row fields
    :param str data_type: one of "categorical" or "continuous"
    :return: Hail MatrixTable with phenotypes as entries
    :rtype: MatrixTable
    """
    if data_type == 'categorical':
        filter_type = {hl.tbool}
        value_type = hl.bool
    else:
        filter_type = {hl.tint, hl.tfloat}
        value_type = hl.float

    special_fields_to_include = []
    fields = set(pheno_ht.row_value)
    for field in special_fields:
        if field in fields:
            fields.remove(field)
            special_fields_to_include.append(field)
    select_fields = {x: value_type(pheno_ht[x]) for x in fields if pheno_ht[x].dtype in filter_type}
    pheno_ht = pheno_ht.select(*special_fields_to_include, **select_fields)

    mt = pheno_ht.to_matrix_table_row_major(
        columns=list(select_fields), entry_field_name='value', col_field_name='phesant_pheno'
    )
    return mt


def format_entries(field, sex_field):
    return dict(
        both_sexes=hl.float64(field),
        females=hl.float64(hl.or_missing(sex_field == 0, field)),
        males=hl.float64(hl.or_missing(sex_field == 1, field))
    )


def load_custom_pheno_with_covariates(data_path, trait_type, modifier, 
                                      cov_folder, custom = None,
                                      sex: str = 'both_sexes', sample_col='s'):
    print(f'Loading {data_path}...')
    extension = os.path.splitext(data_path)[1]
    if extension == 'ht':
        ht = hl.read_table(data_path)
    else:
        if extension == 'tsv.gz':
            ht = hl.import_table(data_path, impute=True, force=True)
        else:
            ht = hl.import_table(data_path, impute=True)
        ht = ht.key_by(userId=ht[sample_col])
        if sample_col != 'userId':
            ht = ht.drop(sample_col)
        if trait_type == 'categorical':
            ht = ht.annotate(**{x: hl.bool(ht[x]) for x in list(ht.row_value)})

    mt = pheno_ht_to_mt(ht, trait_type).annotate_cols(data_type=trait_type)
    mt = mt.key_cols_by(trait_type=trait_type, phenocode=mt.phesant_pheno, pheno_sex=sex, modifier=modifier).drop('phesant_pheno')

    print(f'Now loading covariate table...')
    cov_ht, cust_covar_list = get_covariates_with_custom(cov_folder, custom)
    cov_ht = cov_ht.persist()
    
    mt_this = mt.select_rows(**cov_ht[mt.row_key])
    mt_this = mt_this.select_entries(**format_entries(mt_this.value, mt_this.sex))
    mt_this = mt_this.select_cols(**{f'n_cases_{sex}': hl.agg.count_where(
            hl.cond(mt_this.trait_type == 'categorical', mt_this[sex] == 1.0, hl.is_defined(mt_this[sex]))
         ) for sex in SEXES})
    full_mt = mt_this.unfilter_entries()
    
    return full_mt, cust_covar_list


def summarize_data(pheno_folder, suffix, overwrite):
    mt = hl.read_matrix_table(get_custom_ukb_pheno_mt_path(pheno_folder, suffix))
    ht = mt.group_rows_by('pop').aggregate(
        stats=hl.agg.stats(mt.both_sexes),
        n_cases_by_pop=hl.if_else(hl.set({'continuous', 'biomarkers'}).contains(mt.trait_type),
                                  hl.agg.count_where(hl.is_defined(mt.both_sexes)),
                                  hl.int64(hl.agg.sum(mt.both_sexes)))
    ).entries()
    ht = ht.key_by('pop', *PHENO_KEY_FIELDS)
    ht = ht.checkpoint(get_custom_phenotype_summary_path(pheno_folder, suffix), overwrite=overwrite, _read_if_exists=not overwrite)
    ht.flatten().export(get_custom_phenotype_summary_path(pheno_folder, suffix, 'tsv'))