#!/usr/bin/env python3
import hail as hl
import os, re


######### CONSTANTS ##########
PHENO_KEY_FIELDS = ('trait_type', 'phenocode', 'pheno_sex', 'modifier')
MIN_CASES = 50
CHROMOSOMES = ['chr' + x for x in list(map(str, range(1, 23))) + ['X', 'Y']]
AUTOSOMES = ['chr' + x for x in list(map(str, range(1, 23)))]
SEXES = ('both_sexes', 'females', 'males')
POPS = ['mid', 'eas', 'sas', 'amr', 'afr', 'eur']
BASE_NONPC_COVARS = ['sex','age','age2','age_sex','age2_sex','site_id_bcm','site_id_uw']

PHENO_SAMPLE_ID = 'userId'

SHORT_READ_ROOT = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/"

PHENO_DESCRIPTION_FIELDS = ('description', 'description_more', 'coding_description', 'category')
PHENO_COLUMN_FIELDS = ('n_cases_both_sexes', 'n_cases_females', 'n_cases_males', *PHENO_DESCRIPTION_FIELDS)

SAIGE_PHENO_TYPES = {
    'continuous': 'quantitative',
    'biomarkers': 'quantitative',
    'categorical': 'binary',
    'icd': 'binary',
    'icd_first_occurrence': 'binary',
    'icd_all': 'binary',
    'phecode': 'binary',
    'prescriptions': 'binary'
}

HLA_LOCUS = "chr6:28510120-33480577"
INVERSION_LOCUS = "chr8:8198267-12123140"

######### PATHING AND MUNGING ##########
def get_aou_util_path(util):
    raise NotImplementedError('get_aou_util_path not implemented.')


# Genotypes
def get_plink_for_null_path(geno_folder, pop, sample_qc, analysis_type, ld_pruned,
                            n_common, n_maf, n_mac,
                            use_drc_ancestry_data=False, use_array_for_variant=False):
    # recall, we always produce this with array data
    path_finder = lambda x: get_sites_for_null_path(geno_folder=geno_folder,
                                                    pop=pop, sample_qc=sample_qc, analysis_type=analysis_type,
                                                    ld_pruned=ld_pruned, n_common=n_common,
                                                    n_mac=n_mac, n_maf=n_maf,
                                                    use_drc_ancestry_data=use_drc_ancestry_data,
                                                    use_array_for_variant=use_array_for_variant,
                                                    extension=x)
    return path_finder('bed'), path_finder('bim'), path_finder('fam')


def get_sparse_grm_path(geno_folder, pop, n_markers, relatedness, sample_qc, use_plink, use_drc_ancestry_data=False, af_cutoff=0.01):
    af = f'_maf{str(af_cutoff)}'
    related = f'_rel{str(relatedness)}'
    drc_string = '_drc' if use_drc_ancestry_data else '_axaou'
    qc = '_sample_qc' if sample_qc else ''
    plink = '_plink' if use_plink else ''
    prefix = os.path.join(geno_folder, f'sparse_grm/aou_ld_pruned_{pop}{qc}{drc_string}{af}{plink}_{str(n_markers)}markers{related}')
    return f'{prefix}.mtx', f'{prefix}.mtx.sampleIDs.txt'


def get_call_stats_ht_path(geno_folder, pop, sample_qc, analysis_type, use_drc_ancestry_data=False, use_array_for_variant=False):
    if analysis_type == 'variant':
        source_str = '_array' if use_array_for_variant else '_wgs'
    else:
        source_str = '_exome'
    prune_str = '_sample_qc' if sample_qc else ''
    drc_string = '_drc' if use_drc_ancestry_data else '_axaou'
    return os.path.join(geno_folder, f'call_stats/call_stats{source_str}_{pop}{prune_str}{drc_string}.ht')


def get_sites_for_null_path(geno_folder, pop, sample_qc, analysis_type, ld_pruned,
                            n_common, n_maf, n_mac, extension,
                            use_drc_ancestry_data=False, use_array_for_variant=False):
    if analysis_type == 'variant':
        source_str = '_array' if use_array_for_variant else '_wgs'
    elif analysis_type == 'gene':
        source_str = '_exome'
    elif analysis_type == 'both':
        source_str = '_array_and_exome' if use_array_for_variant else '_wgs_and_exome'
    prune_str = '_ldpruned' if ld_pruned else ''
    drc_string = '_drc' if use_drc_ancestry_data else '_axaou'
    qc = '_sample_qc' if sample_qc else ''
    varct = f'_N{str(n_common)}_{str(n_maf)}_{str(n_mac)}'
    return os.path.join(geno_folder, f'subsampled/sites_for_grm{source_str}{prune_str}_{pop}{qc}{drc_string}{varct}.{extension}')


def get_ld_pruned_array_data_path(geno_folder, pop, extension, sample_qc, use_plink, use_drc_ancestry_data=False, window='1e7', af_cutoff=0.01):
    drc_string = '_drc' if use_drc_ancestry_data else '_axaou'
    qc = '_sample_qc' if sample_qc else ''
    af = f'_af{str(af_cutoff)}'
    plink = '_plink' if use_plink else ''
    return os.path.join(geno_folder, f'ld_prune/ld_pruned_sites_array_{pop}{qc}{drc_string}{af}{plink}_window{window}.{extension}')


def get_plink_inputs_ld_prune(geno_folder, pop, chr, extension, sample_qc, pruned=None, use_drc_ancestry_data=False, af_cutoff=0.01):
    drc_string = '_drc' if use_drc_ancestry_data else '_axaou'
    qc = '_sample_qc' if sample_qc else ''
    af = f'_af{str(af_cutoff)}'
    prun = f'_pruned{pruned}' if pruned is not None else ''
    return os.path.join(geno_folder, f'ld_prune/plink_chr/variants_for_ld_pruning_array_{str(chr)}_{pop}{qc}{drc_string}{af}{prun}.{extension}')


def get_wildcard_path_genotype_bgen(analysis_type):
    data_type = 'exome_v7.1' if analysis_type == 'gene' else 'acaf_threshold_v7.1'
    file = 'exome' if analysis_type == 'gene' else 'acaf_threshold'
    return os.path.join(SHORT_READ_ROOT,f'{data_type}/bgen/{file}@')


# Samples
def get_n_samples_per_pop_path(geno_folder, sample_qc, analysis_type, use_drc_ancestry_data=False, use_array_for_variant=False):
    if analysis_type == 'variant':
        source_str = '_array' if use_array_for_variant else '_wgs'
    else:
        source_str = '_exome'
    prune_str = '_sample_qc' if sample_qc else ''
    drc_string = '_drc' if use_drc_ancestry_data else '_axaou'
    return os.path.join(geno_folder, f'sample_counts/counts_by_pop{source_str}{prune_str}{drc_string}.tsv') 


def get_aou_samples_file_path(geno_folder, pop, sample_qc, use_plink=True, use_drc_ancestry_data=False):
    qc = '_sample_qc' if sample_qc else ''
    drc_string = '_drc' if use_drc_ancestry_data else '_axaou'
    plink = '_plink' if use_plink else ''
    return os.path.join(geno_folder, f'sample_counts/sample_list_{pop}{qc}{drc_string}{plink}.tsv') 


# Phenotypes
def get_custom_ukb_pheno_mt_path(pheno_folder, suffix):
    return os.path.join(pheno_folder, f'mt/phenotype_{suffix}.mt')


def get_custom_phenotype_summary_backup_path(pheno_folder, suffix, curdate):
    return os.path.join(pheno_folder, f'summary/all_pheno_summary_{suffix}_before_{curdate}.txt.bgz')


def get_custom_phenotype_summary_path(pheno_folder, suffix, extension = 'ht'):
    return os.path.join(pheno_folder, f'summary/phenotype_{suffix}.{extension}')


def get_custom_munged_pheno_path(pheno_folder, suffix):
    return os.path.join(pheno_folder, f'mt/munged/munged_raw_phenotype_{suffix}.mt')


def get_pheno_export_dir(pheno_folder, suffix, pop):
    return os.path.join(pheno_folder, f'exported/{suffix}/{pop}')


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


# Covariates
def get_base_covariates_path(cov_folder, drc):
    drc_string = '_drc' if drc else '_axaou'
    return os.path.join(cov_folder, f'base/ht/baseline_covariates{drc_string}.ht')


def get_demographics_path(cov_folder, drc):
    drc_string = '_drc' if drc else '_axaou'
    return os.path.join(cov_folder, f'base/ht/all_covariates{drc_string}.ht')


# Null model
def get_null_model_path(gs_output_path, suffix, pop):
    return os.path.join(gs_output_path, f'null_glmm/{suffix}/{pop}')


def get_null_model_file_paths(null_model_path, pheno_dict, analysis_type):
    root = get_pheno_output_path(null_model_path, pheno_dict, '')
    return f'{root}.rda', f'{root}.{analysis_type}.varianceRatio.txt'


# GWAS results
def get_result_path(gs_output_path, suffix, pop):
    return os.path.join(gs_output_path, f'result/{suffix}/{pop}')


def get_results_prefix(pheno_results_dir, pheno_key_dict, chromosome):
    prefix = os.path.join(pheno_results_dir, f'result_') + pheno_dict_to_str(pheno_key_dict, True)
    return f'{prefix}_{chromosome}'


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


def get_covariates_with_custom(cov_folder, custom=None, drc=False):
    new_covariates = hl.read_table(get_base_covariates_path(cov_folder, drc=drc))
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
                                      cov_folder, custom=None, drc=False,
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
        ht = ht.key_by(**{PHENO_SAMPLE_ID: hl.str(ht[sample_col])})
        if sample_col != PHENO_SAMPLE_ID:
            ht = ht.drop(sample_col)
        if trait_type == 'categorical':
            ht = ht.annotate(**{x: hl.bool(ht[x]) for x in list(ht.row_value)})

    mt = pheno_ht_to_mt(ht, trait_type).annotate_cols(data_type=trait_type)
    mt = mt.key_cols_by(trait_type=trait_type, phenocode=mt.phesant_pheno, pheno_sex=sex, modifier=modifier).drop('phesant_pheno')
    mt.describe()

    print(f'Now loading covariate table...')
    cov_ht, cust_covar_list = get_covariates_with_custom(cov_folder, custom, drc=drc)
    cov_ht = cov_ht.persist()
    
    mt_this = mt.select_rows(**cov_ht[mt.row_key])
    mt_this = mt_this.annotate_rows(**{k: hl.int(v)  for k, v in mt_this.row.items() if v.dtype == hl.tbool})
    mt_this = mt_this.select_entries(**format_entries(mt_this.value, mt_this.sex))
    mt_this = mt_this.select_cols(**{f'n_cases_{sex}': hl.agg.count_where(
            hl.cond(mt_this.trait_type == 'categorical', mt_this[sex] == 1.0, hl.is_defined(mt_this[sex]))
         ) for sex in SEXES})
    full_mt = mt_this.unfilter_entries()
    
    return full_mt, cust_covar_list


def get_custom_ukb_pheno_mt(pheno_folder, cov_folder, custom_covariates, suffix, pop: str = 'all', drc = False):
    mt = hl.read_matrix_table(get_custom_ukb_pheno_mt_path(pheno_folder, suffix))
    covars, _ = get_covariates_with_custom(cov_folder=cov_folder, custom=custom_covariates, drc=drc)
    mt = mt.annotate_rows(**covars[mt.row_key])
    if pop != 'all':
        mt = mt.filter_rows(mt.pop == pop)
    
    return mt


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