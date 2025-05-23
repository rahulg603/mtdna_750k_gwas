#!/usr/bin/env python3
import hail as hl
import os, re
from datetime import date
from copy import deepcopy


######### CONSTANTS ##########
PHENO_KEY_FIELDS = ('trait_type', 'phenocode', 'pheno_sex', 'modifier')
GENE_COL_KEY_FIELDS = ('group', 'max_MAF')
GENE_ROW_KEY_FIELDS = ['gene_symbol']
MIN_CASES = 50
CHROMOSOMES = ['chr' + x for x in list(map(str, range(1, 23))) + ['X', 'Y']]
AUTOSOMES = ['chr' + x for x in list(map(str, range(1, 23)))]
SEXES = ('both_sexes', 'females', 'males')
POPS = ['mid', 'eas', 'sas', 'amr', 'afr', 'eur']
BASE_NONPC_COVARS = ['sex','age','age2','age_sex','age2_sex','site_id_bcm','site_id_uw']
CHUNK_SIZE = {'mid': 2.5e7, 'sas': 2.5e7, 'eas': 2.5e7, 'amr': 1e7, 'afr': 1e7, 'eur': 1e7}

PHENO_SAMPLE_ID = 'userId'

SHORT_READ_ROOT = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/"

PHENO_DESCRIPTION_FIELDS = ('description', 'description_more', 'coding_description', 'category')
PHENO_COLUMN_FIELDS = ('n_cases_both_sexes', 'n_cases_females', 'n_cases_males', *PHENO_DESCRIPTION_FIELDS)
PHENO_GWAS_FIELDS = ('n_cases', 'n_controls', 'heritability', 'saige_version', 'inv_normalized')

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

# fields specific to each category of trait
quant_meta_fields = ['AF_Allele2']
quant_fields = ['AF_Allele2']

binary_meta_fields = ['AF_Cases','AF_Controls']
binary_fields = ['AF.Cases','AF.Controls']

# dictionaries for renaming fields
quant_meta_field_rename_dict = {'AF_Allele2': 'af_meta',
                                'BETA': 'beta_meta',
                                'SE': 'se_meta',
                                'Pvalue': 'pval_meta',
                                'Pvalue_het': 'pval_heterogeneity',
                                'N': 'N_meta'}
quant_field_rename_dict = {'AF_Allele2': 'af',
                           'BETA': 'beta',
                           'SE': 'se',
                           'Pvalue': 'pval',
                           'low_confidence': 'low_confidence'}

binary_meta_field_rename_dict = {'BETA': 'beta_meta',
                                 'SE': 'se_meta',
                                 'Pvalue': 'pval_meta',
                                 'AF_Cases': 'af_cases_meta',
                                 'AF_Controls': 'af_controls_meta',
                                 'Pvalue_het': 'pval_heterogeneity',
                                 'N': 'N_meta'}
binary_field_rename_dict = {'AF.Cases': 'af_cases',
                            'AF.Controls': 'af_controls',
                            'BETA': 'beta',
                            'SE': 'se',
                            'Pvalue': 'pval',
                            'low_confidence': 'low_confidence'}

# consequences
CSQ_CODING_HIGH_IMPACT = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
]
CSQ_CODING_MEDIUM_IMPACT = [
    "start_lost",  # considered high impact in v105, previously medium
    "initiator_codon_variant",  # deprecated
    "transcript_amplification",  # considered high impact in v105, previously medium
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",  # new in v79
]
CSQ_CODING_LOW_IMPACT = [
    "splice_donor_5th_base_variant",  # new in v105
    "splice_region_variant",  # Considered low impact in v105, previously medium
    "splice_donor_region_variant",  # new in v105
    "splice_polypyrimidine_tract_variant",  # new in v105
    "incomplete_terminal_codon_variant",
    "start_retained_variant",  # new in v92
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "coding_transcript_variant",
]
CSQ_NON_CODING = [
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "non_coding_exon_variant",  # deprecated
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "nc_transcript_variant",  # deprecated
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "feature_elongation",
    "regulatory_region_variant",
    "feature_truncation",
    "intergenic_variant",
    "sequence_variant",
]
CSQ_ORDER = (
    CSQ_CODING_HIGH_IMPACT
    + CSQ_CODING_MEDIUM_IMPACT
    + CSQ_CODING_LOW_IMPACT
    + CSQ_NON_CODING
)
LOF_CSQ_SET = {
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
}

PLOF_CSQS = ["transcript_ablation", "splice_acceptor_variant",
             "splice_donor_variant", "stop_gained", "frameshift_variant"]

MISSENSE_CSQS = ["stop_lost", "start_lost", "transcript_amplification",
                 "inframe_insertion", "inframe_deletion", "missense_variant"]

SYNONYMOUS_CSQS = ["stop_retained_variant", "synonymous_variant"]

OTHER_CSQS = ["mature_miRNA_variant", "5_prime_UTR_variant",
              "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
              "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
              "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
              "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
              "regulatory_region_variant", "feature_truncation", "intergenic_variant"]


######### PATHING AND MUNGING ##########
def get_aou_util_path(util):
    raise NotImplementedError('get_aou_util_path not implemented.')


def update_suffix(suffix, use_drc_pop, use_custom_pcs):
    return suffix + ('_drcpop' if use_drc_pop else '') + ('_custompcs' if use_custom_pcs == 'custom' else ('_axaoupcs' if use_custom_pcs == 'axaou' else ''))


def rename_dict_for_log10(dct, legacy_exp_p_values=False):
    dct_out = deepcopy(dct)
    if not legacy_exp_p_values:
        for k, v in dct_out.items():
            if re.search('^Pvalue', k) and re.search('^pval', v):
                dct_out.update({k: re.sub('^pval', 'neglog10_pval',v)})
    return dct_out


def get_processed_vat_path(annot_folder):
    return os.path.join(annot_folder, 'vat_processed_final.ht')


def get_gene_annotation_path(annot_folder, chrom, pop, min_cr):
    if min_cr == 0:
        cr = '_noCRfilter'
    else:
        cr = ''
    return os.path.join(annot_folder, f'gene_annotations2/formatted_exome_variant_consequence_lof_mis_syn_{chrom}_{pop}{cr}.txt')


# Genotypes
def get_plink_for_null_path(geno_folder, pop, sample_qc, analysis_type, ld_pruned,
                            n_common, n_maf, n_mac,
                            use_drc_pop=False, use_array_for_variant=False):
    # recall, we always produce this with array data
    path_finder = lambda x: get_sites_for_null_path(geno_folder=geno_folder,
                                                    pop=pop, sample_qc=sample_qc, analysis_type=analysis_type,
                                                    ld_pruned=ld_pruned, n_common=n_common,
                                                    n_mac=n_mac, n_maf=n_maf,
                                                    use_drc_pop=use_drc_pop,
                                                    use_array_for_variant=use_array_for_variant,
                                                    extension=x)
    return path_finder('bed'), path_finder('bim'), path_finder('fam')


def get_sparse_grm_path(geno_folder, pop, n_markers, relatedness, sample_qc, use_plink, use_drc_pop=False, use_array_data='', af_cutoff=0.01):
    # use_array_data is expected to be A_B_C or blank, where A is n_common, B is n_maf, C is n_mac
    af = f'_maf{str(af_cutoff)}'
    related = f'_rel{str(relatedness)}'
    drc_string = '_drc' if use_drc_pop else '_axaou'
    qc = '_sample_qc' if sample_qc else ''
    plink = '_plink' if use_plink else ''
    use_array_data = '' if use_array_data == '' else f'_wgs_and_exome_N{use_array_data}'
    prefix = os.path.join(geno_folder, f'sparse_grm/aou_ld_pruned_{pop}{qc}{drc_string}{af}{plink}{use_array_data}_{str(n_markers)}markers{related}')
    return f'{prefix}.mtx', f'{prefix}.mtx.sampleIDs.txt'


def get_call_stats_ht_path(geno_folder, pop, sample_qc, analysis_type, use_drc_pop=False, use_array_for_variant=False):
    if analysis_type == 'variant':
        source_str = '_array' if use_array_for_variant else '_wgs'
    else:
        source_str = '_exome'
    prune_str = '_sample_qc' if sample_qc else ''
    drc_string = '_drc' if use_drc_pop else '_axaou'
    return os.path.join(geno_folder, f'call_stats/call_stats{source_str}_{pop}{prune_str}{drc_string}.ht')


def get_sites_for_null_path(geno_folder, pop, sample_qc, analysis_type, ld_pruned,
                            n_common, n_maf, n_mac, extension,
                            use_drc_pop=False, use_array_for_variant=False):
    if analysis_type == 'variant':
        source_str = '_array' if use_array_for_variant else '_wgs'
    elif analysis_type == 'gene':
        source_str = '_exome'
    elif analysis_type == 'both':
        source_str = '_array_and_exome' if use_array_for_variant else '_wgs_and_exome'
    prune_str = '_ldpruned' if ld_pruned else ''
    drc_string = '_drc' if use_drc_pop else '_axaou'
    qc = '_sample_qc' if sample_qc else ''
    varct = f'_N{str(n_common)}_{str(n_maf)}_{str(n_mac)}'
    return os.path.join(geno_folder, f'subsampled/sites_for_grm{source_str}{prune_str}_{pop}{qc}{drc_string}{varct}.{extension}')


def get_ld_pruned_array_data_path(geno_folder, pop, extension, sample_qc, use_plink, use_drc_pop=False, window='1e7', af_cutoff=0.01):
    drc_string = '_drc' if use_drc_pop else '_axaou'
    qc = '_sample_qc' if sample_qc else ''
    af = f'_af{str(af_cutoff)}'
    plink = '_plink' if use_plink else ''
    return os.path.join(geno_folder, f'ld_prune/ld_pruned_sites_array_{pop}{qc}{drc_string}{af}{plink}_window{window}.{extension}')


def get_plink_inputs_ld_prune(geno_folder, pop, chr, extension, sample_qc, pruned=None, use_drc_pop=False, af_cutoff=0.01):
    drc_string = '_drc' if use_drc_pop else '_axaou'
    qc = '_sample_qc' if sample_qc else ''
    af = f'_af{str(af_cutoff)}'
    prun = f'_pruned{pruned}' if pruned is not None else ''
    return os.path.join(geno_folder, f'ld_prune/plink_chr/variants_for_ld_pruning_array_{str(chr)}_{pop}{qc}{drc_string}{af}{prun}.{extension}')


def get_wildcard_path_genotype_bgen(analysis_type):
    data_type = 'exome_v7.1' if analysis_type == 'gene' else 'acaf_threshold_v7.1'
    file = 'exome' if analysis_type == 'gene' else 'acaf_threshold'
    return os.path.join(SHORT_READ_ROOT,f'{data_type}/bgen/{file}.@')


def get_wildcard_path_intervals_bgen(geno_folder, pop, use_drc_pop, encoding='additive', analysis_type='variant'):
    drc_string = '_drc' if use_drc_pop else '_axaou'
    fold = 'split_acaf_bgen' if analysis_type == 'variant' else 'split_exome_bgen'
    intermediate = 'wgs_acaf' if analysis_type == 'variant' else 'wes'
    if analysis_type == 'gene':
        if encoding != 'additive':
            raise ValueError('If doing a gene analysis, encoding must be additive.')
    return os.path.join(geno_folder, f'{fold}/aou_{intermediate}_qc{drc_string}_{pop}.@.#.?_{encoding}')


def get_saige_interval_path(geno_folder, pop, analysis_type):
    return os.path.join(geno_folder, f'intervals/aou_{analysis_type}_{pop}.tsv')


def read_variant_intervals(geno_folder, pop, analysis_type):
    variant_interval_path = get_saige_interval_path(geno_folder, analysis_type=analysis_type, pop=pop)
    df = hl.import_table(variant_interval_path, 
                         types={'chrom': hl.tstr, 'start': hl.tint, 'end': hl.tint}).to_pandas()
    print(pop)
    print(df.head(5))
    print(df.shape)
    return(df)


def stringify_interval(chr, start, end):
    return f'{chr}.{str(start)}.{str(end)}'


# Samples
def get_n_samples_per_pop_path(geno_folder, sample_qc, analysis_type, use_drc_pop=False, use_array_for_variant=False):
    if analysis_type == 'variant':
        source_str = '_array' if use_array_for_variant else '_wgs'
    else:
        source_str = '_exome'
    prune_str = '_sample_qc' if sample_qc else ''
    drc_string = '_drc' if use_drc_pop else '_axaou'
    return os.path.join(geno_folder, f'sample_counts/counts_by_pop{source_str}{prune_str}{drc_string}.tsv') 


def get_aou_samples_file_path(geno_folder, pop, sample_qc, use_plink=True, use_drc_pop=False):
    qc = '_sample_qc' if sample_qc else ''
    drc_string = '_drc' if use_drc_pop else '_axaou'
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


def get_pheno_output_suffix(pheno_coding_trait):
    extended_suffix = f'{pheno_coding_trait["pheno_sex"]}-{pheno_coding_trait["modifier"]}'
    return f'{pheno_coding_trait["trait_type"]}-{format_pheno_dir(pheno_coding_trait["phenocode"])}-{extended_suffix}'


def get_pheno_output_path(pheno_export_dir, pheno_coding_trait, extension = '.tsv'):
    return os.path.join(pheno_export_dir, f'{get_pheno_output_suffix(pheno_coding_trait)}{extension}')


# Covariates
def parse_drc_custom_inputs(use_drc_pop, use_custom_pcs):
    if use_custom_pcs not in ['custom', 'axaou', None]:
        raise NotImplementedError('ERROR: use_custom_pcs must be custom, axaou, or None (drc stock PCs).')
    if not use_drc_pop and (use_custom_pcs == 'custom'):
        raise NotImplementedError('ERROR: both custom PCs and axaou population definitions is not implemented.')

    if use_drc_pop:
        drc_string = '_drc'
    else:
        drc_string = '_axaou'
    
    if use_custom_pcs == 'custom':
        pcs = '_custompcs'
    elif use_custom_pcs == 'axaou':
        pcs = '_axaoupcs'
    elif use_custom_pcs is None:
        pcs = ''

    return drc_string, pcs


def get_base_covariates_path(cov_folder, use_drc_pop, use_custom_pcs):
    drc_string, pcs = parse_drc_custom_inputs(use_drc_pop, use_custom_pcs)
    path_out = os.path.join(cov_folder, f'base/ht/baseline_covariates{drc_string}{pcs}.ht')
    print(f'Accessing {path_out}...')
    return path_out


def get_demographics_path(cov_folder, use_drc_pop, use_custom_pcs):
    drc_string, pcs = parse_drc_custom_inputs(use_drc_pop, use_custom_pcs)
    path_out = os.path.join(cov_folder, f'base/ht/all_covariates{drc_string}{pcs}.ht')
    print(f'Accessing {path_out}...')
    return path_out


def get_custom_pc_path(cov_folder, iteration):
    return f'{cov_folder}/PCA/recomputed_pca_gnomad{make_iteration_suffix(iteration)}.ht'


# Null model
def get_null_model_path(gs_output_path, suffix, pop):
    return os.path.join(gs_output_path, f'null_glmm/{suffix}/{pop}')


def get_null_model_file_paths(null_model_path, pheno_dict, analysis_type):
    root = get_pheno_output_path(null_model_path, pheno_dict, '')
    return f'{root}.rda', f'{root}.{analysis_type}.varianceRatio.txt', f'{root}.null.log'


# GWAS results
def get_result_path(gs_output_path, suffix, pop, encoding):
    return os.path.join(gs_output_path, f'result/{suffix}/{encoding}/{pop}')


def get_results_prefix(pheno_results_dir, pheno_key_dict, chromosome):
    prefix = os.path.join(pheno_results_dir, f'result_') + pheno_dict_to_str(pheno_key_dict, True)
    return f'{prefix}_{chromosome}'


def get_results_files(results_pre, analysis_type):
    if analysis_type == 'variant':
        return f'{results_pre}.single_variant.txt', None, f'{results_pre}.log'
    else:
        return f'{results_pre}.singleAssoc.txt', f'{results_pre}.geneAssoc.txt', f'{results_pre}.log'


def get_merged_ht_path(gs_output_path, suffix, pop, pheno_dct, encoding, gene_analysis=False):
    result_dir = get_result_path(gs_output_path, suffix, pop, encoding)
    path = get_pheno_output_path(result_dir, pheno_dct, "")
    if gene_analysis:
        return f'{path}/gene_results_singleAssoc.ht', f'{path}/gene_results.ht'
    else:
        return f'{path}/variant_results.ht', ''


def get_merged_flat_path(gs_output_path, suffix, pop, pheno_dct, encoding, gene_analysis=False):
    result_dir = get_result_path(gs_output_path, suffix, pop, encoding)
    path = get_pheno_output_path(result_dir, pheno_dct, "")
    if gene_analysis:
        return f'{path}/gene_results_singleAssoc.tsv.bgz', f'{path}/gene_results.tsv.bgz'
    else:
        return f'{path}/variant_results.tsv.bgz', ''


# Merged sumstats
def get_hail_sumstats_path(gs_hail_gwas_path, model, fold):
    return os.path.join(gs_hail_gwas_path, f'sumstats/{model}/{fold}')


def get_saige_sumstats_root_folder(gs_gwas_path, suffix, encoding, gene_analysis):
    return os.path.join(gs_gwas_path, f'{"gene_based_sumstats" if gene_analysis else "sumstats"}/{suffix}/{encoding}')


def get_saige_sumstats_mt_folder(gs_gwas_path, suffix, encoding, gene_analysis):
    return os.path.join(get_saige_sumstats_root_folder(gs_gwas_path, suffix, encoding, gene_analysis), 'mt')


def get_saige_sumstats_mt_path(gs_gwas_path, suffix, encoding, gene_analysis, pop):
    return os.path.join(get_saige_sumstats_mt_folder(gs_gwas_path, suffix, encoding, gene_analysis), f'results_{pop}.mt')


def get_saige_meta_mt_path(gs_gwas_path, suffix, encoding, gene_analysis):
    return os.path.join(get_saige_sumstats_mt_folder(gs_gwas_path, suffix, encoding, gene_analysis), 'meta_analysis.mt')


def get_saige_cross_biobank_meta_mt_path(gs_gwas_path, suffix, encoding, gene_analysis, pop=None):
    if pop is not None:
        return os.path.join(get_saige_sumstats_mt_folder(gs_gwas_path, suffix, encoding, gene_analysis), f'cross_biobank_meta_analysis_{pop}_only.mt')
    else:
        return os.path.join(get_saige_sumstats_mt_folder(gs_gwas_path, suffix, encoding, gene_analysis), 'cross_biobank_meta_analysis.mt')


def get_saige_sumstats_tsv_folder(gs_gwas_path, suffix, encoding, gene_analysis):
    return os.path.join(get_saige_sumstats_root_folder(gs_gwas_path, suffix, encoding, gene_analysis), 'flat_files')


def get_saige_cross_biobank_meta_sumstats_tsv_folder(gs_gwas_path, suffix, encoding, gene_analysis, pop=None):
    if pop is None:
        return os.path.join(get_saige_sumstats_root_folder(gs_gwas_path, suffix, encoding, gene_analysis), 'cross_biobank_meta_flat_files')
    else:
        return os.path.join(get_saige_sumstats_root_folder(gs_gwas_path, suffix, encoding, gene_analysis), f'cross_biobank_meta_{pop}_only_flat_files')


######### IMPORT UTILS ##########
def parse_bucket(gs_bucket):
    if re.search('^gs://', gs_bucket):
        return gs_bucket
    else:
        return f'gs://{gs_bucket}'


def get_covariates_with_custom(cov_folder, custom=None, use_drc_pop=False, use_custom_pcs=None):
    new_covariates = hl.read_table(get_base_covariates_path(cov_folder, use_drc_pop=use_drc_pop, use_custom_pcs=use_custom_pcs))
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
                                      cov_folder, custom=None, use_drc_pop=True, use_custom_pcs='custom',
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
    cov_ht, cust_covar_list = get_covariates_with_custom(cov_folder, custom, use_drc_pop=use_drc_pop, use_custom_pcs=use_custom_pcs)
    cov_ht = cov_ht.persist()
    
    mt_this = mt.select_rows(**cov_ht[mt.row_key])
    mt_this = mt_this.annotate_rows(**{k: hl.int(v)  for k, v in mt_this.row.items() if v.dtype == hl.tbool})
    mt_this = mt_this.select_entries(**format_entries(mt_this.value, mt_this.sex))
    mt_this = mt_this.select_cols(**{f'n_cases_{sex}': hl.agg.count_where(
            hl.cond(mt_this.trait_type == 'categorical', mt_this[sex] == 1.0, hl.is_defined(mt_this[sex]))
         ) for sex in SEXES})
    full_mt = mt_this.unfilter_entries()
    
    return full_mt, cust_covar_list


def get_custom_ukb_pheno_mt(pheno_folder, cov_folder, custom_covariates, suffix, pop: str = 'all', 
                            use_drc_pop=True, use_custom_pcs='custom'):
    mt = hl.read_matrix_table(get_custom_ukb_pheno_mt_path(pheno_folder, suffix))
    covars, _ = get_covariates_with_custom(cov_folder=cov_folder, custom=custom_covariates, 
                                           use_drc_pop=use_drc_pop, use_custom_pcs=use_custom_pcs)
    mt = mt.annotate_rows(**covars[mt.row_key])
    mt = mt.annotate_rows(**{k: hl.int(v)  for k, v in mt.row.items() if v.dtype == hl.tbool})
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


def process_phenotype_table(phenotype_flat_file, trait_type, modifier, suffix,
                            pheno_path, covar_path, num_pcs,
                            sample_col='s', include_base_covars=True, 
                            addl_cov=None, use_drc_pop=True, use_custom_pcs='custom',
                            overwrite=False, append=False):

    curdate = date.today().strftime("%y%m%d")

    _, custom_suff = parse_drc_custom_inputs(use_drc_pop=use_drc_pop, use_custom_pcs=use_custom_pcs)

    if use_drc_pop:
        suffix = suffix + '_drcpop'
    suffix = suffix + custom_suff

    kwargs = {'data_path': phenotype_flat_file,
              'trait_type': trait_type,
              'modifier': modifier,
              'cov_folder': covar_path,
              'custom': addl_cov,
              'sample_col': sample_col,
              'use_drc_pop': use_drc_pop,
              'use_custom_pcs': use_custom_pcs}
    mt, cust_covar_list = load_custom_pheno_with_covariates(**kwargs)

    basic_covars = BASE_NONPC_COVARS if include_base_covars else []
    covariates = ','.join(basic_covars + cust_covar_list + [f'PC{x}' for x in range(1, num_pcs + 1)])

    mt_path = get_custom_ukb_pheno_mt_path(pheno_path, suffix)

    if not hl.hadoop_exists(f'{mt_path}/_SUCCESS') or (overwrite or append):
        mt_this = mt.group_rows_by('pop').aggregate(
            n_cases=hl.agg.count_where(mt.both_sexes == 1.0),
            n_controls=hl.agg.count_where(mt.both_sexes == 0.0),
            n_defined=hl.agg.count_where(hl.is_defined(mt.both_sexes))
        ).entries()
        mt_this.drop(*[x for x in PHENO_COLUMN_FIELDS if x != 'description' and x in mt_this.row]).show(100, width=180)
        
        # save table
        if append and hl.hadoop_exists(f'{mt_path}/_SUCCESS'):
            original_mt = hl.read_matrix_table(mt_path)
            original_mt = original_mt.checkpoint(get_custom_ukb_pheno_mt_path(pheno_path, f'{suffix}_before_{curdate}'), overwrite=overwrite)
            original_mt.cols().export(get_custom_phenotype_summary_backup_path(pheno_path, suffix, curdate))
            original_mt.union_cols(mt, row_join_type='outer').write(mt_path, overwrite=overwrite)
        else:
            mt.write(mt_path, overwrite=overwrite)
        
        summarize_data(pheno_path, suffix, overwrite=overwrite)
    
    return covariates


def get_cases_controls_from_logs(log_list):
    cases = controls = -1
    for log in log_list:
        try:
            with open(log, 'r') as f:
                for line in f:
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


def get_heritability_from_log(null_log, trait_type):
    import math
    heritability = -1
    with open(null_log, 'r') as f:
        for line in f:
            if line.startswith('Final'):
                fields = line.strip().split()
                if len(fields) == 4:
                    try:
                        tau = float(fields[2])
                        if trait_type == 'quantitative':
                            tau1 = float(fields[1])
                            heritability = tau / (tau1 + tau)
                        else:
                            heritability = tau / (tau + math.pi ** 2 / 3)
                        break
                    except:
                        print(f'Could not load heritability from {line}.')
    return heritability


def get_inverse_normalize_status(null_log):
    status = 'Unknown'
    with hl.hadoop_open(null_log) as f:
        for line in f:
            if line.startswith('$invNormalize'):
                try:
                    status = f.readline().strip().split()[1]
                except:
                    print(f'Could not load inv_norm status from {line} in {null_log}.')
    return status.capitalize()


def get_saige_version_from_log(null_log):
    version = 'NA'
    with open(null_log, 'r') as f:
        for line in f:
            if line.startswith('other attached packages:'):
                try:
                    line2 = f.readline()
                    packages = line2.strip().split()
                    version = [x for x in packages if 'SAIGE' in x][0]
                except:
                    print(f'Could not load version number from {line2} in {null_log}.')
    return version


def load_variant_data(output_ht_path, temp_path, paths, extension, trait_type, pheno_dict,
                      null_log, test_logs):
    
    n_cases, n_controls = get_cases_controls_from_logs(test_logs)
    heritability = get_heritability_from_log(null_log, trait_type)
    inv_normalized = get_inverse_normalize_status(null_log)
    saige_version = get_saige_version_from_log(null_log)

    print(f'Loading variant data...')
    ht = hl.import_table(paths, delimiter='\t', impute=True, min_partitions=300).checkpoint(temp_path, overwrite=True)

    print(f'Case/control counts: {str(n_cases)} cases, {str(n_controls)} controls.')
    print(f'Heritability: {str(heritability)}.')
    print(f'Inverse normalized: {str(inv_normalized)}.')
    print(f'Saige version: {str(saige_version)}.')

    if extension == 'single.txt':
        # this is common variant analysis
        marker_id_col = 'MarkerID'  
        alleles = (ht[marker_id_col].split('_')[1]).split('/')
        ht = ht.key_by(locus=hl.locus(contig=ht[marker_id_col].split(':')[0], 
                                      pos=hl.int32((ht[marker_id_col].split(':')[1]).split('_')[0]),
                                      reference_genome='GRCh38'), 
                       alleles=alleles,
                       **pheno_dict).distinct().naive_coalesce(150)
    else:
        # this is rare analysis
        marker_id_col = 'MarkerID'
        cond = ht.CHR == 'UR'
        ht = ht.key_by(gene_symbol=hl.or_missing(cond, ht[marker_id_col].split(':')[0]), 
                       group=hl.or_missing(cond, ht[marker_id_col].split(':')[1]),
                       max_MAF=hl.or_missing(cond, ht[marker_id_col].split(':')[2]),
                       locus=hl.or_missing(~cond, hl.locus(contig=ht[marker_id_col].split(':')[0], 
                                                           pos=hl.int32((ht[marker_id_col].split(':')[1]).split('_')[0]),
                                                           reference_genome='GRCh38')), 
                       alleles=hl.or_missing(~cond, (ht[marker_id_col].split('_')[1]).split('/')),
                       **pheno_dict).distinct().naive_coalesce(50)
    
    ht = ht.drop('CHR', 'POS', 'MarkerID', 'Allele1', 'Allele2')

    if n_cases == -1: n_cases = hl.null(hl.tint)
    if n_controls == -1: n_controls = hl.null(hl.tint)
    if heritability == -1.0: heritability = hl.null(hl.tfloat)
    if saige_version == 'NA': saige_version = hl.null(hl.tstr)
    if inv_normalized == 'NA': inv_normalized = hl.null(hl.tstr)
    
    ht = ht.transmute(Pvalue=ht['p.value']).annotate_globals(
        n_cases=n_cases, n_controls=n_controls, heritability=heritability, saige_version=saige_version,
        inv_normalized=inv_normalized)
    ht = ht.drop('Tstat')
    return ht.checkpoint(output_ht_path, overwrite=True)


def load_gene_data(output_ht_path, paths, trait_type, pheno_dict,
                      null_log, test_logs):
    
    n_cases, n_controls = get_cases_controls_from_logs(test_logs)
    heritability = get_heritability_from_log(null_log, trait_type)
    inv_normalized = get_inverse_normalize_status(null_log)
    saige_version = get_saige_version_from_log(null_log)

    print(f'Loading gene data ...')
    types = {x: hl.tfloat64 for x in ('Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT', 'BETA_Burden', 'SE_Burden')}
    types.update({x: hl.tint64 for x in ('Number_rare', 'Number_ultra_rare')})
    ht = hl.import_table(paths, delimiter='\t', impute=True, types=types)
    if n_cases == -1: n_cases = hl.null(hl.tint)
    if n_controls == -1: n_controls = hl.null(hl.tint)
    if heritability == -1.0: heritability = hl.null(hl.tfloat)
    if saige_version == 'NA': saige_version = hl.null(hl.tstr)
    if inv_normalized == 'NA': inv_normalized = hl.null(hl.tstr)

    ht = ht.transmute(gene_symbol = ht.Region, group = ht.Group)
    ht = ht.key_by(gene_symbol=ht.gene_symbol, group=ht.group, max_MAF=ht.max_MAF, **pheno_dict).naive_coalesce(20).annotate_globals(
        n_cases=n_cases, n_controls=n_controls, heritability=heritability, saige_version=saige_version, inv_normalized=inv_normalized)
    ht = ht.checkpoint(output_ht_path, overwrite=True)
    return ht


def make_iteration_suffix(iteration):
    if iteration == 0:
        return ''
    else:
        return f'_iter{str(iteration)}'
    

def remove_bucket(path):
    this_bucket = os.getenv('WORKSPACE_BUCKET')
    return re.sub('^'+this_bucket, '', path)


def get_pheno_dict(gs_phenotype_path, suffix, pop, min_cases=50, sex_stratified=''):
    raise NotImplementedError('Not implemented correctly.')
    # sex_stratified can be 'all', 'only', ''.
    ht = hl.read_table(get_custom_phenotype_summary_path(gs_phenotype_path, suffix))
    ht = ht.filter(ht.pop == pop)

    criteria = True
    criteria &= (ht.n_cases_by_pop >= min_cases)
    ht = ht.filter(criteria).key_by()

    if len(sex_stratified) > 0:
        ht_sex_specific = ht.annotate(pheno_sex='males').union(ht.annotate(pheno_sex='females'))
        if sex_stratified == 'all':
            ht = ht.union(ht_sex_specific)
        else:
            ht = ht_sex_specific

    out = set([tuple(x[field] for field in PHENO_KEY_FIELDS) for x in ht.select(*PHENO_KEY_FIELDS).collect()])
    pheno_key_dict = [dict(zip(PHENO_KEY_FIELDS, x)) for x in out]

    return pheno_key_dict


def get_hail_pheno_dict(gs_phenotype_path, suffix):
    pheno_ht = hl.read_matrix_table(get_custom_ukb_pheno_mt_path(gs_phenotype_path, suffix)).cols()
    key = pheno_ht.key
    ht = key._indices.source
    value = ht.row_value
    return hl.dict(ht.aggregate(hl.agg.collect((key, value)), _localize=False))


def get_all_merged_ht_paths(gs_output_path, suffix, pop, encoding, gene_analysis=False, sex_stratified=''):
    items = [x['path'] for x in hl.hadoop_ls(get_result_path(gs_output_path, suffix, pop=pop, encoding=encoding)) if x['is_dir']]

    paths_list = []
    for item in items:
        this_pheno_dict = pheno_str_to_dict(os.path.basename(item))
        this_ht, this_ht_gene = get_merged_ht_path(gs_output_path, suffix, pop, this_pheno_dict, encoding=encoding, gene_analysis=gene_analysis)
        if gene_analysis:
            this_ht = this_ht_gene
        if hl.hadoop_exists(os.path.join(this_ht, '_SUCCESS')):
            paths_list.append(this_ht)

    return paths_list


def get_modified_key(mt, enforce_pheno_cols=False):
    key = mt.col_key.annotate(phenocode=format_pheno_dir(mt.phenocode),
                              modifier=hl.case(missing_false=True)
                              #.when(mt.trait_type == "biomarkers", "")
                              .default(mt.modifier))
    if enforce_pheno_cols:
        key = key.select(*PHENO_KEY_FIELDS)
    
    return key


def get_ukb_b37_b38_liftover():
    ht = hl.import_table('gs://mito-wgs-public-2023/ukb_b37_b38_lifted_variants.tsv.bgz', min_partitions=100)
    ht = ht.select(locus_b37 = hl.locus(contig=ht.chr_b37, pos=hl.int(ht.pos_b37), reference_genome='GRCh37'),
                   alleles_b37 = [ht.ref_b37, ht.alt_b37],
                   locus_b38 = hl.locus(contig=ht.chr_b38, pos=hl.int(ht.pos_b38), reference_genome='GRCh38'),
                   alleles_b38 = [ht.ref_b38, ht.alt_b38],
                   rsid = ht.rsid).key_by('locus_b37', 'alleles_b37')
    return ht


def mwzj_hts_by_tree(all_hts, temp_dir, globals_for_col_key, debug=False, inner_mode = 'overwrite',
                     repartition_final: int = None, gene_analysis=False):


    def get_n_even_intervals(n):
        ref = hl.default_reference()
        genome_size = sum(ref.lengths.values())
        partition_size = int(genome_size / n) + 1
        return list(map(
            lambda x: hl.Interval(hl.eval(hl.locus_from_global_position(x * partition_size)),
                                hl.eval(hl.locus_from_global_position(min(x * partition_size + partition_size, genome_size - 1)))),
            range(n)))


    chunk_size = int(len(all_hts) ** 0.5) + 1
    outer_hts = []

    checkpoint_kwargs = {inner_mode: True}
    if repartition_final is not None and not gene_analysis:
        intervals = get_n_even_intervals(repartition_final)
        checkpoint_kwargs['_intervals'] = intervals

    if debug: print(f'Running chunk size {chunk_size}...')
    for i in range(chunk_size):
        if i * chunk_size >= len(all_hts): break
        hts = all_hts[i * chunk_size:(i + 1) * chunk_size]
        if debug: print(f'Going from {i * chunk_size} to {(i + 1) * chunk_size} ({len(hts)} HTs)...')
        try:
            if isinstance(hts[0], str):
                hts = list(map(lambda x: hl.read_table(x), hts))
            ht = hl.Table.multi_way_zip_join(hts, 'row_field_name', 'global_field_name')
        except:
            if debug:
                print(f'problem in range {i * chunk_size}-{i * chunk_size + chunk_size}')
                _ = [ht.describe() for ht in hts]
            raise
        outer_hts.append(ht.checkpoint(f'{temp_dir}/temp_output_{i}.ht', **checkpoint_kwargs))
    ht = hl.Table.multi_way_zip_join(outer_hts, 'row_field_name_outer', 'global_field_name_outer')
    
    if debug:
        ht = ht.checkpoint(f'{temp_dir}/temp_output_final_combined.ht', **checkpoint_kwargs)
    
    ht = ht.transmute(inner_row=hl.flatmap(lambda i:
                                           hl.cond(hl.is_missing(ht.row_field_name_outer[i].row_field_name),
                                                   hl.range(0, hl.len(ht.global_field_name_outer[i].global_field_name))
                                                   .map(lambda _: hl.null(ht.row_field_name_outer[i].row_field_name.dtype.element_type)),
                                                   ht.row_field_name_outer[i].row_field_name),
                                           hl.range(hl.len(ht.global_field_name_outer))))
    ht = ht.transmute_globals(inner_global=hl.flatmap(lambda x: x.global_field_name, ht.global_field_name_outer))
    
    mt = ht._unlocalize_entries('inner_row', 'inner_global', globals_for_col_key)
    return mt


def get_ucsc_gene_table(gs_bucket, gs_annot, overwrite):
    ht_gene_size = hl.import_table(os.path.join(gs_bucket, 'reference_files/gencode.v37.annotation.gff3.gz'), force=True, no_header=True, impute=True, comment='#').repartition(100)
    ht_gene_size = ht_gene_size.filter(ht_gene_size.f2 == 'gene')
    ht_gene_size = ht_gene_size.transmute(chrom = ht_gene_size.f0,
                                        source = ht_gene_size.f1,
                                        type = ht_gene_size.f2,
                                        start = ht_gene_size.f3,
                                        end = ht_gene_size.f4,
                                        strand = ht_gene_size.f6,
                                        size = ht_gene_size.f4 - ht_gene_size.f3,
                                        info = ht_gene_size.f8.split(';'))
    ht_gene_size = ht_gene_size.annotate(gene_symbol = ht_gene_size.info.filter(lambda x: x.startswith('gene_name=')))
    ht_gene_size = ht_gene_size.filter(hl.len(ht_gene_size.gene_symbol) > 0)
    ht_gene_size = ht_gene_size.annotate(gene_symbol = ht_gene_size.gene_symbol[0].replace('^gene_name=', '')).drop('info')
    ht_gene_size = ht_gene_size.key_by('gene_symbol')
    ht_mult = ht_gene_size.group_by(ht_gene_size.gene_symbol).aggregate(has_multiple = hl.agg.count() > 1)
    ht_gene_size = ht_gene_size.annotate(has_multiple = ht_mult[ht_gene_size.key].has_multiple)
    ht_gene_size = ht_gene_size.checkpoint(os.path.join(gs_annot, 'ucsc_gencode.v37.gene_information.ht'), _read_if_exists=not overwrite, overwrite=overwrite)

    return ht_gene_size