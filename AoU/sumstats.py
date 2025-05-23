import os
import math
import secrets
import json
import hail as hl
import pandas as pd
from pprint import pprint
from collections import Counter
from tqdm import tqdm
from AoU.paths import *
from AoU.munge_genotypes import get_call_stats_ht, get_n_samples_per_pop_vec
from utils.SaigeImporters import *
from cromwell.classes import CromwellManager

hl.init(default_reference='GRCh38', log='combine_results.log', branching_factor=8)
hl._set_flags(no_whole_stage_codegen="1")
MAX_LAMBDA = 1.5
MIN_LAMBDA = 0

### ADAPTED FROM UPDATED PAN UKBB REPO
def all_and_leave_one_out(x, pop_array, all_f=hl.sum, loo_f=lambda i, x: hl.sum(x) - hl.or_else(x[i], 0)):
    """
    Applies a function to an input array for all populations, and for each of leave-one-out populations.

    :param x: Input array
    :param pop_array: Population array
    :param all_f: Function for all populations. It takes the input array and returns a new value
    :param loo_f: Function for each of leave-one-out populations. It takes an index of leave-one-out
                  population and the input array, and returns an array of new values.
    ...
    :return: Array of new values for all populations and for each of leave-one-out populations.
    :rtype: ArrayExpression
    """
    arr = hl.array([all_f(x)])
    arr = arr.extend(hl.map(lambda i: loo_f(i, x), hl.range(hl.len(pop_array))))
    return hl.or_missing(hl.any(hl.is_defined, x), arr)


### ADAPTED FROM UPDATED PAN UKBB REPO
def _run_inverse_variance_meta(mt, beta_col, se_col, loo_operator, suffix):
    mt = mt.annotate_entries(
        unnorm_beta=mt.summary_stats[beta_col] / (mt.summary_stats[se_col] ** 2), inv_se2=1 / (mt.summary_stats[se_col] ** 2)
    )
    mt = mt.annotate_entries(
        sum_unnorm_beta=all_and_leave_one_out(mt.unnorm_beta, mt.pheno_data[loo_operator]),
        sum_inv_se2=all_and_leave_one_out(mt.inv_se2, mt.pheno_data[loo_operator]),
    )
    mt = mt.transmute_entries(**{
        f'META_BETA{suffix}': mt.sum_unnorm_beta / mt.sum_inv_se2, 
        f'META_SE{suffix}': hl.map(lambda x: hl.sqrt(1 / x), mt.sum_inv_se2)
    })
    mt = mt.annotate_entries(**{
        f'META_Pvalue{suffix}': hl.map(lambda x: hl.log(2) + hl.pnorm(x, log_p=True), -hl.abs(mt[f'META_BETA{suffix}'] / mt[f'META_SE{suffix}']))
    })

    # Run heterogeneity test (Cochran's Q)
    mt = mt.annotate_entries(**{
        f'META_Q{suffix}': hl.map(lambda x: hl.sum((mt.summary_stats[beta_col] - x) ** 2 * mt.inv_se2), mt[f'META_BETA{suffix}']),
        'variant_exists': hl.map(lambda x: ~hl.is_missing(x), mt.summary_stats[beta_col]),
    })
    mt = mt.annotate_entries(**{f'META_N_pops': all_and_leave_one_out(mt.variant_exists, mt.pheno_data[loo_operator])})
    # filter Q-values when N_pops == 1
    mt = mt.annotate_entries(**{
        f'META_Q{suffix}': hl.map(lambda i: hl.or_missing(mt[f'META_N_pops'][i] > 1, mt[f'META_Q{suffix}'][i]), hl.range(hl.len(mt[f'META_Q{suffix}'])))
    })
    mt = mt.annotate_entries(**{
        f'META_Pvalue_het{suffix}': hl.map(
            lambda i: hl.pchisqtail(mt[f'META_Q{suffix}'][i], mt[f'META_N_pops'][i] - 1, log_p=True), hl.range(hl.len(mt[f'META_Q{suffix}']))
        )
    })
    return mt


def _run_stouffers_meta(mt, p_col, beta_col, loo_operator):
    print(f'Meta analyzing {p_col} results...')

    def _edit_pvalue(p):
        return hl.map(lambda x: hl.if_else(hl.exp(x) > 0.99, hl.log(0.99), x), p)

    two_tail = p_col == 'Pvalue_Burden'
    test_lookup = {'Pvalue_Burden': 'Burden', 'Pvalue_SKAT': 'SKAT', 'Pvalue': 'SKATO', 'Pvalue_SKATO': 'SKATO'}
    test = test_lookup[p_col]

    mt = mt.annotate_entries(**{p_col: _edit_pvalue(mt.summary_stats[p_col])})

    if two_tail:
        mt = mt.annotate_entries(**{p_col: hl.map(lambda x: x + hl.log(hl.literal(1/2)), mt.summary_stats[p_col])},
                                 **{beta_col: mt.summary_stats[beta_col]})
    else:
        mt = mt.annotate_entries(**{p_col: mt.summary_stats[p_col],
                                    beta_col: hl.map(lambda x: hl.int(hl.is_defined(x)), mt.summary_stats[p_col])})

    mt = mt.annotate_entries(**{f'weighted_Z_numerator_{test}': 
                                    hl.map(lambda x, y, z: hl.sqrt(x) * (-1 * hl.qnorm(y, log_p=True)) *  hl.sign(z),
                                           mt.summary_stats.N, mt[p_col], mt[beta_col])})
    mt = mt.annotate_entries(**{f'sum_weighted_Z_numerator_{test}': all_and_leave_one_out(mt[f'weighted_Z_numerator_{test}'], mt.pheno_data[loo_operator]),
                                f'META_N_{test}': all_and_leave_one_out(mt.summary_stats.N, mt.pheno_data[loo_operator])})

    mt = mt.annotate_entries(**{f'META_Stats_{test}': mt[f'sum_weighted_Z_numerator_{test}'] / (hl.map(lambda x: hl.sqrt(x), mt[f'META_N_{test}']))})

    if two_tail:
        mt = mt.annotate_entries(**{f'META_Pvalue_{test}': hl.map(lambda x: hl.log(hl.literal(2)) + hl.pnorm(hl.abs(x), lower_tail=False, log_p=True), mt[f'META_Stats_{test}'])})
    else:
        mt = mt.annotate_entries(**{f'META_Pvalue_{test}': hl.map(lambda x: hl.pnorm(x, lower_tail=False, log_p=True), mt[f'META_Stats_{test}'])})
    
    mt = mt.drop(f'weighted_Z_numerator_{test}', f'sum_weighted_Z_numerator_{test}', p_col, beta_col, f'META_N_{test}')
    
    return mt


def run_meta_analysis(mt, saige=True, remove_low_confidence=True, cross_biobank_meta=False, single_pop=False):
    """
    Run inverse-variance fixed-effect meta-analysis for a given MatrixTable. 
    Modified now to remove entries that are low confidence prior to meta-analysis.

    :param mt: Input MatrixTable, formatted similar to `load_final_sumstats_mt()`
    ...
    :return: Result MatrixTable with `meta_analysis` entries and `meta_analysis_data` columns.
    :rtype: MatrixTable
    """
    # Annotate per-entry sample size
    def get_n(pheno_data, i):
        return pheno_data[i].n_cases + hl.or_else(pheno_data[i].n_controls, 0)
    
    if remove_low_confidence:
        print('Removing low_confidence variants prior to meta-analysis...')
        mt = mt.annotate_entries(summary_stats=mt.summary_stats.map(
            lambda x: hl.or_missing(~x.low_confidence, x)
        ))
    mt = mt.filter_entries(~mt.summary_stats.all(lambda x: hl.is_missing(x.Pvalue)))

    if not cross_biobank_meta:
        loo_operator = 'pop'
        mt = mt.annotate_entries(
            summary_stats=hl.map(
                lambda x: x[1].annotate(N=hl.or_missing(hl.is_defined(x[1]), get_n(mt.pheno_data, x[0]))),
                hl.enumerate(mt.summary_stats),
            )
        )
    else:
        loo_operator = 'cohort'

    # Run fixed-effect meta-analysis (all + leave-one-out)
    mt = mt.annotate_entries(
        unnorm_beta=mt.summary_stats.BETA / (mt.summary_stats.SE ** 2), inv_se2=1 / (mt.summary_stats.SE ** 2)
    )
    mt = mt.annotate_entries(
        sum_unnorm_beta=all_and_leave_one_out(mt.unnorm_beta, mt.pheno_data[loo_operator]),
        sum_inv_se2=all_and_leave_one_out(mt.inv_se2, mt.pheno_data[loo_operator]),
    )
    mt = mt.transmute_entries(
        META_BETA=mt.sum_unnorm_beta / mt.sum_inv_se2, META_SE=hl.map(lambda x: hl.sqrt(1 / x), mt.sum_inv_se2)
    )
    mt = mt.annotate_entries(
        META_Pvalue=hl.map(lambda x: hl.log(2) + hl.pnorm(x, log_p=True), -hl.abs(mt.META_BETA / mt.META_SE))
    )

    # Run heterogeneity test (Cochran's Q)
    mt = mt.annotate_entries(
        META_Q=hl.map(lambda x: hl.sum((mt.summary_stats.BETA - x) ** 2 * mt.inv_se2), mt.META_BETA),
        variant_exists=hl.map(lambda x: ~hl.is_missing(x), mt.summary_stats.BETA),
    )
    mt = mt.annotate_entries(META_N_pops=all_and_leave_one_out(mt.variant_exists, mt.pheno_data[loo_operator]))
    # filter Q-values when N_pops == 1
    mt = mt.annotate_entries(
        META_Q=hl.map(lambda i: hl.or_missing(mt.META_N_pops[i] > 1, mt.META_Q[i]), hl.range(hl.len(mt.META_Q)))
    )
    mt = mt.annotate_entries(
        META_Pvalue_het=hl.map(
            lambda i: hl.pchisqtail(mt.META_Q[i], mt.META_N_pops[i] - 1, log_p=True), hl.range(hl.len(mt.META_Q))
        )
    )

    # Add other annotations
    if saige:
        if cross_biobank_meta and not single_pop:
            mt = mt.annotate_entries(
                ac_cases=hl.map(lambda x: x["AF_Cases"] * x.N, mt.summary_stats),
                ac_controls=hl.map(lambda x: x["AF_Controls"] * x.N, mt.summary_stats),
                META_AC_Allele2=all_and_leave_one_out(mt.summary_stats.AF_Allele2 * mt.summary_stats.N, mt.pheno_data[loo_operator]),
                META_N=all_and_leave_one_out(mt.summary_stats.N, mt.pheno_data[loo_operator]),
            )
        else:
            mt = mt.annotate_entries(
                ac_cases=hl.map(lambda x: x["AF.Cases"] * x.N, mt.summary_stats),
                ac_controls=hl.map(lambda x: x["AF.Controls"] * x.N, mt.summary_stats),
                META_AC_Allele2=all_and_leave_one_out(mt.summary_stats.AF_Allele2 * mt.summary_stats.N, mt.pheno_data[loo_operator]),
                META_N=all_and_leave_one_out(mt.summary_stats.N, mt.pheno_data[loo_operator]),
            )
        mt = mt.annotate_entries(
            META_AF_Allele2=mt.META_AC_Allele2 / mt.META_N,
            META_AF_Cases=all_and_leave_one_out(mt.ac_cases, mt.pheno_data[loo_operator]) / mt.META_N,
            META_AF_Controls=all_and_leave_one_out(mt.ac_controls, mt.pheno_data[loo_operator]) / mt.META_N,
        )
        mt = mt.drop(
            "unnorm_beta", "inv_se2", "variant_exists", "ac_cases", "ac_controls", "summary_stats", "META_AC_Allele2"
        )

        meta_fields = ["BETA", "SE", "Pvalue", "Q", "Pvalue_het", "N", "N_pops", "AF_Allele2", "AF_Cases", "AF_Controls"]
    
    else:
        mt = mt.annotate_entries(
            META_AC_Allele2=all_and_leave_one_out(mt.summary_stats.AF * mt.summary_stats.N, mt.pheno_data.pop),
            META_N=all_and_leave_one_out(mt.summary_stats.N, mt.pheno_data.pop),
        )
        mt = mt.annotate_entries(
            META_AF=mt.META_AC_Allele2 / mt.META_N,
        )
        mt = mt.drop(
            "unnorm_beta", "inv_se2", "variant_exists", "summary_stats", "META_AC_Allele2"
        )

        meta_fields = ["BETA", "SE", "Pvalue", "Q", "Pvalue_het", "N", "N_pops", "AF"]
    

    # Format everything into array<struct>
    def is_finite_or_missing(x):
        return hl.or_missing(hl.is_finite(x), x)

    mt = mt.transmute_entries(
        meta_analysis=hl.map(
            lambda i: hl.struct(**{field: is_finite_or_missing(mt[f"META_{field}"][i]) for field in meta_fields}),
            hl.range(hl.len(mt.META_BETA)),
        )
    )

    col_fields = ["n_cases", "n_controls"]
    mt = mt.annotate_cols(
        **{field: all_and_leave_one_out(mt.pheno_data[field], mt.pheno_data[loo_operator]) for field in col_fields}
    )
    col_fields += ["pop"]
    mt = mt.annotate_cols(
        pop=all_and_leave_one_out(
            mt.pheno_data.pop,
            mt.pheno_data[loo_operator],
            all_f=lambda x: x,
            loo_f=lambda i, x: hl.filter(lambda y: y != x[i], x),
        )
    )
    if cross_biobank_meta:
        col_fields += ["cohort"]
        mt = mt.annotate_cols(
            cohort=all_and_leave_one_out(
                mt.pheno_data.cohort,
                mt.pheno_data[loo_operator],
                all_f=lambda x: x,
                loo_f=lambda i, x: hl.filter(lambda y: y != x[i], x),
            )
        )
    mt = mt.transmute_cols(
        meta_analysis_data=hl.map(
            lambda i: hl.struct(**{field: mt[field][i] for field in col_fields}), hl.range(hl.len(mt[loo_operator]))
        )
    )

    return mt


def run_gene_meta_analysis(mt, remove_low_confidence=True, cross_biobank_meta=False):
    """
    Run meta-analysis for a given MatrixTable. 
    Modified now to remove entries that are low confidence prior to meta-analysis.
    This also now runs it on a per-gene basis, for RVAS results.

    RVAS results comprise burden test results, SKAT results, and a SKAT-O set of results.
    Burden test results output effect sizes and for these an inverse-variance scheme can be used.
    SKAT/SKAT-O do not output effect sizes and thus for these we have to use a p-value based method.

    :param mt: Input MatrixTable, formatted similar to `load_final_sumstats_mt()`
    ...
    :return: Result MatrixTable with `meta_analysis` entries and `meta_analysis_data` columns.
    :rtype: MatrixTable
    """
    # Annotate per-entry sample size
    def get_n(pheno_data, i):
        return pheno_data[i].n_cases + hl.or_else(pheno_data[i].n_controls, 0)
    
    def get_n_eff(pheno_data, i):
        n_controls = hl.or_else(pheno_data[i].n_controls, 0)
        n_cases = pheno_data[i].n_cases
        return (4 * n_cases * n_controls) / (n_cases + n_controls)

    if remove_low_confidence:
        print('Removing low_confidence genes prior to meta-analysis...')
        mt = mt.annotate_entries(summary_stats=mt.summary_stats.map(
            lambda x: hl.or_missing(~x.low_confidence, x)
        ))
    mt = mt.filter_entries(~mt.summary_stats.all(lambda x: hl.is_missing(x.Pvalue_SKATO)))

    if not cross_biobank_meta:
        loo_operator = 'pop'
        mt = mt.annotate_entries(
            summary_stats=hl.map(
                lambda x: x[1].annotate(N=hl.or_missing(hl.is_defined(x[1]), get_n(mt.pheno_data, x[0]))),
                hl.enumerate(mt.summary_stats),
            )
        )

    else:
        loo_operator = 'cohort'

    mt = mt.annotate_entries(
        summary_stats=hl.map(
            lambda x: x[1].annotate(N_eff=hl.or_missing(hl.is_defined(x[1]), get_n_eff(mt.pheno_data, x[0]))),
            hl.enumerate(mt.summary_stats),
        )
    )

    if not cross_biobank_meta:
        # Run fixed-effect meta-analysis (all + leave-one-out)
        mt = _run_inverse_variance_meta(mt, 'BETA_Burden', 'SE_Burden', loo_operator=loo_operator, suffix='_IV_Burden')

        # Run p-value meta-analysis
        mt = _run_stouffers_meta(mt, 'Pvalue_Burden', 'BETA_Burden', loo_operator=loo_operator)
        mt = _run_stouffers_meta(mt, 'Pvalue', 'BETA_SKATO', loo_operator=loo_operator)
        mt = _run_stouffers_meta(mt, 'Pvalue_SKAT', 'BETA_SKAT', loo_operator=loo_operator)

        mt = mt.annotate_entries(**{
            'variant_exists': hl.map(lambda x: ~hl.is_missing(x), mt.summary_stats.Pvalue),
        })

    else:
        # Run fixed-effect meta-analysis (all + leave-one-out) -- now with cross-biobank!
        mt = _run_inverse_variance_meta(mt, 'BETA_IV_Burden', 'SE_IV_Burden', loo_operator=loo_operator, suffix='_IV_Burden')

        # Run p-value meta-analysis
        mt = _run_stouffers_meta(mt, 'Pvalue_Burden', 'BETA_IV_Burden', loo_operator=loo_operator)
        mt = _run_stouffers_meta(mt, 'Pvalue_SKATO', 'BETA_SKATO', loo_operator=loo_operator)
        mt = _run_stouffers_meta(mt, 'Pvalue_SKAT', 'BETA_SKAT', loo_operator=loo_operator)

        mt = mt.annotate_entries(**{
            'variant_exists': hl.map(lambda x: ~hl.is_missing(x), mt.summary_stats.Pvalue_SKATO),
        })

    # Add other annotations
    mt = mt.annotate_entries(
        META_MAC=all_and_leave_one_out(mt.summary_stats.MAC, mt.pheno_data[loo_operator]),
        META_N=all_and_leave_one_out(mt.summary_stats.N, mt.pheno_data[loo_operator]),
        META_N_pops=all_and_leave_one_out(mt.variant_exists, mt.pheno_data[loo_operator])
    )
        
    mt = mt.drop(
        "unnorm_beta", "inv_se2", "variant_exists", "summary_stats"
    )

    meta_fields = ["BETA_IV_Burden", "SE_IV_Burden", "Pvalue_IV_Burden", "Q_IV_Burden", "Pvalue_het_IV_Burden", 
                   "Pvalue_Burden", "Pvalue_SKAT", "Pvalue_SKATO", 
                   "Stats_Burden", "Stats_SKAT", "Stats_SKATO", 
                   "N", "N_pops", "MAC"]
    

    # Format everything into array<struct>
    def is_finite_or_missing(x):
        return hl.or_missing(hl.is_finite(x), x)

    mt = mt.transmute_entries(
        meta_analysis=hl.map(
            lambda i: hl.struct(**{field: is_finite_or_missing(mt[f"META_{field}"][i]) for field in meta_fields}),
            hl.range(hl.len(mt.META_Pvalue_SKATO)),
        )
    )

    col_fields = ["n_cases", "n_controls"]
    mt = mt.annotate_cols(
        **{field: all_and_leave_one_out(mt.pheno_data[field], mt.pheno_data[loo_operator]) for field in col_fields}
    )
    col_fields += ["pop"]
    mt = mt.annotate_cols(
        pop=all_and_leave_one_out(
            mt.pheno_data.pop,
            mt.pheno_data[loo_operator],
            all_f=lambda x: x,
            loo_f=lambda i, x: hl.filter(lambda y: y != x[i], x),
        )
    )
    if cross_biobank_meta:
        col_fields += ["cohort"]
        mt = mt.annotate_cols(
            cohort=all_and_leave_one_out(
                mt.pheno_data.cohort,
                mt.pheno_data[loo_operator],
                all_f=lambda x: x,
                loo_f=lambda i, x: hl.filter(lambda y: y != x[i], x),
            )
        )
    mt = mt.transmute_cols(
        meta_analysis_data=hl.map(
            lambda i: hl.struct(**{field: mt[field][i] for field in col_fields}), hl.range(hl.len(mt[loo_operator]))
        )
    )

    return mt


def get_hail_lambdas_path(suffix, pop, encoding, extn):
    return os.path.join(HAIL_GWAS_PATH, f'lambdas/{suffix}/{encoding}/lambda_export_{pop}.{extn}')


def get_saige_lambdas_path(suffix, pop, encoding, extn):
    return os.path.join(GWAS_PATH, f'lambdas/{suffix}/{encoding}/lambda_export_{pop}.{extn}')


def aou_generate_final_lambdas(mt, suffix, encoding, overwrite, cross_biobank=False, saige=True, exp_p=False, gene_analysis=False, specific_pop=None):
    # enable specific_pop to generate lambdas for a specific population. Only impacts the naming scheme for the output table.
    #   NOTE: does not impact the usage of low_confidence of array operations!
    if exp_p:
        transf = lambda x: hl.exp(x)
    else:
        transf = lambda x: x
    pval = 'Pvalue_SKATO' if gene_analysis and cross_biobank else 'Pvalue'
    keyword = 'genes' if gene_analysis else 'variants'
    
    if cross_biobank and specific_pop is None:
        mt = mt.annotate_cols(
            pheno_data=hl.zip(mt.pheno_data, hl.agg.array_agg(
                lambda ss: hl.struct(**{'lambda_gc': hl.methods.statgen._lambda_gc_agg(transf(ss[pval]), approximate=True),
                                        f'n_{keyword}': hl.agg.count_where(hl.is_defined(ss[pval])),
                                        f'n_sig_{keyword}': hl.agg.count_where(transf(ss[pval]) < 5e-8)}),
                mt.summary_stats)).map(lambda x: x[0].annotate(**x[1]))
        )
    else:
        mt = mt.annotate_cols(
            pheno_data=hl.zip(mt.pheno_data, hl.agg.array_agg(
                lambda ss: hl.agg.filter(~ss.low_confidence,
                    hl.struct(**{'lambda_gc': hl.methods.statgen._lambda_gc_agg(transf(ss[pval]), approximate=True),
                                 f'n_{keyword}': hl.agg.count_where(hl.is_defined(ss[pval])),
                                 f'n_sig_{keyword}': hl.agg.count_where(transf(ss[pval]) < 5e-8)})),
                mt.summary_stats)).map(lambda x: x[0].annotate(**x[1]))
        )
    ht = mt.cols()

    if saige:
        get_lambdas_path = get_saige_lambdas_path
    else:
        get_lambdas_path = get_hail_lambdas_path
    
    pop = specific_pop if specific_pop is not None else 'full'
    term = f'{pop}_cross_biobank' if cross_biobank else pop
    term = f'gene_{term}' if gene_analysis else term

    # trial to show variability
    #ht.explode('pheno_data').pheno_data.lambda_gc.show(50)
    #ht.explode('pheno_data').pheno_data.lambda_gc.show(50)
    
    ht = ht.checkpoint(get_lambdas_path(suffix, term, encoding, 'ht'), overwrite=overwrite, _read_if_exists=not overwrite)
    ht.explode('pheno_data').flatten().export(get_lambdas_path(suffix, term, encoding, 'txt.bgz'))
    
    return mt


def custom_patch_mt_keys(mt, gene_analysis=False):
    key_set = PHENO_KEY_FIELDS + GENE_COL_KEY_FIELDS if gene_analysis else PHENO_KEY_FIELDS
    mt = mt.key_cols_by(**{x: hl.case(missing_false=True)
                        .default(mt[x])
                           for x in key_set})
    return mt


def check_and_annotate_with_dict(mt, input_dict, dict_key_from_mt, axis='cols'):
    direction = mt.col if axis == 'cols' else mt.row
    annotation = {}
    for key, value in list(input_dict.values()[0].items()):
        if key in list(direction):
            annotation[key] = hl.case().when(
                hl.is_defined(mt[key]), mt[key]).when(
                hl.is_defined(input_dict.get(dict_key_from_mt)), input_dict.get(dict_key_from_mt)[key]).default(
                hl.null(value.dtype)
            )
        else:
            annotation[key] = hl.if_else(hl.is_defined(input_dict.get(dict_key_from_mt)),
                                         input_dict.get(dict_key_from_mt)[key],
                                         hl.null(value.dtype))
    if axis == 'cols':
        mt = mt.annotate_cols(**annotation)
    else:
        mt = mt.annotate_rows(**annotation)
    return mt


def unify_saige_ht_schema(ht, path, tmp, row_keys, col_keys):
    """

    :param Table ht:
    :param str patch_case_control_count: Path to file (hack to get cases and controls back if loading later)
    :return:
    :rtype: Table
    """
    #assert ht.head(1).annotation.collect()[0] is None, f'failed at {patch_case_control_count}'
    if 'gene' not in list(ht.row):
        ht = ht.annotate(gene = hl.missing('tstr'))
    if 'annotation' not in list(ht.row):
        ht = ht.annotate(annotation = hl.missing('tstr'))
    
    if 'AF_case' not in list(ht.row):
        ht = ht.select('AC_Allele2', 'AF_Allele2', 'MissingRate', 'N', 'BETA', 'SE',
                       **{'p.value.NA': hl.null(hl.tfloat64), 'Is.SPA.converge': hl.null(hl.tint32),
                          'var': ht.var, 'AF.Cases': hl.null(hl.tfloat64),
                          'AF.Controls': hl.null(hl.tfloat64), 'Pvalue': ht.Pvalue,
                          'gene': hl.or_else(ht.gene, ''), 'annotation': hl.or_else(ht.annotation, '')})
    else:
        ht = ht.rename({'Is.SPA': 'Is.SPA.converge', 'AF_case':'AF.Cases', 'AF_ctrl':'AF.Controls'})
        ht = ht.annotate_globals(n_cases = ht.head(1).N_case.collect()[0], 
                                 n_controls = ht.head(1).N_ctrl.collect()[0])
        ht = ht.select('AC_Allele2', 'AF_Allele2', 'MissingRate', 
                       **{'N': ht.N_case + ht.N_ctrl, 'BETA':ht.BETA, 'SE':ht.SE,
                          'p.value.NA':ht['p.value.NA'], 'Is.SPA.converge':hl.null(hl.tint32),#'Is.SPA.converge':ht['Is.SPA.converge'], 
                          'var':ht.var, 'AF.Cases':ht['AF.Cases'], 'AF.Controls':ht['AF.Controls'], 
                          'Pvalue':ht.Pvalue, 'gene': hl.or_else(ht.gene, ''), 'annotation':hl.or_else(ht.annotation, '')})
    
    ht = ht.annotate(BETA = hl.float64(ht.BETA), SE = hl.float64(ht.SE))

    if 'heritability' not in list(ht.globals):
        ht = ht.annotate_globals(heritability=hl.null(hl.tfloat64))
    if 'saige_version' not in list(ht.globals):
        ht = ht.annotate_globals(saige_version=hl.null(hl.tstr))
    
    ht2 = ht.head(1)
    glob = ht2.aggregate(hl.agg.take(hl.struct(**{x: ht2[x] for x in col_keys}), 1)[0], _localize=False)
    ht = ht.key_by(*row_keys).drop(*col_keys).annotate_globals(**glob)
    
    ht = ht.checkpoint(os.path.join(tmp, 'unified_schema_individual_tables', secrets.token_urlsafe(12), os.path.basename(os.path.dirname(path))))
    
    return ht


def unify_saige_gene_ht_schema(ht, path, tmp, row_keys, col_keys):
    """

    :param Table ht:
    :param str patch_case_control_count: Path to file (hack to get cases and controls back if loading later)
    :return:
    :rtype: Table
    """
    #assert ht.head(1).annotation.collect()[0] is None, f'failed at {patch_case_control_count}'

    if 'AF_case' not in list(ht.row):
        ht = ht.select('Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT', 'BETA_Burden', 'SE_Burden', 'MAC', 'Number_rare', 'Number_ultra_rare')
    else:
        raise NotImplementedError('ERROR: case/control RVAS importing is not yet implemented.')
        ht = ht.rename({'Is.SPA': 'Is.SPA.converge', 'AF_case':'AF.Cases', 'AF_ctrl':'AF.Controls'})
        ht = ht.annotate_globals(n_cases = ht.head(1).N_case.collect()[0], 
                                 n_controls = ht.head(1).N_ctrl.collect()[0])
        ht = ht.select('AC_Allele2', 'AF_Allele2', 'MissingRate', 
                       **{'N': ht.N_case + ht.N_ctrl, 'BETA':ht.BETA, 'SE':ht.SE,
                          'p.value.NA':ht['p.value.NA'], 'Is.SPA.converge':hl.null(hl.tint32),#'Is.SPA.converge':ht['Is.SPA.converge'], 
                          'var':ht.var, 'AF.Cases':ht['AF.Cases'], 'AF.Controls':ht['AF.Controls'], 
                          'Pvalue':ht.Pvalue, 'gene': hl.or_else(ht.gene, ''), 'annotation':hl.or_else(ht.annotation, '')})
    
    if 'heritability' not in list(ht.globals):
        ht = ht.annotate_globals(heritability=hl.null(hl.tfloat64))
    if 'saige_version' not in list(ht.globals):
        ht = ht.annotate_globals(saige_version=hl.null(hl.tstr))
    
    ht2 = ht.head(1)
    glob = ht2.aggregate(hl.agg.take(hl.struct(**{x: ht2[x] for x in col_keys}), 1)[0], _localize=False)
    ht = ht.key_by()
    ht = ht.annotate(max_MAF = hl.if_else(ht.group == 'Cauchy', 'NA', hl.str(ht.max_MAF)))
    ht = ht.key_by(*row_keys).drop(*col_keys).annotate_globals(**glob)
    
    ht = ht.checkpoint(os.path.join(tmp, 'gene_unified_schema_individual_tables', secrets.token_urlsafe(12), os.path.basename(os.path.dirname(path))))
    
    return ht


def saige_apply_qc(mt, filter_sumstats, this_pop_N, case_ac_threshold: int = 3, overall_mac_threshold: int = 20, min_case_count: int = 50,
                   min_call_rate: float = CALLRATE_CUTOFF):
    
    if filter_sumstats:
        mt = mt.filter_cols(mt.n_cases >= min_case_count)
    
    ac_cases = mt.n_cases * 2 * mt['AF.Cases']
    an_controls = mt.n_controls * 2 * mt['AF.Controls']

    maf_total = 0.5 - hl.abs(0.5 - mt['AF_Allele2'])
    an_total = (mt.n_cases + hl.or_else(mt.n_controls, 0)) * 2
    mac_total = maf_total * an_total

    if min_call_rate is not None:
        return mt.annotate_entries(
            low_confidence=hl.case(missing_false=True)
                .when(ac_cases <= case_ac_threshold, True)
                .when(an_controls <= case_ac_threshold, True)
                .when(mac_total <= overall_mac_threshold, True)
                .when(mt.overall_AN < 2 * this_pop_N * min_call_rate, True)
                .default(False)
        ) 
    else:
        return mt.annotate_entries(
            low_confidence=hl.case(missing_false=True)
                .when(ac_cases <= case_ac_threshold, True)
                .when(an_controls <= case_ac_threshold, True)
                .when(mac_total <= overall_mac_threshold, True)
                .default(False)
        ) 


def saige_apply_gene_qc(mt, filter_sumstats, case_ac_threshold: int = 3, overall_mac_threshold: int = 20, min_case_count: int = 50, bypass_missing=False):
    
    if filter_sumstats:
        if bypass_missing:
            mt = mt.filter_cols(hl.is_missing(mt.n_cases) | (mt.n_cases >= min_case_count))
        else:
            mt = mt.filter_cols(mt.n_cases >= min_case_count)
    
    #ac_cases = mt.n_cases * 2 * mt['AF.Cases']
    #an_controls = mt.n_controls * 2 * mt['AF.Controls']
    #maf_total = 0.5 - hl.abs(0.5 - mt['AF_Allele2'])
    #an_total = (mt.n_cases + hl.or_else(mt.n_controls, 0)) * 2
    mac_total = mt['MAC']

    return mt.annotate_entries(
        low_confidence=hl.case(missing_false=True)
            #.when(ac_cases <= case_ac_threshold, True)
            #.when(an_controls <= case_ac_threshold, True)
            .when(mac_total <= overall_mac_threshold, True)
            #.when(mt.overall_AN < 2 * this_pop_N * min_call_rate, True)
            .default(False)
    )


def split_initial_saige_gene_ht(ht, new_row_keys, new_col_keys):
    htk = ht.key_by(*new_col_keys).collect_by_key()
    htk_indexed = htk.add_index("idx")
    n = htk_indexed.count()
    final_holder = []
    for i in range(n):
        this_ht = htk_indexed.filter(htk_indexed.idx == i).explode('values').drop('idx')
        this_ht = this_ht.transmute(**this_ht.values.flatten())
        ht2 = this_ht.head(1)
        glob = ht2.aggregate(hl.agg.take(hl.struct(**{x: ht2[x] for x in new_col_keys}), 1)[0], _localize=False)
        this_ht = this_ht.key_by(*new_row_keys).drop(*new_col_keys).annotate_globals(**glob)
        final_holder.append(this_ht)
    
    return final_holder


def saige_generate_sumstats_mt(all_variant_outputs, pheno_dict, temp_dir, inner_mode, checkpoint, n_partitions, gene_analysis):
    
    if gene_analysis:        
        row_keys = ['gene_symbol']
        initial_col_keys = PHENO_KEY_FIELDS
        initial_hts = [unify_saige_gene_ht_schema(hl.read_table(x), x, temp_dir, row_keys, initial_col_keys) for x in tqdm(all_variant_outputs)]
        new_col_keys = GENE_COL_KEY_FIELDS
        row_keys = GENE_ROW_KEY_FIELDS
        col_keys = initial_col_keys + new_col_keys
        all_hts = [y for x in initial_hts for y in split_initial_saige_gene_ht(x, row_keys, new_col_keys)]
    
    else:
        col_keys = PHENO_KEY_FIELDS
        row_keys = ['locus', 'alleles', 'gene', 'annotation']
        all_hts = [unify_saige_ht_schema(hl.read_table(x), x, temp_dir, row_keys, col_keys) for x in tqdm(all_variant_outputs)]
    
    print('Schemas unified. Starting joining...')

    mt = mwzj_hts_by_tree(all_hts, temp_dir, col_keys, debug=True, gene_analysis=gene_analysis,
                          inner_mode=inner_mode, repartition_final=n_partitions)
    
    print(f'Unioned MTs...')
    print('After merge schema...')
    mt.describe()

    if checkpoint:
        mt = mt.checkpoint(f'{temp_dir}/staging.mt', **{inner_mode: True})
    
    mt = custom_patch_mt_keys(mt, gene_analysis=gene_analysis)
    key = mt.col_key.annotate(phenocode=format_pheno_dir(mt.phenocode)).select(*PHENO_KEY_FIELDS)
    mt = check_and_annotate_with_dict(mt, pheno_dict, key)
    if mt.inv_normalized.dtype == hl.tstr:
        mt = mt.annotate_cols(inv_normalized=hl.bool(mt.inv_normalized))
    key = mt.col_key.annotate(phenocode=format_pheno_dir(mt.phenocode))

    if gene_analysis:
        mt = mt.filter_cols(mt.phenocode != "")
        mt = mt.key_rows_by(*row_keys)
    else:
        mt = mt.filter_cols(mt.phenocode != "").drop('AC_Allele2', 'Is.SPA.converge')
        mt = mt.key_rows_by('locus', 'alleles')
    
    print('Prior to output schema...')
    mt.describe()
    
    return mt


def saige_merge_raw_sumstats(suffix, encoding, use_drc_pop, use_custom_pcs, pops=POPS, read_previous=False, overwrite=True, gene_analysis=False, n_partitions=1500):
    
    inner_mode = 'overwrite' if overwrite else '_read_if_exists'
    suffix = update_suffix(suffix, use_drc_pop, use_custom_pcs)

    for pop in pops:

        merged_mt_path = get_saige_sumstats_mt_path(GWAS_PATH, suffix, encoding, gene_analysis, pop)
        print(f'Final path: {merged_mt_path}')

        if read_previous and hl.hadoop_exists(f'{merged_mt_path}/_SUCCESS'):
            continue

        all_variant_outputs = get_all_merged_ht_paths(RESULTS_PATH, suffix, pop, encoding, gene_analysis)
        pheno_dict = get_hail_pheno_dict(PHENO_PATH, suffix)

        print(f'For {suffix}, pop {pop}, {encoding}, found {str(hl.len(pheno_dict).collect()[0])} phenos with {str(len(all_variant_outputs))} valid per-pheno HTs.')

        if len(all_variant_outputs) > 0:
            mt = saige_generate_sumstats_mt(all_variant_outputs, pheno_dict, temp_dir=f'{TEMP_PATH}/{suffix}/{encoding}/{pop}/{"variant" if not gene_analysis else "gene"}', 
                                            inner_mode=inner_mode, checkpoint=True, n_partitions=n_partitions, gene_analysis=gene_analysis)
            mt.write(merged_mt_path, overwrite=overwrite)

    return None


def saige_append_merged_sumstats():
    return None


def saige_combine_per_pop_gene_mt(suffix, encoding, use_drc_pop, use_custom_pcs, pops=POPS, overwrite=False, 
                                  skip_producing_lambdas=False, filter_sumstats=True, _debug_read=False):
    temp_dir = f'{TEMP_PATH}/{suffix}/{encoding}/'
    staging_full = f'{temp_dir}/staging_full.mt'
    suffix = update_suffix(suffix, use_drc_pop, use_custom_pcs)

    def reannotate_cols(mt, suffix):
        pheno_dict = get_hail_pheno_dict(PHENO_PATH, suffix)
        key = get_modified_key(mt, enforce_pheno_cols=True)
        mt = check_and_annotate_with_dict(mt, pheno_dict, key)
        return mt

    def re_colkey_mt(mt):
        mt = mt.key_cols_by().select_cols(*PHENO_KEY_FIELDS, *GENE_COL_KEY_FIELDS, *(set(PHENO_COLUMN_FIELDS) - set(PHENO_DESCRIPTION_FIELDS)), *PHENO_GWAS_FIELDS, 'pop')
        return mt.key_cols_by(*PHENO_KEY_FIELDS, *GENE_COL_KEY_FIELDS)

    if _debug_read and hl.hadoop_exists(staging_full + '/_SUCCESS'):
        full_mt = hl.read_matrix_table(staging_full)
    else:
        mts = []

        for pop in pops:
            mt = hl.read_matrix_table(get_saige_sumstats_mt_path(GWAS_PATH, suffix, encoding, gene_analysis=True, pop=pop)).annotate_cols(pop=pop)
            mt = mt.key_rows_by(*GENE_ROW_KEY_FIELDS)
            mt = mt.annotate_cols(_logged=hl.agg.any(mt.Pvalue < 0))
            mt = mt.annotate_entries(Pvalue=hl.if_else(mt._logged, mt.Pvalue, hl.log(mt.Pvalue)),
                                     Pvalue_Burden=hl.if_else(mt._logged, mt.Pvalue_Burden, hl.log(mt.Pvalue_Burden)),
                                     Pvalue_SKAT=hl.if_else(mt._logged, mt.Pvalue_SKAT, hl.log(mt.Pvalue_SKAT))).drop('_logged')

            mt = saige_apply_gene_qc(mt, filter_sumstats)
            mt = custom_patch_mt_keys(mt, gene_analysis=True)
            mt = reannotate_cols(mt, suffix)
            mt = re_colkey_mt(mt)
            mt = mt.select_cols(pheno_data=mt.col_value)
            mt = mt.select_entries(summary_stats=mt.entry)
            mts.append(mt)

        full_mt = mts[0]
        for mt in mts[1:]:
            full_mt = full_mt.union_cols(mt, row_join_type='outer')
        
        full_mt = full_mt.checkpoint(f'{temp_dir}/staging_full.mt', overwrite=True)

    full_mt = full_mt.collect_cols_by_key()
    full_mt = full_mt.annotate_cols(
        pheno_data=full_mt.pheno_data.map(lambda x: x.drop(*(set(PHENO_COLUMN_FIELDS) - set(PHENO_DESCRIPTION_FIELDS)))),
        **{f'n_cases_full_cohort_{sex}': full_mt.pheno_data[f'n_cases_{sex}'][0]
           for sex in ('both_sexes', 'females', 'males')}
    )

    staging_lambda = f'{temp_dir}/staging_lambdas.mt'
    if not skip_producing_lambdas:
        if hl.hadoop_exists(staging_lambda + '/_SUCCESS') and _debug_read:
            full_mt = hl.read_matrix_table(staging_lambda)
        else:
            full_mt = full_mt.checkpoint(f'{temp_dir}/staging_lambdas.mt', overwrite=True)
        full_mt = aou_generate_final_lambdas(full_mt, suffix, encoding=encoding, overwrite=True, exp_p=True, gene_analysis=True)
        if filter_sumstats:
            full_mt = full_mt.annotate_entries(summary_stats = hl.zip(full_mt.summary_stats, full_mt.pheno_data.lambda_gc
                                                                     ).filter(lambda x: (x[1] < MAX_LAMBDA) & (x[1] >= MIN_LAMBDA)
                                                                     ).map(lambda x: x[0]))
            full_mt = full_mt.annotate_cols(pheno_data = full_mt.pheno_data.filter(lambda x: (x.lambda_gc < MAX_LAMBDA) & (x.lambda_gc >= MIN_LAMBDA)))
            
    full_mt = full_mt.checkpoint(get_saige_sumstats_mt_path(GWAS_PATH, suffix, encoding, gene_analysis=True, pop='full'), overwrite)

    full_mt_meta = run_gene_meta_analysis(full_mt)
    full_mt_meta = full_mt_meta.checkpoint(get_saige_meta_mt_path(GWAS_PATH, suffix, encoding, gene_analysis=True), overwrite)
    
    print('Pops per pheno:')
    pprint(dict(Counter(full_mt_meta.aggregate_cols(hl.agg.counter(hl.len(full_mt_meta.pheno_data))))))

    return None


def saige_combine_aou_ukb_gene_mt(suffix, encoding, ukb_path, filter_sumstats, 
                                  overwrite=True, use_drc_pop=True, use_custom_pcs='custom', append_ukb_modifier='_irnt'):
    # NOTE: this method assumes that the UKB analysis is EUR only. This is due to extremely limited sample size in this cohort.

    def munge_mt_for_merge(mt, cohort):
        mt = mt.select_cols(cohort = cohort, 
                            n_cases = mt.meta_analysis_data[0].n_cases,
                            n_controls = mt.meta_analysis_data[0].n_controls,
                            saige_version = hl.array(hl.set(mt.pheno_data.saige_version))[0],
                            inv_normalized = hl.array(hl.set(mt.pheno_data.inv_normalized))[0],
                            pop = mt.meta_analysis_data[0].pop,
                            n_genes = hl.agg.count_where(hl.is_defined(mt.meta_analysis[0])),
                            n_sig_genes = hl.agg.count_where(mt.meta_analysis[0].Pvalue_SKATO <  math.log(5e-8)))
        mt = mt.select_entries(**mt.meta_analysis[0], low_confidence = False)

        mt = mt.select_cols(pheno_data=mt.col_value)
        mt = mt.select_entries(summary_stats=mt.entry)
        return mt

    def re_colkey_mt(mt):
        mt = mt.key_cols_by().select_cols(*PHENO_KEY_FIELDS, *GENE_COL_KEY_FIELDS, *PHENO_GWAS_FIELDS, 'pop')
        return mt.key_cols_by(*PHENO_KEY_FIELDS, *GENE_COL_KEY_FIELDS)

    suffix = update_suffix(suffix, use_drc_pop, use_custom_pcs)
    temp_dir = f'{TEMP_PATH}/{suffix}/{encoding}/'
    meta_path = get_saige_cross_biobank_meta_mt_path(GWAS_PATH, suffix, encoding, gene_analysis=True)

    if overwrite or not hl.hadoop_exists(os.path.join(meta_path, '_SUCCESS')):

        # load mts
        mt_aou = hl.read_matrix_table(get_saige_meta_mt_path(GWAS_PATH, suffix, encoding, gene_analysis=True))
        mt_ukb = hl.read_matrix_table(ukb_path)
        
        # modify ukb to be compatible with aou
        mt_ukb = mt_ukb.key_cols_by()
        mt_ukb = mt_ukb.annotate_cols(modifier = mt_ukb.modifier + append_ukb_modifier,
                                      pop=['EUR'])
        mt_ukb = mt_ukb.key_cols_by(*PHENO_KEY_FIELDS, *GENE_COL_KEY_FIELDS)

        mt_ukb = mt_ukb.annotate_cols(_logged=hl.agg.any(mt_ukb.Pvalue < 0))
        mt_ukb = mt_ukb.annotate_entries(Pvalue=hl.if_else(mt_ukb._logged, mt_ukb.Pvalue, hl.log(mt_ukb.Pvalue)),
                                         Pvalue_Burden=hl.if_else(mt_ukb._logged, mt_ukb.Pvalue_Burden, hl.log(mt_ukb.Pvalue_Burden)),
                                         Pvalue_SKAT=hl.if_else(mt_ukb._logged, mt_ukb.Pvalue_SKAT, hl.log(mt_ukb.Pvalue_SKAT))).drop('_logged')

        mt_ukb = saige_apply_gene_qc(mt_ukb, filter_sumstats, bypass_missing=False)
        mt_ukb = custom_patch_mt_keys(mt_ukb, gene_analysis=True)
        mt_ukb = re_colkey_mt(mt_ukb)
        mt_ukb = mt_ukb.select_cols(cohort = 'ukb', 
                                    n_cases = mt_ukb.n_cases,
                                    n_controls = mt_ukb.n_controls,
                                    saige_version = mt_ukb.saige_version,
                                    inv_normalized = mt_ukb.inv_normalized,
                                    pop=mt_ukb.pop,
                                    n_genes=hl.agg.count_where(hl.is_defined(mt_ukb.Pvalue)),
                                    n_sig_genes=hl.agg.count_where(hl.exp(mt_ukb.Pvalue) < 5e-8))

        # we now create new fields to use to perform meta-analysis. The burden P can be used with both Stouffer's and the IV test.
        mt_ukb = mt_ukb.select_entries(BETA_IV_Burden = mt_ukb.BETA_Burden,
                                       SE_IV_Burden = mt_ukb.SE_Burden,
                                       Pvalue_IV_Burden = mt_ukb.Pvalue_Burden,
                                       Q_IV_Burden = hl.missing(hl.tfloat64),
                                       Pvalue_het_IV_Burden = hl.missing(hl.tfloat64),
                                       Pvalue_Burden = mt_ukb.Pvalue_Burden,
                                       Pvalue_SKAT = mt_ukb.Pvalue_SKAT,
                                       Pvalue_SKATO = mt_ukb.Pvalue,
                                       Stats_Burden = hl.missing(hl.tfloat64),
                                       Stats_SKAT = hl.missing(hl.tfloat64),
                                       Stats_SKATO = hl.missing(hl.tfloat64),
                                       N = mt_ukb.n_cases + hl.or_else(mt_ukb.n_controls, 0),
                                       N_pops = 1,
                                       MAC = mt_ukb.MAC,
                                       low_confidence = mt_ukb.low_confidence)
        mt_ukb = mt_ukb.select_cols(pheno_data=mt_ukb.col_value)
        mt_ukb = mt_ukb.select_entries(summary_stats=mt_ukb.entry)

        # combine the mts so that each row has 2 items
        mt_aou_munged = munge_mt_for_merge(mt_aou, cohort='aou')
        full_mt = mt_ukb.union_cols(mt_aou_munged, row_join_type='outer').naive_coalesce(2000).checkpoint(os.path.join(TEMP_PATH, 'staging_gene_cross_biobank_joint.mt'), overwrite=True)
        full_mt = full_mt.collect_cols_by_key()
        staging_lambda = os.path.join(temp_dir, 'staging_cross_biobank_before_lambda_meta.mt')

        if hl.hadoop_exists(staging_lambda + '/_SUCCESS') and not overwrite:
            full_mt = hl.read_matrix_table(staging_lambda)
        else:
            full_mt = full_mt.checkpoint(staging_lambda, overwrite=True)
        full_mt = aou_generate_final_lambdas(full_mt, suffix, encoding=encoding, cross_biobank=True, overwrite=True, exp_p=True, gene_analysis=True)
        if filter_sumstats:
            full_mt = full_mt.annotate_entries(summary_stats = hl.zip(full_mt.summary_stats, full_mt.pheno_data.lambda_gc
                                                                    ).filter(lambda x: (x[1] < MAX_LAMBDA) & (x[1] >= MIN_LAMBDA)
                                                                    ).map(lambda x: x[0]))
            full_mt = full_mt.annotate_cols(pheno_data = full_mt.pheno_data.filter(lambda x: (x.lambda_gc < MAX_LAMBDA) & (x.lambda_gc >= MIN_LAMBDA)))

        # feed into meta analysis generator
        full_mt_meta = run_gene_meta_analysis(full_mt, remove_low_confidence=True, cross_biobank_meta=True)
        full_mt_meta = full_mt_meta.checkpoint(meta_path, overwrite=overwrite)


    print('Meta analysis complete.')
    return None


def saige_combine_per_pop_sumstats_mt(suffix, encoding, use_drc_pop, use_custom_pcs, sample_qc=True, pops=POPS, overwrite=False,
                                      skip_producing_lambdas=False, min_call_rate=CALLRATE_CUTOFF, filter_sumstats=True, _debug_read=False):
    temp_dir = f'{TEMP_PATH}/{suffix}/{encoding}/'
    staging_full = f'{temp_dir}/staging_full.mt'
    suffix = update_suffix(suffix, use_drc_pop, use_custom_pcs)

    def reannotate_cols(mt, suffix):
        pheno_dict = get_hail_pheno_dict(PHENO_PATH, suffix)
        key = get_modified_key(mt)
        mt = check_and_annotate_with_dict(mt, pheno_dict, key)
        return mt

    def re_colkey_mt(mt):
        mt = mt.key_cols_by().select_cols(*PHENO_KEY_FIELDS, *(set(PHENO_COLUMN_FIELDS) - set(PHENO_DESCRIPTION_FIELDS)), *PHENO_GWAS_FIELDS, 'pop')
        return mt.key_cols_by(*PHENO_KEY_FIELDS)

    if _debug_read and hl.hadoop_exists(staging_full + '/_SUCCESS'):
        full_mt = hl.read_matrix_table(staging_full)
    else:
        mts = []

        per_pop_N = get_n_samples_per_pop_vec(sample_qc=sample_qc, use_drc_pop=use_drc_pop,
                                              analysis_type='variant', 
                                              use_array_for_variant=False)
        for pop in pops:
            mt = hl.read_matrix_table(get_saige_sumstats_mt_path(GWAS_PATH, suffix, encoding, gene_analysis=False, pop=pop)).annotate_cols(pop=pop)
            mt = mt.key_rows_by('locus', 'alleles')
            mt = mt.annotate_cols(_logged=hl.agg.any(mt.Pvalue < 0))
            mt = mt.annotate_entries(Pvalue=hl.if_else(mt._logged, mt.Pvalue, hl.log(mt.Pvalue))).drop('_logged')

            call_stats = get_call_stats_ht(pop=pop, sample_qc=sample_qc, use_drc_pop=use_drc_pop,
                                           analysis_type='variant', use_array_for_variant=False,
                                           overwrite=False)
            mt = mt.annotate_rows(overall_AN = call_stats[mt.locus, mt.alleles].call_stats.AN)

            mt = saige_apply_qc(mt, filter_sumstats, min_call_rate=min_call_rate, this_pop_N=per_pop_N[pop])
            mt = custom_patch_mt_keys(mt, gene_analysis=False)
            mt = reannotate_cols(mt, suffix)
            mt = re_colkey_mt(mt)
            mt = mt.select_cols(pheno_data=mt.col_value)
            mt = mt.select_entries(summary_stats=mt.entry)
            mts.append(mt)

        full_mt = mts[0]
        for mt in mts[1:]:
            full_mt = full_mt.union_cols(mt, row_join_type='outer')
        
        full_mt = full_mt.checkpoint(f'{temp_dir}/staging_full.mt', overwrite=True)

    full_mt = full_mt.collect_cols_by_key()
    full_mt = full_mt.annotate_cols(
        pheno_data=full_mt.pheno_data.map(lambda x: x.drop(*(set(PHENO_COLUMN_FIELDS) - set(PHENO_DESCRIPTION_FIELDS)))),
        **{f'n_cases_full_cohort_{sex}': full_mt.pheno_data[f'n_cases_{sex}'][0]
           for sex in ('both_sexes', 'females', 'males')}
    )

    staging_lambda = f'{temp_dir}/staging_lambdas.mt'
    if not skip_producing_lambdas:
        if hl.hadoop_exists(staging_lambda + '/_SUCCESS') and _debug_read:
            full_mt = hl.read_matrix_table(staging_lambda)
        else:
            full_mt = full_mt.checkpoint(f'{temp_dir}/staging_lambdas.mt', overwrite=True)
        full_mt = aou_generate_final_lambdas(full_mt, suffix, encoding=encoding, overwrite=True, exp_p=True, gene_analysis=False)

        if filter_sumstats:
            full_mt = full_mt.annotate_entries(summary_stats = hl.zip(full_mt.summary_stats, full_mt.pheno_data.lambda_gc
                                                                     ).filter(lambda x: (x[1] < MAX_LAMBDA) & (x[1] > MIN_LAMBDA)
                                                                     ).map(lambda x: x[0]))
            full_mt = full_mt.annotate_cols(pheno_data = full_mt.pheno_data.filter(lambda x: (x.lambda_gc < MAX_LAMBDA) & (x.lambda_gc > MIN_LAMBDA)))

    full_mt = full_mt.checkpoint(get_saige_sumstats_mt_path(GWAS_PATH, suffix, encoding, gene_analysis=False, pop='full'), overwrite)

    full_mt_meta = run_meta_analysis(full_mt, saige=True)
    full_mt_meta = full_mt_meta.checkpoint(get_saige_meta_mt_path(GWAS_PATH, suffix, encoding, gene_analysis=False), overwrite)
    
    print('Pops per pheno:')
    pprint(dict(Counter(full_mt_meta.aggregate_cols(hl.agg.counter(hl.len(full_mt_meta.pheno_data))))))

    return None


def saige_combine_aou_ukb_sumstats_mt(suffix, encoding, gene_analysis, ukb_meta_path, filter_sumstats, 
                                      overwrite=True, use_drc_pop=True, use_custom_pcs='custom', append_ukb_modifier='_irnt'):


    def munge_mt_for_merge(mt, cohort):
        mt = mt.select_cols(cohort = cohort, 
                            n_cases = mt.meta_analysis_data[0].n_cases,
                            n_controls = mt.meta_analysis_data[0].n_controls,
                            saige_version = hl.array(hl.set(mt.pheno_data.saige_version))[0],
                            inv_normalized = hl.array(hl.set(mt.pheno_data.inv_normalized))[0],
                            pop = mt.meta_analysis_data[0].pop,
                            n_variants = hl.agg.count_where(hl.is_defined(mt.meta_analysis[0])),
                            n_sig_variants = hl.agg.count_where(mt.meta_analysis[0].Pvalue <  math.log(5e-8)))
        mt = mt.select_entries(**mt.meta_analysis[0])

        mt = mt.select_cols(pheno_data=mt.col_value)
        mt = mt.select_entries(summary_stats=mt.entry)
        return mt


    suffix = update_suffix(suffix, use_drc_pop, use_custom_pcs)
    temp_dir = f'{TEMP_PATH}/{suffix}/{encoding}/'
    meta_path = get_saige_cross_biobank_meta_mt_path(GWAS_PATH, suffix, encoding, gene_analysis)

    if overwrite or not hl.hadoop_exists(os.path.join(meta_path, '_SUCCESS')):

        # load mts
        mt_aou = hl.read_matrix_table(get_saige_meta_mt_path(GWAS_PATH, suffix, encoding, gene_analysis))
        mt_ukb = hl.read_matrix_table(ukb_meta_path)
        # 'gs://rgupta-assoc/saige_gwas/sumstats/241215_case_only_heteroplasmy_500k_fixed/mt/meta_analysis.mt'
        
        # modify ukb to be compatible with aou
        ht_liftover = get_ukb_b37_b38_liftover().checkpoint(os.path.join(BUCKET, 'ukb_500k', 'ukb_lift_b37_b38.ht'), _read_if_exists=True)
        mt_ukb = mt_ukb.key_cols_by()
        mt_ukb = mt_ukb.annotate_cols(modifier = mt_ukb.modifier + append_ukb_modifier)
        mt_ukb = mt_ukb.key_cols_by(*PHENO_KEY_FIELDS)
        mt_ukb = mt_ukb.annotate_rows(**ht_liftover[mt_ukb.row_key]).key_rows_by()
        mt_ukb = mt_ukb.transmute_rows(locus_grch37 = mt_ukb.locus, alleles_grch37 = mt_ukb.alleles)
        mt_ukb = mt_ukb.transmute_rows(locus = mt_ukb.locus_b38, alleles = mt_ukb.alleles_b38).key_rows_by('locus','alleles')
        mt_ukb = mt_ukb.checkpoint(os.path.join(temp_dir, 'ukb_heteroplasmy_meta_modified.mt'), _read_if_exists=not overwrite, overwrite=overwrite)

        # confirm that the number of shared variants and traits make sense
        #mt_aou.semi_join_cols(mt_ukb.cols()).count_cols() # 65
        #mt_aou.semi_join_rows(mt_ukb.rows()).count() # 24657342 (of 28152304)

        # find the highest power analysis for each trait
        # it actually looks like the meta table pulls in the highest power pop automatically
        # we assume that low confidence variants are pre-filtered

        # combine the mts so that each row has 2 items
        mt_ukb_munged = munge_mt_for_merge(mt_ukb, cohort='ukb')
        mt_aou_munged = munge_mt_for_merge(mt_aou, cohort='aou')
        full_mt = mt_ukb_munged.union_cols(mt_aou_munged, row_join_type='outer').naive_coalesce(4000).checkpoint(os.path.join(TEMP_PATH, 'staging_cross_biobank_joint.mt'), overwrite=True)
        full_mt = full_mt.collect_cols_by_key()
        staging_lambda = os.path.join(temp_dir, 'staging_cross_biobank_before_lambda_meta.mt')

        if hl.hadoop_exists(staging_lambda + '/_SUCCESS') and not overwrite:
            full_mt = hl.read_matrix_table(staging_lambda)
        else:
            full_mt = full_mt.checkpoint(staging_lambda, overwrite=True)
        full_mt = aou_generate_final_lambdas(full_mt, suffix, encoding=encoding, cross_biobank=True, overwrite=True, exp_p=True)
        if filter_sumstats:
            full_mt = full_mt.annotate_entries(summary_stats = hl.zip(full_mt.summary_stats, full_mt.pheno_data.lambda_gc
                                                                    ).filter(lambda x: (x[1] < MAX_LAMBDA) & (x[1] > MIN_LAMBDA)
                                                                    ).map(lambda x: x[0]))
            full_mt = full_mt.annotate_cols(pheno_data = full_mt.pheno_data.filter(lambda x: (x.lambda_gc < MAX_LAMBDA) & (x.lambda_gc > MIN_LAMBDA)))

        # feed into meta analysis generator
        full_mt_meta = run_meta_analysis(full_mt, saige=True, remove_low_confidence=False, cross_biobank_meta=True)
        full_mt_meta = full_mt_meta.checkpoint(meta_path, overwrite=overwrite)


    print('Meta analysis complete.')
    return None


def saige_combine_aou_ukb_single_pop_sumstats_mt(pop, suffix, encoding, ukb_mt_path, filter_sumstats, overwrite=True, sample_qc=True, 
                                                 skip_producing_lambdas=False, min_call_rate=CALLRATE_CUTOFF,
                                                 use_drc_pop=True, use_custom_pcs='custom', append_ukb_modifier='_irnt'):

    temp_dir = f'{TEMP_PATH}/{suffix}/{encoding}/ukbaou/{pop}/'
    staging_full = f'{temp_dir}/staging_full.mt'
    suffix = update_suffix(suffix, use_drc_pop, use_custom_pcs)
    meta_path = get_saige_cross_biobank_meta_mt_path(GWAS_PATH, suffix, encoding, gene_analysis=False, pop=pop)

    def reannotate_cols(mt, suffix):
        pheno_dict = get_hail_pheno_dict(PHENO_PATH, suffix)
        key = get_modified_key(mt)
        mt = check_and_annotate_with_dict(mt, pheno_dict, key)
        return mt

    def re_colkey_mt(mt):
        mt = mt.key_cols_by().select_cols(*PHENO_KEY_FIELDS, *(set(PHENO_COLUMN_FIELDS) - set(PHENO_DESCRIPTION_FIELDS)), *PHENO_GWAS_FIELDS, 'pop', 'cohort')
        return mt.key_cols_by(*PHENO_KEY_FIELDS)

    if not overwrite and hl.hadoop_exists(staging_full + '/_SUCCESS'):
        full_mt = hl.read_matrix_table(staging_full)
    else:
        mts = []

        # start with AoU
        per_pop_N = get_n_samples_per_pop_vec(sample_qc=sample_qc, use_drc_pop=use_drc_pop,
                                              analysis_type='variant', 
                                              use_array_for_variant=False)
        mt = hl.read_matrix_table(get_saige_sumstats_mt_path(GWAS_PATH, suffix, encoding, gene_analysis=False, pop=pop)).annotate_cols(pop=pop, cohort='aou')
        mt = mt.key_rows_by('locus', 'alleles')
        mt = mt.annotate_cols(_logged=hl.agg.any(mt.Pvalue < 0))
        mt = mt.annotate_entries(Pvalue=hl.if_else(mt._logged, mt.Pvalue, hl.log(mt.Pvalue))).drop('_logged')

        call_stats = get_call_stats_ht(pop=pop, sample_qc=sample_qc, use_drc_pop=use_drc_pop,
                                        analysis_type='variant', use_array_for_variant=False,
                                        overwrite=False)
        mt = mt.annotate_rows(overall_AN = call_stats[mt.locus, mt.alleles].call_stats.AN)

        mt = saige_apply_qc(mt, filter_sumstats, min_call_rate=min_call_rate, this_pop_N=per_pop_N[pop])
        mt = custom_patch_mt_keys(mt, gene_analysis=False)
        mt = reannotate_cols(mt, suffix)
        mt = re_colkey_mt(mt)
        mt = mt.annotate_entries(MissingRate = hl.float64(mt.MissingRate))
        mt = mt.select_cols(pheno_data=mt.col_value)
        mt = mt.select_entries(summary_stats=mt.entry)
        mt = mt.select_rows('gene', 'annotation')
        mts.append(mt)

        # now UKB
        mt_ukb = hl.read_matrix_table(ukb_mt_path).annotate_cols(pop=pop, cohort='ukb')
        mt_ukb = mt_ukb.key_rows_by('locus', 'alleles')
        mt_ukb = mt_ukb.annotate_cols(_logged=hl.agg.any(mt_ukb.Pvalue < 0))
        mt_ukb = mt_ukb.annotate_entries(Pvalue=hl.if_else(mt_ukb._logged, mt_ukb.Pvalue, hl.log(mt_ukb.Pvalue))).drop('_logged')

        ht_liftover = get_ukb_b37_b38_liftover().checkpoint(os.path.join(BUCKET, 'ukb_500k', 'ukb_lift_b37_b38.ht'), _read_if_exists=True)
        mt_ukb = mt_ukb.key_cols_by()
        mt_ukb = mt_ukb.annotate_cols(modifier = mt_ukb.modifier + append_ukb_modifier)
        mt_ukb = mt_ukb.key_cols_by(*PHENO_KEY_FIELDS)
        mt_ukb = mt_ukb.annotate_rows(**ht_liftover[mt_ukb.row_key]).key_rows_by()
        mt_ukb = mt_ukb.transmute_rows(locus_grch37 = mt_ukb.locus, alleles_grch37 = mt_ukb.alleles)
        mt_ukb = mt_ukb.transmute_rows(locus = mt_ukb.locus_b38, alleles = mt_ukb.alleles_b38).key_rows_by('locus','alleles')
        mt_ukb = mt_ukb.checkpoint(os.path.join(temp_dir, f'ukb_gwas_{pop}_modified.mt'), _read_if_exists=not overwrite, overwrite=overwrite)

        mt_ukb = saige_apply_qc(mt_ukb, filter_sumstats, min_call_rate=None, this_pop_N=None)
        mt_ukb = custom_patch_mt_keys(mt_ukb, gene_analysis=False)
        mt_ukb = re_colkey_mt(mt_ukb)
        mt_ukb = mt_ukb.annotate_entries(MissingRate = hl.float64(mt_ukb.MissingRate))
        mt_ukb = mt_ukb.select_cols(pheno_data=mt_ukb.col_value)
        mt_ukb = mt_ukb.select_entries(summary_stats=mt_ukb.entry)
        mt_ukb = mt_ukb.select_rows('gene', 'annotation')
        mts.append(mt_ukb)

        full_mt = mts[0]
        for mt in mts[1:]:
            full_mt = full_mt.union_cols(mt, row_join_type='outer')
        
        full_mt = full_mt.checkpoint(f'{temp_dir}/staging_full.mt', overwrite=True)

    full_mt = full_mt.collect_cols_by_key()
    full_mt = full_mt.annotate_cols(
        pheno_data=full_mt.pheno_data.map(lambda x: x.drop(*(set(PHENO_COLUMN_FIELDS) - set(PHENO_DESCRIPTION_FIELDS)))),
        **{f'n_cases_full_cohort_{sex}': full_mt.pheno_data[f'n_cases_{sex}'][0]
           for sex in ('both_sexes', 'females', 'males')}
    )

    staging_lambda = f'{temp_dir}/staging_{pop}_cross_biobank_before_lambdas.mt'
    if not skip_producing_lambdas:
        if hl.hadoop_exists(staging_lambda + '/_SUCCESS') and not overwrite:
            full_mt = hl.read_matrix_table(staging_lambda)
        else:
            full_mt = full_mt.checkpoint(staging_lambda, overwrite=True)
        # single pop indicates a single population cross-biobank meta
        full_mt = aou_generate_final_lambdas(full_mt, suffix, encoding=encoding, cross_biobank=True, specific_pop=pop, overwrite=True, exp_p=True, gene_analysis=False)
        if filter_sumstats:
            full_mt = full_mt.annotate_entries(summary_stats = hl.zip(full_mt.summary_stats, full_mt.pheno_data.lambda_gc
                                                                     ).filter(lambda x: (x[1] < MAX_LAMBDA) & (x[1] > MIN_LAMBDA)
                                                                     ).map(lambda x: x[0]))
            full_mt = full_mt.annotate_cols(pheno_data = full_mt.pheno_data.filter(lambda x: (x.lambda_gc < MAX_LAMBDA) & (x.lambda_gc > MIN_LAMBDA)))
            
    full_mt_meta = run_meta_analysis(full_mt, saige=True, remove_low_confidence=True, cross_biobank_meta=True, single_pop=True)
    full_mt_meta = full_mt_meta.checkpoint(meta_path, overwrite)
    
    print('Cohorts per pheno:')
    pprint(dict(Counter(full_mt_meta.aggregate_cols(hl.agg.counter(hl.len(full_mt_meta.pheno_data))))))

    return None


def make_manhattan_plots(wdl_path, 
                         sumstat_paths: list, phenotypes: list, pops: list,
                         suffix, p_col, af_col, conf_col=None,
                         run_name='aou_manhattan_plotting',
                         wid=1300, hei=640, cex=1.3, point_size=18,
                         hq_file=None,
                         exponentiate_p=False,
                         keep_x=True,
                         af_filter=None,
                         var_as_rsid=True,
                         mem=30,
                         no_wait=False):
    """
    The manhattan plotter expects a single per-ancestry summary statistic.

    PARAMETERS
    ----------
    wdl_path: local path to ManhattanPlotter.wdl
    sumstat_paths: paths to summary statistic flat files
    phenotypes: names of each phenotype
    pop: population of each analysis
    suffix: suffix of analysis
    p_col: column name of p-value column
    af_col: column name of AF column
    conf_col: column name of "low_confidence" column. Can be None.
    hq_file: path to a file linking variants to "high quality" measure. Optional.
    exponentiate_p: if True, will perform exp(p-value)
    keep_x: if True, will keep the x-chromosome. Otherwise, will drop it
    af_filter: If provided, will filter summary stats based on AF. Optional; will not filter if not provided
    var_as_rsid: if True, will name "variant" based on combination of "chr" "pos" "ref" "alt" columns
    mem: memory
    no_wait: if true will not wait for jobs to finish
    """
    # make json
    this_run = {'sumstats': sumstat_paths, 'pheno': phenotypes, 'pop': pops}
    df = pd.DataFrame(this_run)
    df.to_csv(os.path.abspath('./this_run.tsv'), index=False, sep='\t') # ADD THIS

    baseline = {'ManhattanPlotter.suffix': suffix,
                'ManhattanPlotter.p_col': p_col,
                'ManhattanPlotter.af_col': af_col,
                'ManhattanPlotter.wid': wid,
                'ManhattanPlotter.hei': hei,
                'ManhattanPlotter.cex': cex,
                'ManhattanPlotter.point_size': point_size,
                'ManhattanPlotter.exponentiate_p': exponentiate_p,
                'ManhattanPlotter.keep_x': keep_x,
                'ManhattanPlotter.var_as_rsid': var_as_rsid,
                'ManhattanPlotter.mem': mem}
    if conf_col is not None:
        baseline.update({'ManhattanPlotter.conf_col': conf_col})
    if hq_file is not None:
        baseline.update({'ManhattanPlotter.hq_file': hq_file})
    if af_filter is not None:
        baseline.update({'ManhattanPlotter.af_filter': af_filter})

    with open(os.path.abspath('./saige_template.json'), 'w') as j:
        json.dump(baseline, j)

    # run manhattan plotting
    print('Starting manhattan plotting...')
    print('This stage will use Cromwell.')
    manager = CromwellManager(run_name=run_name,
                              inputs_file=df,
                              json_template_path=os.path.abspath('./saige_template.json'),
                              wdl_path=wdl_path,
                              batch=None, limit=None, n_parallel_workflows=99, 
                              add_requester_pays_parameter=False,
                              restart=True, batches_precomputed=False, 
                              submission_sleep=0, check_freq=120, quiet=False)
    manager.run_pipeline(submission_retries=0, cromwell_timeout=60, skip_waiting=no_wait)
    
    return manager


def proximity_clump_sumstats(wdl_path,
                             sumstat_paths: list, 
                             phenotypes: list,
                             suffix: str, 
                             p_col: str,
                             p_thresh=5e-5,
                             window_kb=100,
                             chr='chr', pos='pos', ref='ref', alt='alt',
                             conf_col=None,
                             n_partitions=20,
                             n_cpu=16,
                             exponentiate_p=True,
                             build='GRCh38',
                             run_name='aou_clumping',
                             no_wait=False):
    """

    PARAMETERS
    ----------
    wdl_path: local path to clump.wdl
    sumstat_paths: paths to summary statistic flat files
    phenotypes: names of each phenotype
    suffix: suffix of analysis
    p_col: column name of p-value column
    conf_col: column name of "low_confidence" column. Can be None.
    exponentiate_p: if True, will perform exp(p-value)
    no_wait: if true will not wait for jobs to finish
    """
    # make json
    this_run = {'sumstats': sumstat_paths, 'pheno': phenotypes}
    df = pd.DataFrame(this_run)
    df.to_csv(os.path.abspath('./this_run.tsv'), index=False, sep='\t') # ADD THIS

    if build=='GRCh38':
        gene_path = 'gs://mito-wgs-public-free/NCBI38_ensembl.gene.loc'
    else:
        gene_path = 'gs://mito-wgs-public-free/NCBI37_ensembl.gene.loc'

    baseline = {'clump_sumstats.suffix': suffix,
                'clump_sumstats.p_col': p_col,
                'clump_sumstats.p_thresh': p_thresh,
                'clump_sumstats.window_kb': window_kb,
                'clump_sumstats.chr': chr,
                'clump_sumstats.pos': pos,
                'clump_sumstats.ref': ref,
                'clump_sumstats.alt': alt,
                'clump_sumstats.n_partitions': n_partitions,
                'clump_sumstats.exp_p': exponentiate_p,
                'clump_sumstats.gene_file': gene_path,
                'clump_sumstats.ref_genome': build,
                'clump_sumstats.cpu': n_cpu}
    if conf_col is not None:
        baseline.update({'clump_sumstats.conf_col': conf_col})

    with open(os.path.abspath('./saige_template.json'), 'w') as j:
        json.dump(baseline, j)

    print('Starting clumping-by-proximity...')
    print('This stage will use Cromwell.')
    manager = CromwellManager(run_name=run_name,
                              inputs_file=df,
                              json_template_path=os.path.abspath('./saige_template.json'),
                              wdl_path=wdl_path,
                              batch=None, limit=None, n_parallel_workflows=99, 
                              add_requester_pays_parameter=False,
                              restart=True, batches_precomputed=False, 
                              submission_sleep=0, check_freq=120, quiet=False)
    manager.run_pipeline(submission_retries=0, cromwell_timeout=60, skip_waiting=no_wait)
    
    return manager