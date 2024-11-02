import os
import json
import hail as hl
import pandas as pd
from AoU.paths import *
from cromwell.classes import CromwellManager

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
def run_meta_analysis(mt, remove_low_confidence=True):
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
        mt = mt.annotate_entries(summary_stats = mt.summary_stats.map(lambda x: hl.if_else(x.low_confidence, hl.missing(x.dtype), x)))

    mt = mt.annotate_entries(
        summary_stats=hl.map(
            lambda x: x[1].annotate(N=hl.or_missing(hl.is_defined(x[1]), get_n(mt.pheno_data, x[0]))),
            hl.enumerate(mt.summary_stats),
        )
    )

    # Run fixed-effect meta-analysis (all + leave-one-out)
    mt = mt.annotate_entries(
        unnorm_beta=mt.summary_stats.BETA / (mt.summary_stats.SE ** 2), inv_se2=1 / (mt.summary_stats.SE ** 2)
    )
    mt = mt.annotate_entries(
        sum_unnorm_beta=all_and_leave_one_out(mt.unnorm_beta, mt.pheno_data.pop),
        sum_inv_se2=all_and_leave_one_out(mt.inv_se2, mt.pheno_data.pop),
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
    mt = mt.annotate_entries(META_N_pops=all_and_leave_one_out(mt.variant_exists, mt.pheno_data.pop))
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
    mt = mt.annotate_entries(
        #ac_cases=hl.map(lambda x: x["AF.Cases"] * x.N, mt.summary_stats),
        #ac_controls=hl.map(lambda x: x["AF.Controls"] * x.N, mt.summary_stats),
        META_AC_Allele2=all_and_leave_one_out(mt.summary_stats.AF * mt.summary_stats.N, mt.pheno_data.pop),
        META_N=all_and_leave_one_out(mt.summary_stats.N, mt.pheno_data.pop),
    )
    mt = mt.annotate_entries(
        META_AF=mt.META_AC_Allele2 / mt.META_N,
        #META_AF_Allele2=mt.META_AC_Allele2 / mt.META_N,
        #META_AF_Cases=all_and_leave_one_out(mt.ac_cases, mt.pheno_data.pop) / mt.META_N,
        #META_AF_Controls=all_and_leave_one_out(mt.ac_controls, mt.pheno_data.pop) / mt.META_N,
    )
    mt = mt.drop(
        #"unnorm_beta", "inv_se2", "variant_exists", "ac_cases", "ac_controls", "summary_stats", "META_AC_Allele2"
        "unnorm_beta", "inv_se2", "variant_exists", "summary_stats", "META_AC_Allele2"
    )

    # Format everything into array<struct>
    def is_finite_or_missing(x):
        return hl.or_missing(hl.is_finite(x), x)

    meta_fields = ["BETA", "SE", "Pvalue", "Q", "Pvalue_het", "N", "N_pops", "AF"]
    #meta_fields = ["BETA", "SE", "Pvalue", "Q", "Pvalue_het", "N", "N_pops", "AF_Allele2", "AF_Cases", "AF_Controls"]
    mt = mt.transmute_entries(
        meta_analysis=hl.map(
            lambda i: hl.struct(**{field: is_finite_or_missing(mt[f"META_{field}"][i]) for field in meta_fields}),
            hl.range(hl.len(mt.META_BETA)),
        )
    )

    col_fields = ["n_cases", "n_controls"]
    mt = mt.annotate_cols(
        **{field: all_and_leave_one_out(mt.pheno_data[field], mt.pheno_data.pop) for field in col_fields}
    )
    col_fields += ["pop"]
    mt = mt.annotate_cols(
        pop=all_and_leave_one_out(
            mt.pheno_data.pop,
            mt.pheno_data.pop,
            all_f=lambda x: x,
            loo_f=lambda i, x: hl.filter(lambda y: y != x[i], x),
        )
    )
    mt = mt.transmute_cols(
        meta_analysis_data=hl.map(
            lambda i: hl.struct(**{field: mt[field][i] for field in col_fields}), hl.range(hl.len(mt.pop))
        )
    )

    return mt


def get_lambdas_path(suffix, pop, extn):
    return os.path.join(HAIL_GWAS_PATH, f'lambdas/{suffix}/lambda_export_{pop}.{extn}')


def aou_generate_final_lambdas(mt, suffix, overwrite):
    mt = mt.annotate_cols(
        pheno_data=hl.zip(mt.pheno_data, hl.agg.array_agg(
            lambda ss: hl.agg.filter(~ss.low_confidence,
                hl.struct(lambda_gc=hl.methods.statgen._lambda_gc_agg(ss.Pvalue),
                          n_variants=hl.agg.count_where(hl.is_defined(ss.Pvalue)),
                          n_sig_variants=hl.agg.count_where(ss.Pvalue < 5e-8))),
            mt.summary_stats)).map(lambda x: x[0].annotate(**x[1]))
    )
    ht = mt.cols()
    ht = ht.checkpoint(get_lambdas_path(suffix, 'full', 'ht'), overwrite=overwrite, _read_if_exists=not overwrite)
    ht.explode('pheno_data').flatten().export(get_lambdas_path(suffix, 'full', 'txt.bgz'))
    return mt


def get_hail_sumstats_path(model, fold):
    return os.path.join(HAIL_GWAS_PATH, f'sumstats/{model}/{fold}')


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
    
    return manager.workflow_status