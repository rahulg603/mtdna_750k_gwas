import os
import secrets
import json
import hail as hl
import pandas as pd
from tqdm import tqdm
from AoU.paths import *
from utils.SaigeImporters import *
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
        print('Removing low_confidence variants prior to meta-analysis...')
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


def get_saige_sumstats_mt_folder(suffix, encoding, gene_analysis):
    return os.path.join(GWAS_PATH, f'{"gene_based_sumstats" if gene_analysis else "sumstats"}/{suffix}/{encoding}/mt')


def get_saige_sumstats_mt_path(suffix, encoding, gene_analysis, pop):
    return os.path.join(get_saige_sumstats_mt_folder(suffix, encoding, gene_analysis), f'results_{pop}.mt')


def get_saige_meta_mt_path(suffix, encoding, gene_analysis):
    return os.path.join(get_saige_sumstats_mt_folder(suffix, encoding, gene_analysis), 'meta_analysis.mt')


def custom_patch_mt_keys(mt):
    mt = mt.key_cols_by(**{x: hl.case(missing_false=True)
                        .default(mt[x])
                           for x in PHENO_KEY_FIELDS})
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


def custom_unify_saige_ht_schema(ht, path, tmp):
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
    
    ht = ht.checkpoint(os.path.join(tmp, 'unified_schema_individual_tables', secrets.token_urlsafe(12), os.path.basename(os.path.dirname(path))))
    
    return ht


def join_pheno_hts_to_mt(all_hts, row_keys, col_keys, temp_dir = None, inner_mode: str = 'overwrite',
                         repartition_final: int = None):
    
    def pull_out_col_keys(all_hts, row_keys, col_keys):
        rekeyed_hts = []
        for ht in all_hts:
            ht2 = ht.head(1)
            glob = ht2.aggregate(hl.agg.take(hl.struct(**{x: ht2[x] for x in col_keys}), 1)[0], _localize=False)
            rekeyed_hts.append(ht.key_by(*row_keys).drop(*col_keys).annotate_globals(**glob))
        return rekeyed_hts

    rekeyed_hts = pull_out_col_keys(all_hts, row_keys, col_keys)
    mt = mwzj_hts_by_tree(rekeyed_hts, temp_dir, col_keys, debug=True,
                          inner_mode=inner_mode, repartition_final=repartition_final)
    print(f'Unioned MTs...')
    return mt


def saige_generate_sumstats_mt(all_variant_outputs, pheno_dict, temp_dir, inner_mode, checkpoint, n_partitions):
    row_keys = ['locus', 'alleles', 'gene', 'annotation']
    col_keys = PHENO_KEY_FIELDS

    all_hts = [custom_unify_saige_ht_schema(hl.read_table(x), x, temp_dir) for x in tqdm(all_variant_outputs)]
    print('Schemas unified. Starting joining...')
    mt = join_pheno_hts_to_mt(all_hts, row_keys, col_keys, temp_dir=temp_dir,
                              inner_mode=inner_mode, repartition_final=n_partitions)
    print('After merge schema...')
    mt.describe()

    if checkpoint:
        mt = mt.checkpoint(f'{temp_dir}/staging.mt', **{inner_mode: True})
    
    mt = custom_patch_mt_keys(mt)
    key = mt.col_key.annotate(phenocode=format_pheno_dir(mt.phenocode))
    mt = check_and_annotate_with_dict(mt, pheno_dict, key)
    if mt.inv_normalized.dtype == hl.tstr:
        mt = mt.annotate_cols(inv_normalized=hl.bool(mt.inv_normalized))
    key = mt.col_key.annotate(phenocode=format_pheno_dir(mt.phenocode))
    mt = mt.filter_cols(mt.phenocode != "").drop('AC_Allele2', 'Is.SPA.converge')

    mt = mt.key_rows_by('locus', 'alleles')
    print('Prior to output schema...')
    mt.describe()
    
    return mt


def saige_merge_raw_sumstats(suffix, encoding, use_drc_pop, use_custom_pcs, pops=POPS, read_previous=False, overwrite=True, gene_analysis=False, n_partitions=5000):
    
    hl._set_flags(no_whole_stage_codegen="1")
    inner_mode = 'overwrite' if overwrite else '_read_if_exists'
    suffix = suffix + ('_drcpop' if use_drc_pop else '') + ('_custompcs' if use_custom_pcs == 'custom' else ('_axaoupcs' if use_custom_pcs == 'axaou' else ''))

    for pop in pops:

        merged_mt_path = get_saige_sumstats_mt_path(suffix, encoding, gene_analysis, pop)
        print(f'Final path: {merged_mt_path}')

        if read_previous and hl.hadoop_exists(f'{merged_mt_path}/_SUCCESS'):
            continue

        all_variant_outputs = get_all_merged_ht_paths(RESULTS_PATH, PHENO_PATH, suffix, pop, encoding)
        pheno_dict = get_pheno_dict(PHENO_PATH, suffix, pop, min_cases=50, sex_stratified='')

        print(f'For {suffix}, pop {pop}, {encoding}, found {str(len(pheno_dict))} phenos with {str(len(all_variant_outputs))} valid per-pheno HTs.')

        if len(all_variant_outputs) > 0:
            mt = saige_generate_sumstats_mt(all_variant_outputs, pheno_dict, temp_dir=f'{TEMP_PATH}/{suffix}/{encoding}/{pop}/variant', inner_mode=inner_mode, checkpoint=True, n_partitions=n_partitions)
            mt.write(merged_mt_path, overwrite=overwrite)

            hl.read_matrix_table(merged_mt_path).describe()

            return mt

        else:
            return None


def saige_append_merged_sumstats():
    return None


def saige_combine_per_pop_sumstats_mt():
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
    
    return manager.workflow_status