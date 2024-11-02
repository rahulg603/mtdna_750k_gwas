import hail as hl
import pandas as pd
from AoU.paths import *
from AoU.phenotypes import *
from AoU.sumstats import *


def apply_irnt(ht, cols):
    # similar to Neale lab round 2 approach:
    # https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/0.2/irnt.biomarkers.py
    
    df = ht.select(*cols).to_pandas()
    df.index = df['s']
    df = df.drop('s', axis=1)

    dfp = df.rank()
    dfp = (dfp - 0.5) / (~dfp.isnull()).sum()
    dfp.columns = [x + '_prob' for x in dfp.columns]
    dfp.loc[:, 's'] = dfp.index

    ht = hl.Table.from_pandas(dfp, key='s')
    ht = ht.annotate(**{x.replace('_prob', ''): hl.qnorm(ht[x])
                        for x in ht.row_value if x.endswith('_prob')})
    ht = ht.annotate(**{x: hl.or_missing(~hl.is_nan(ht[x]), ht[x]) for x in ht.row_value})
    ht = ht.drop(*[x for x in ht.row_value if x.endswith('_prob')])
    
    return ht


def run_regressions(mt, phenos, covars, pass_through, gwas_name, thresh=0, model='additive', overwrite=False, n_partition=5000):
    """
    Much of the heavy lifiting for converting the ht to mt is pulled from UKB round 2.
    Performs covariate correction for all elements in the covariates struct.
    """
    mt_dir = os.path.join(HAIL_GWAS_PATH, 'mt')
    cutoff = '' if thresh == 0 else f'_filt_prob_{thresh}'
    filename = f'{gwas_name}_gwas_{model}{cutoff}.mt'
    filename_raw = f'{gwas_name}_gwas_{model}{cutoff}_raw.mt'
    if (not overwrite) & hl.hadoop_exists(os.path.join(mt_dir,filename)):
        mt = hl.read_matrix_table(os.path.join(mt_dir,filename))

    else:
        if model == 'additive':
            entry = mt.GT.n_alt_alleles()
        elif model == 'recessive':
            entry = hl.if_else(mt.GT.n_alt_alleles() == 2, 1, 0)

        ht_dir = os.path.join(TEMP_PATH, f'{gwas_name}_gwas_{model}{cutoff}_res.ht')
        if hl.hadoop_exists(ht_dir) and not overwrite:
            ht_out = hl.read_table(ht_dir)
        else:
            ht_out = hl.linear_regression_rows(
                y=[[mt['phenotypes'][y]] for y in phenos],
                x=entry,
                covariates=[1, *[mt['covariates'][item] for item in covars]],
                pass_through=pass_through)
            ht_out = ht_out.repartition(n_partition)
            ht_out = ht_out.checkpoint(ht_dir, overwrite=True)

        ht_res = ht_out.annotate_globals(columns=hl.map(lambda i: hl.struct(phenotype=hl.literal(phenos)[i]), 
                                                        hl.range(0, hl.len(phenos))))
        ht_res = ht_res.annotate(entries=hl.map(
            lambda i: hl.struct(
                n=ht_res['n'][i],
                sum_x=ht_res['sum_x'][i],
                y_transpose_x=ht_res['y_transpose_x'][i][0],
                beta=ht_res['beta'][i][0],
                standard_error=ht_res['standard_error'][i][0],
                t_stat=ht_res['t_stat'][i][0],
                p_value=ht_res['p_value'][i][0]),
            hl.range(0, hl.len(phenos))))
        ht_res = ht_res.select(*(pass_through + ['entries']))
        mt = ht_res._unlocalize_entries('entries', 'columns', ['phenotype'])
        mt = mt.repartition(n_partition)
        mt = mt.checkpoint(os.path.join(mt_dir,filename_raw), overwrite=True)

        mt = mt.select_entries(N = mt.n,
                               AC = mt.sum_x,
                               ytx = mt.y_transpose_x,
                               BETA = mt.beta,
                               SE = mt.standard_error,
                               tstat = mt.t_stat,
                               Pvalue = mt.p_value)
        mt = mt.annotate_entries(AF = mt.AC / (2 * mt.N))
        mt = mt.annotate_entries(minor_AF = hl.cond(mt.AF <= 0.5, mt.AF, 1.0-mt.AF),
                                 minor_AC = hl.cond(mt.AF <= 0.5, mt.AC, (2 * mt.N)-mt.AC))
        mt = mt.annotate_entries(low_confidence = mt.minor_AC <= 20)

        mt = mt.checkpoint(os.path.join(mt_dir,filename), overwrite=True)
    
    return mt


def export_for_manhattan(mt, phenos, entry_keep, model, fold, suffix, overwrite, include_cols_for_mung):
    if type(entry_keep) == str:
        entry_keep = [entry_keep]

    for pheno in phenos:
        file_out = os.path.join(HAIL_GWAS_PATH,f'sumstats/{model}/{fold}/{pheno}{suffix}')
        if include_cols_for_mung:
            extra_cols = ['rsid']
        else:
            extra_cols = ['rsid']
        extra_cols = [x for x in extra_cols if x in mt.row]
        if overwrite or not hl.hadoop_exists(file_out):
            ht_f = mt.filter_cols(mt.phenotype == pheno).entries().select(*(entry_keep + extra_cols)).repartition(100)
            ht_f = ht_f.key_by()
            if include_cols_for_mung:
                ht_f = ht_f.rename({'minor_AF':'MAF', 'rsid':'SNP', 'n':'N', 
                                    'p_value': 'P', 'beta':'BETA', 
                                    'standard_error':'BETA_SE'})
                ht_f = ht_f.annotate(A1 = ht_f.alleles[1], A2 = ht_f.alleles[0])
                ht_f = ht_f.drop(ht_f.locus, ht_f.alleles, ht_f.phenotype).key_by('SNP')
            else:
                ht_f = ht_f.annotate(variant = hl.str(ht_f.locus)+ ':' + hl.delimit(ht_f.alleles, ':'))
                ht_f = ht_f.drop(ht_f.locus, ht_f.alleles, ht_f.phenotype).key_by('variant')
            ht_f.export(file_out)


def export_meta_for_manhattan(this_suffix, export_fold):
    mt = hl.read_matrix_table(get_meta_path(this_suffix))
    mt = mt.select_entries(meta_analysis=mt.meta_analysis[0]).select_cols()
    mt = mt.select_entries(**{x: mt.meta_analysis[x] for x in mt.meta_analysis.keys()})
    mt = mt.annotate_entries(minor_AF = hl.min([mt.AF, 1-mt.AF]), Pvalue = hl.exp(mt.Pvalue))
    pheno = list(mt.aggregate_cols(hl.agg.collect_as_set(mt.phenotype)))
    
    export_for_manhattan(mt=mt, phenos=pheno, entry_keep=['N','Pvalue','BETA','SE','Q','Pvalue_het','N_pops','AF', 'minor_AF'], 
                         model='additive', fold=export_fold, suffix=f'_{this_suffix("meta")}_geno_af_0.01.tsv.bgz', 
                         overwrite=True, include_cols_for_mung=False)


def filter_mt_per_pop_maf(mt, pop, cutoff, overwrite_gt, perform_per_pop_hwe=False):
    """ Expects that MT has a GT and .covariates.pop field.
    Also filters based on p_hwe (two sided) per-population.
    """
    this_pop_path = f'{TEMP_PATH}mt/genotype_mt_filtered_{pop}_maf_pass.mt'
    if hl.hadoop_exists(f'{this_pop_path}/_SUCCESS') and not overwrite_gt:
        mt_af_filt = hl.read_matrix_table(this_pop_path)
    else:
        mt_af_filt = mt.filter_cols(mt.covariates.pop == pop)
        
        if perform_per_pop_hwe:
            print('Running per-pop HWE filtering!')
            mt_af_filt = mt_af_filt.annotate_rows(hwe = hl.agg.hardy_weinberg_test(mt_af_filt.GT))
            mt_af_filt = mt_af_filt.filter_rows(mt_af_filt.hwe.p_value > 1e-10, keep = True).drop('hwe')
        
        mt_af_filt = mt_af_filt.annotate_rows(gt_stats = hl.agg.call_stats(mt_af_filt.GT, mt_af_filt.alleles))
        mt_af_filt = mt_af_filt.annotate_rows(minor_AF = hl.min(mt_af_filt.gt_stats.AF))
        mt_af_filt = mt_af_filt.filter_rows(mt_af_filt.minor_AF > cutoff).drop('gt_stats')
        mt_af_filt = mt_af_filt.checkpoint(this_pop_path, overwrite=True)
    return mt_af_filt