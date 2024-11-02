"""
Setup:
pip install git+https://github.com/broadinstitute/cromshell;

pip install --upgrade pip;
pip install --upgrade aiohttp;
pip install --upgrade gnomad;

git clone https://github.com/broadinstitute/gnomad_methods.git;

cd ~;
mkdir mtSwirl_fork;
cd mtSwirl_fork;
git clone https://github.com/rahulg603/mtSwirl.git;
cd mtSwirl;
git switch refactor_merging;
cd ~;

git clone https://github.com/rahulg603/saige_aou_wdl.git
"""

import os, re
from tempfile import tempdir
import pandas as pd
import hail as hl

from AoU.hail_gwas_helpers import *
from AoU.paths import *
from AoU.covariates import *
from utils.SaigeImporters import *
from AoU.munge_genotypes import *
from AoU.sumstats import *


ANALYSIS_POP = ['afr','amr','eur','sas','eas']
#ANALYSIS_POP = ['amr']
MIN_CASES = 100
IRNT = True
OVERWRITE_SUMSTATS = False
EXPORT_SUMSTATS = False
overwrite_full_gt = True

hl.init(tmp_dir = f'{TEMP_PATH}/')


def apply_qc_continuous(mt, min_case_count: int = 50):
    mt = mt.filter_cols(mt.n_cases >= min_case_count)
    return mt


def get_this_suffix(naming_insert, irnt_suff):
    return lambda pop: f'241031_{pop}_{naming_insert}_mtdna_variant_qc_hl_case_only{irnt_suff}'


def merge_mt_list(mts, this_suffix):
    if len(mts) > 1:
        print('Starting summary statistics merging...')
        full_mt = mts[0]
        for this_mt in mts[1:]:
            full_mt = full_mt.union_cols(this_mt, row_join_type='outer')
        full_mt = full_mt.checkpoint(os.path.join(TEMP_PATH, 'mt/staging_full.mt'), overwrite=True)
        full_mt = full_mt.collect_cols_by_key()
        full_mt = full_mt.checkpoint(os.path.join(TEMP_PATH, 'mt/staging_lambdas.mt'), overwrite=True)
        full_mt = aou_generate_final_lambdas(full_mt, this_suffix('full'), overwrite=True)
        full_mt = full_mt.checkpoint(get_all_pop_mt_path(this_suffix), overwrite=True)

        full_mt.describe()

        full_mt_meta = run_meta_analysis(full_mt)
        full_mt_meta = full_mt_meta.checkpoint(get_meta_path(this_suffix))

        return full_mt
    else:
        return mts[0]


def only_merge_gwas(naming_insert):
    mts = []
    irntsuff = '_irnt' if IRNT else ''
    this_suffix = get_this_suffix(naming_insert=naming_insert, irnt_suff=irntsuff)
    for pop in ANALYSIS_POP:
        res = run_regressions(hl.MatrixTable, [], [], [], 
                              this_suffix(pop), overwrite=False)
        res = res.annotate_cols(n_cases = hl.array(hl.agg.collect_as_set(res.N)).filter(lambda x: hl.is_defined(x))[0],
                                n_controls = hl.missing(hl.tint32),
                                pop = pop,
                                inv_normalized = IRNT,
                                log_pvalue = False)
        res = apply_qc_continuous(res)
        res = res.select_cols(pheno_data=res.col_value)
        res = res.select_entries(summary_stats=res.entry)
        mts.append(res)

    # Collapse into single MT
    return merge_mt_list(mts, this_suffix)


def run_full_gwas(sample_covariates, ht_pheno, num_PC, naming_insert, fold, pheno, min_cases):
    
    # Prepare mt for regressions
    covar_full_list = list(sample_covariates.row)
    covars_base = BASE_NONPC_COVARS + \
        [f'PC{str(x)}'for x in range(1,num_PC+1)]
    covars = covars_base + [x for x in covar_full_list if re.search('^hap_', x)]
    print('Using covariates:')
    print(covars)
    irntsuff = '_irnt' if IRNT else ''
    this_suffix = lambda pop: f'241031_{pop}_{naming_insert}_mtdna_variant_qc_hl_case_only{irntsuff}'

    # Run per-population GWAS
    mts = []
    for pop in ANALYSIS_POP:
        
        # filter table by per-pop MAF (and to specific pop)
        mt_a = get_filtered_genotype_mt('variant', pop,
                                        filter_samples=True, filter_variants=True,
                                        use_array_for_variant=False,
                                        use_drc_ancestry_data=True,
                                        remove_related=True)
        mt_a = mt_a.annotate_cols(phenotypes = ht_pheno[mt_a.s])
        mt_a = mt_a.annotate_cols(covariates = sample_covariates[mt_a.s])
        pheno_f = [x for x in pheno if mt_a.aggregate_cols(hl.agg.count_where(hl.is_defined(mt_a.phenotypes[x]))) > min_cases]

        # Run variant HL GWAS and export sumstats
        res = run_regressions(mt_a, pheno_f, covars, [], 
                              this_suffix(pop), overwrite=OVERWRITE_SUMSTATS)

        if EXPORT_SUMSTATS:
            export_for_manhattan(mt=res, phenos=pheno_f, entry_keep=['N','Pvalue','BETA','SE','tstat','ytx','AC','minor_AC','AF', 'minor_AF', 'low_confidence'], 
                                 model='additive', fold=fold, suffix=f'_{this_suffix(pop)}_geno_af_0.01.tsv.bgz', 
                                 overwrite=OVERWRITE_SUMSTATS, include_cols_for_mung=False)
        
        # finish formatting MT
        res = res.annotate_cols(n_cases = hl.array(hl.agg.collect_as_set(res.N)).filter(lambda x: hl.is_defined(x))[0],
                                n_controls = hl.missing(hl.tint32),
                                pop = pop,
                                inv_normalized = IRNT,
                                log_pvalue = False)
        res = apply_qc_continuous(res)
        res = res.select_cols(pheno_data=res.col_value)
        res = res.select_entries(summary_stats=res.entry)
        mts.append(res)
    
    # Collapse into single MT
    return merge_mt_list(mts, this_suffix)


# Import covariates
gwas_covariates = get_gwas_covariates(overwrite=False, use_drc_ancestry_data=False, use_custom_data=True)
hap_covariates = get_hap_covariates('v6andv7', 'wide')
covariates = gwas_covariates.annotate(**hap_covariates[gwas_covariates.key])

# Import phenotypes
ht_pheno = get_case_only_mtdna_callset(num_to_keep=300, overwrite=False, version='v6andv7')
ht_snv = get_snv_count_phenotype(version='v6andv7')
ht_positive_control = hl.import_table(get_path_raw_positive_control(), impute=True, types={'s':hl.tstr}, key='s')

ht_pheno = ht_pheno.annotate(**ht_snv[ht_pheno.key])
ht_pheno = ht_pheno.annotate(**ht_positive_control[ht_pheno.key])
pheno_irnt = [x for x in ht_pheno.row if x in ['chrM_302_A_AC', 'chrM_302_A_ACC', 'height']]
pheno_non_irnt = ['snv_count_qcpass']
ht_pheno_for_analysis = ht_pheno.select(*pheno_irnt)
ht_pheno_non_irnt = ht_pheno.select(*pheno_non_irnt)

# IRNT phenotypes
if IRNT:
    ht_pheno_for_analysis = apply_irnt(ht_pheno_for_analysis, pheno_irnt)
    ht_pheno_for_analysis = ht_pheno_for_analysis.annotate(**ht_pheno_non_irnt[ht_pheno_for_analysis.key])
    ht_pheno_for_analysis = ht_pheno_for_analysis.checkpoint(f'{TEMP_PATH}/phenotypes_after_irnt_checkpoint.ht', overwrite=True)
    ht_pheno_for_analysis.export(os.path.join(BUCKET, f'analyses/241031_gwas_302_snvcount/Data/phenotypes_post_irnt.tsv'))

# Run GWAS using new PCs, no iterations (raw)
fold = '241031_selected_variant_gwas_heteroplasmies'
run_full_gwas(covariates, ht_pheno_for_analysis, num_PC=20,
              naming_insert='newPCs_iter0_hail', fold=fold, pheno=pheno_irnt + pheno_non_irnt, min_cases=MIN_CASES)

# export sumstats for meta analysis
irntsuff = '_irnt' if IRNT else ''
this_suffix = get_this_suffix(naming_insert='newPCs_iter0_hail', irnt_suff=irntsuff)
export_meta_for_manhattan(this_suffix, fold)


make_manhattan_plots(wdl_path='saige_aou_wdl/WDL/ManhattanPlotter.wdl', 
                     sumstat_paths = [], 
                     phenotypes = [], 
                     pops = [],
                     suffix=get_this_suffix('meta'), p_col='', af_col='', conf_col=None,
                     wid=1300, hei=640, cex=1.3, point_size=18,
                     hq_file=None,
                     exponentiate_p=False,
                     keep_x=False,
                     af_filter=None,
                     var_as_rsid=True,
                     mem=20,
                     no_wait=False)