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
IRNT = False
OVERWRITE_SUMSTATS = False
EXPORT_SUMSTATS = False

hl.init(tmp_dir = f'{TEMP_PATH}/')


def apply_qc_continuous(mt, min_case_count: int = 50):
    mt = mt.filter_cols(mt.n_cases >= min_case_count)
    return mt


def get_this_suffix(naming_insert, irnt_suff):
    return lambda pop: f'241103_{pop}_{naming_insert}_mtdna_heteroplasmic_snv_count_trunc{irnt_suff}'


def merge_mt_list(mts, this_suffix):
    if len(mts) > 1:
        print('Starting summary statistics merging...')
        print(f'Will output table to {get_all_pop_mt_path(this_suffix)}...')
        full_mt = mts[0]
        for this_mt in mts[1:]:
            full_mt = full_mt.union_cols(this_mt, row_join_type='outer')
        full_mt = full_mt.checkpoint(os.path.join(TEMP_PATH, 'mt/staging_full.mt'), overwrite=True)
        full_mt = full_mt.collect_cols_by_key()
        full_mt = full_mt.checkpoint(os.path.join(TEMP_PATH, 'mt/staging_lambdas.mt'), overwrite=True)
        full_mt = aou_generate_final_lambdas(full_mt, this_suffix('full'), overwrite=True)
        full_mt = full_mt.checkpoint(get_all_pop_mt_path(this_suffix), overwrite=True)

        full_mt.describe()
        print('Running meta analysis...')
        full_mt_meta = run_meta_analysis(full_mt)
        full_mt_meta = full_mt_meta.checkpoint(get_meta_path(this_suffix), overwrite=True)

        return full_mt
    else:
        return mts[0]


def only_merge_gwas(naming_insert, _quick_fix_for_call_rate=True):
    mts = []
    irntsuff = '_irnt' if IRNT else ''

    # we messed up the initial GWAS and overwrote the per-pop sumstats from a previous analysis. So read these in here...
    this_suffix = lambda pop: f'241031_{pop}_{naming_insert}_mtdna_variant_qc_hl_case_only{irntsuff}'

    n_samples = get_n_samples_per_pop_vec(analysis_type='variant', sample_qc=True, use_array_for_variant=False, use_drc_ancestry_data=True)
    
    for pop in ANALYSIS_POP:
        # filter table by call rate
        ht = get_call_stats_ht(pop=pop, sample_qc=True, analysis_type='variant',
                               use_drc_ancestry_data=True, 
                               use_array_for_variant=False,
                               overwrite=False)
        ht = ht.filter(
            (ht.call_stats.AN >= (n_samples[pop] * 2 * CALLRATE_CUTOFF))
            & (ht.call_stats.AC[1] > 0)
            & (ht.call_stats.AC[0] > 0)
        )

        res = run_regressions(hl.MatrixTable, [], [], [], 
                              this_suffix(pop), overwrite=False)
        
        if _quick_fix_for_call_rate:
            print('Applying hotfix call rate filter *after* performing GWAS...')
            res = res.semi_join_rows(ht)

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
    # ...but now for the merged data we will use a new prefix / suffix.
    this_suffix = get_this_suffix(naming_insert=naming_insert, irnt_suff=irntsuff)
    return merge_mt_list(mts, this_suffix)


def run_full_gwas(sample_covariates, ht_pheno, num_PC, naming_insert, fold, pheno, min_cases,
                  _quick_fix_for_call_rate=False):
    
    # Prepare mt for regressions
    covar_full_list = list(sample_covariates.row)
    covars_base = BASE_NONPC_COVARS + \
        [f'PC{str(x)}'for x in range(1,num_PC+1)]
    covars = covars_base + [x for x in covar_full_list if re.search('^hap_', x)]
    print('Using covariates:')
    print(covars)
    irntsuff = '_irnt' if IRNT else ''
    this_suffix = get_this_suffix(naming_insert=naming_insert, irnt_suff=irntsuff)

    # Run per-population GWAS
    n_samples = get_n_samples_per_pop_vec(analysis_type='variant', sample_qc=True, use_array_for_variant=False, use_drc_ancestry_data=True)
    
    mts = []
    for pop in ANALYSIS_POP:
        print(f'Performing GWAS for pop {pop}...')
        
        # filter table by per-pop MAF (and to specific pop)
        mt_a = get_filtered_genotype_mt('variant', pop,
                                        filter_samples=True, filter_variants=True,
                                        use_array_for_variant=False,
                                        use_drc_ancestry_data=True,
                                        remove_related=True)
        
        # filter table by call rate
        ht = get_call_stats_ht(pop=pop, sample_qc=True, analysis_type='variant',
                               use_drc_ancestry_data=True, 
                               use_array_for_variant=False,
                               overwrite=False)
        ht = ht.filter(
            (ht.call_stats.AN >= (n_samples[pop] * 2 * CALLRATE_CUTOFF))
            & (ht.call_stats.AC[1] > 0)
            & (ht.call_stats.AC[0] > 0)
        )

        if not _quick_fix_for_call_rate:
            mt_a = mt_a.semi_join_rows(ht)

        mt_a = mt_a.annotate_cols(phenotypes = ht_pheno[mt_a.s])
        mt_a = mt_a.annotate_cols(covariates = sample_covariates[mt_a.s])
        pheno_f = [x for x in pheno if mt_a.aggregate_cols(hl.agg.count_where(hl.is_defined(mt_a.phenotypes[x]))) > min_cases]

        # Run variant HL GWAS and export sumstats
        res = run_regressions(mt_a, pheno_f, covars, [], 
                              this_suffix(pop), overwrite=OVERWRITE_SUMSTATS)
        
        if _quick_fix_for_call_rate:
            print('Applying hotfix call rate filter *after* performing GWAS...')
            res = res.semi_join_rows(ht)

        if EXPORT_SUMSTATS:
            print('Exporting summary statistics...')
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

ht_snv = get_snv_count_phenotype(version='v6andv7')
ht_snv = ht_snv.select(snv_count_qcpass_filt4 = hl.if_else(ht_snv.snv_count_qcpass >= 5, 
                                                           hl.missing(hl.tint64), 
                                                           ht_snv.snv_count_qcpass))

pheno_irnt = []
pheno_non_irnt = ['snv_count_qcpass_filt4']

# variables
fold = '241103_gwas_snv_count_truncated'
irntsuff = '_irnt' if IRNT else ''
this_suffix = get_this_suffix(naming_insert='newPCs_iter0_hail', irnt_suff=irntsuff)

# Run GWAS using new PCs, no iterations (raw)
run_full_gwas(covariates, ht_snv, num_PC=20, _quick_fix_for_call_rate=True,
              naming_insert='newPCs_iter0_hail', fold=fold, pheno=pheno_irnt + pheno_non_irnt, min_cases=MIN_CASES)

# code to manually run merging if this wasn't run above
# only_merge_gwas('newPCs_iter0_hail', _quick_fix_for_call_rate=True)

# code to manually run a meta-analysis if this wasn't run above
# full_mt = hl.read_matrix_table(get_all_pop_mt_path(this_suffix))
# full_mt_meta = run_meta_analysis(full_mt)
# full_mt_meta = full_mt_meta.checkpoint(get_meta_path(this_suffix), overwrite=True)

# export sumstats for meta analysis
export_meta_for_manhattan(this_suffix, fold)

# make manhattan plots
# meta:
phenotypes = pheno_non_irnt + pheno_irnt
meta_suffix = this_suffix('meta')
for_paths = '_geno_af_0.01.tsv.bgz'
filenames = [x + '_' + meta_suffix + for_paths for x in phenotypes]
file_path = get_hail_sumstats_path('additive', fold)

output = make_manhattan_plots(run_name='aou_manhattan_meta_snvcount_callrate_x',
                              wdl_path='/home/jupyter/saige_aou_wdl/WDL/ManhattanPlotter.wdl', 
                              sumstat_paths = [os.path.join(file_path, x + "_" + meta_suffix + for_paths) for x in phenotypes], 
                              phenotypes = phenotypes, 
                              pops = ['meta' for _ in phenotypes],
                              suffix=meta_suffix, p_col='Pvalue', af_col='AF', conf_col=None,
                              wid=1300, hei=640, cex=1.0, point_size=18,
                              hq_file=None,
                              exponentiate_p=False,
                              keep_x=True,
                              af_filter=None,
                              var_as_rsid=True,
                              mem=60,
                              no_wait=False)