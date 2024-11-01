import os, re
from tempfile import tempdir
import pandas as pd
import hail as hl

from AoU.hail_gwas_helpers import *
from AoU.paths import *
from AoU.covariates import *
from utils.SaigeImporters import *
from AoU.munge_genotypes import get_filtered_genotype_mt


ANALYSIS_POP = ['afr','amr','eur','sas','eas','mid']
MIN_CASES = 100
IRNT = True
OVERWRITE_SUMSTATS = True
EXPORT_SUMSTATS = True
overwrite_full_gt = True


def apply_qc_continuous(mt, min_case_count: int = 50):
    mt = mt.filter_cols(mt.n_cases >= min_case_count)
    return mt


def run_full_gwas(sample_covariates, ht_pheno, num_PC, overwrite_gt, naming_insert, fold, pheno, min_cases):
    
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
                                        use_drc_ancestry_data=False)
        mt_a = mt_a.annotate_cols(phenotypes = ht_pheno[mt_a.s])
        mt_a = mt_a.annotate_cols(covariates = sample_covariates[mt_a.s])
        pheno_f = [x for x in pheno if mt_a.aggregate_cols(hl.agg.count_where(hl.is_defined(mt_a.phenotypes[x]))) > min_cases]

        # Run variant HL GWAS and export sumstats
        res = run_regressions(mt_a, pheno_f, covars, [], 
                              this_suffix(pop), overwrite=OVERWRITE_SUMSTATS)

        if EXPORT_SUMSTATS:
            export_for_manhattan(mt=res, phenos=pheno_f, entry_keep=['N','Pvalue','BETA','SE','tstat','ytx','AC','minor_AC','AF', 'minor_AF', 'low_confidence'], 
                                 model='additive', fold=fold, suffix=f'_{this_suffix(pop)}_geno_af_0.001.tsv.bgz', 
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
    full_mt = mts[0]
    for this_mt in mts[1:]:
        full_mt = full_mt.union_cols(this_mt, row_join_type='outer')
    full_mt = full_mt.checkpoint(os.path.join(TEMP_PATH, 'mt/staging_full.mt'), overwrite=True)
    full_mt = full_mt.collect_cols_by_key()
    full_mt = full_mt.checkpoint(os.path.join(TEMP_PATH, 'mt/staging_lambdas.mt'), overwrite=True)
    full_mt = aou_generate_final_lambdas(full_mt, this_suffix('full'), overwrite=True)
    full_mt.write(os.path.join(HAIL_GWAS_PATH, f'/all_pop_mt/{this_suffix("full")}.mt'), overwrite=True)


# Import covariates
gwas_covariates = get_gwas_covariates(overwrite=False, use_drc_ancestry_data=False, use_custom_data=True)

# Import phenotypes and IRNT
ht_pheno = get_case_only_mtdna_callset(num_to_keep=300, overwrite=False)
pheno = [x for x in ht_pheno.row if x not in ht_pheno.key]
if IRNT:
    ht_pheno = apply_irnt(ht_pheno, pheno)
    ht_pheno = ht_pheno.checkpoint(f'{TEMP_PATH}phenotypes_after_irnt_checkpoint.ht', overwrite=True)
    ht_pheno.export(f'{BUCKET}/analyses/221206_final_gwas_heteroplasmies_hwefix/Data/phenotypes_post_irnt.tsv')

# Run GWAS using new PCs, no iterations (raw)
fold = '241031_selected_variant_gwas_heteroplasmies'
run_full_gwas(gwas_covariates, ht_pheno, num_PC=20, overwrite_gt=True, 
              naming_insert='newPCs_iter0_hail', fold=fold, pheno=pheno, min_cases=MIN_CASES)

# AFR: 40832109
# AMR: 34642693
# EUR: 22513978
# SAS: 27266217
# EAS: 17856139

