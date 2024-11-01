#!/usr/bin/env python3
from paths import *
from utils.SaigeImporters import *
from covariates import *
import hail as hl
import os

PCA_DIR = os.path.join(GENO_PATH, 'pca/')

# initialize
hl.init(default_reference="GRCh38", tmp_dir=TEMP_PATH)


def run_pca(mt: hl.MatrixTable, out_prefix: str, k=20, overwrite: bool = False):
    """
    Run PCA on a dataset
    :param mt: dataset to run PCA on
    :param out_prefix: directory and filename prefix for where to put PCA output
    :return:
    """
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(mt.GT, k=k, compute_loadings=True)
    pca_mt = mt.annotate_rows(pca_af=hl.agg.mean(mt.GT.n_alt_alleles()) / 2)
    pca_loadings = pca_loadings.annotate(pca_af=pca_mt.rows()[pca_loadings.key].pca_af)

    pca_scores.write(out_prefix + 'scores.ht', overwrite)
    pca_scores = hl.read_table(out_prefix + 'scores.ht')
    pca_scores = add_pcs(pca_scores, k)
    pca_scores.export(out_prefix + 'scores.txt.bgz')  # individual-level PCs

    pca_loadings.write(out_prefix + 'loadings.ht', overwrite)  # PCA loadings


def pc_project(
    mt: hl.MatrixTable,
    loadings_ht: hl.Table,
    loading_location: str = "loadings",
    af_location: str = "pca_af",
) -> hl.Table:
    """
    Project samples in `mt` on pre-computed PCs.

    :param mt: MT containing the samples to project
    :param loadings_ht: HT containing the PCA loadings and allele frequencies used for the PCA
    :param loading_location: Location of expression for loadings in `loadings_ht`
    :param af_location: Location of expression for allele frequency in `loadings_ht`
    :return: Table with scores calculated from loadings in column `scores`
    """
    n_variants = loadings_ht.count()

    mt = mt.annotate_rows(
        pca_loadings=loadings_ht[mt.row_key][loading_location],
        pca_af=loadings_ht[mt.row_key][af_location],
    )

    mt = mt.filter_rows(
        hl.is_defined(mt.pca_loadings)
        & hl.is_defined(mt.pca_af)
        & (mt.pca_af > 0)
        & (mt.pca_af < 1)
    )

    gt_norm = (mt.GT.n_alt_alleles() - 2 * mt.pca_af) / hl.sqrt(
        n_variants * 2 * mt.pca_af * (1 - mt.pca_af)
    )

    mt = mt.annotate_cols(scores=hl.agg.array_sum(mt.pca_loadings * gt_norm))

    return mt.cols().select("scores")


def project_individuals(pca_loadings, project_mt, k):
    """
    Project samples into predefined PCA space
    :param pca_loadings: existing PCA space
    :param project_mt: matrix table of data to project
    :param project_prefix: directory and filename prefix for where to put PCA projection output
    :return:
    """
    ht_projections = pc_project(project_mt, pca_loadings)
    ht_projections = add_pcs(ht_projections, k)
    return ht_projections


def import_ld_pruned_hq_variants(overwrite=False):
    """ Perform per-population LD pruning of the original high quality variant set from gnomAD.
    """
    final_pruned_ht_path = f'{PCA_DIR}aou_hq_variants.ht'
    if hl.hadoop_exists(f'{final_pruned_ht_path}/_SUCCESS') and not overwrite:
        ht = hl.read_table(final_pruned_ht_path)

    else:
        # munge gnomAD pre-pruned hq variants
        ht = hl.import_vcf(get_aou_hq_sites(), n_partitions=50).rows()
        ht = ht.checkpoint(final_pruned_ht_path, overwrite=overwrite)
    
    return ht


def import_aou_gt_for_pca(overwrite=False):
    """ Generates full genotype file for PCA. If MT is not found, will
    draw this file from the raw genotype data.
    """
    final_mt_path = f'{PCA_DIR}final_filtered_gt_for_PCA_qc.mt'
    
    if not overwrite and hl.hadoop_exists(f'{final_mt_path}/_SUCCESS'):
        mt = hl.read_matrix_table(final_mt_path)
    
    else:
        prelim_filt_mt = f'{TEMP_PATH}/preliminary_filt_gt_for_PCA.mt'
        if hl.hadoop_exists(f'{prelim_filt_mt}/_SUCCESS') and not overwrite:
            mt = hl.read_matrix_table(prelim_filt_mt)
        else:
            # in step 1, split multi allelics, remove variants with 0 MAF, remove variants which don't pass
            print('Generating HQ WGS MT...')
            mt = hl.read_matrix_table(WGS_PATH)
            ht_hq = import_ld_pruned_hq_variants(overwrite=overwrite) #150k sites
            mt = mt.semi_join_rows(ht_hq) # filter to HQ variants and add per-ancestry inclusion variables
            mt = mt.naive_coalesce(25000).checkpoint(prelim_filt_mt, overwrite=True) #still 150k sites

        # filter full MT and add ancestry information
        print('Appending demographics and relatedness data to WGS MatrixTable...')
        demodata = get_all_demographics()
        related_samples = hl.import_table(get_aou_related_samples(), key='sample_id')
        mt = mt.filter_cols(demodata[mt.col_key].pass_qc) # remove flagged samples
        mt = mt.annotate_cols(pop = demodata[mt.col_key].ancestry.pop, 
                              related = hl.is_defined(related_samples[mt.col_key]))
        mt = mt.naive_coalesce(10000).checkpoint(final_mt_path, overwrite=True)

    return mt


def add_pcs(ht, k, goal=20):
    """ If ht contains fewer PCs than k, adds columns with 0s to get up to k.
    """
    ht = ht.annotate(**{f'PC{i}': ht.scores[i - 1] for i in range(1, k+1)})
    if k < goal:
        ht = ht.annotate(**{f'PC{i}': 0 for i in range(1, goal+1) if f'PC{i}' not in ht.row})
    ht = ht.select(*[f'PC{i}' for i in range(1, goal+1)])
    return ht


def main(k=20, global_overwrite=False, iteration=0):
    # import full MT
    # this function will produce this PCA set if it isn't found
    mt = import_aou_gt_for_pca(overwrite=global_overwrite)
    iter_suff = make_iteration_suffix(iteration)

    # for each population, remove related samples, run PCA, and then project unrelated samples onto PCs
    hts = []
    for pop in POPS:
        print(f'Starting PCA analysis for {pop}...')
        related_prefix = f'{PCA_DIR}{pop}_results_related_aouqc{iter_suff}'
        unrelated_prefix = f'{PCA_DIR}{pop}_results_unrelated_aouqc{iter_suff}'

        file_found = hl.hadoop_exists(f'{related_prefix}.scores_projected.ht/_SUCCESS') and \
            hl.hadoop_exists(f'{unrelated_prefix}.scores.ht/_SUCCESS')
        
        if global_overwrite or not file_found:
            mtf = mt.filter_cols(mt.pop == pop)
            mtf_unrel = mtf.filter_cols(~mtf.related)
            mtf_rel = mtf.filter_cols(mtf.related)
            print(f'For pop {pop}, {str(mtf_unrel.count_cols())} unrelated samples found with {str(mtf_unrel.count_rows())} variants.')
            run_pca(mtf_unrel, f'{unrelated_prefix}.', k, True)

            pca_loadings = hl.read_table(f'{unrelated_prefix}.loadings.ht')
            ht = project_individuals(pca_loadings, mtf_rel, k)
            ht.write(f'{related_prefix}.scores_projected.ht', overwrite=True)
            hl.read_table(f'{related_prefix}.scores_projected.ht').export(
                f'{related_prefix}.scores_projected.txt.bgz')

        ht_score_rel = hl.read_table(f'{related_prefix}.scores_projected.ht')
        hts.append(ht_score_rel.annotate(pop=pop, related=True))
        ht_score_unrel = hl.read_table(f'{unrelated_prefix}.scores.ht')
        ht_score_unrel = add_pcs(ht_score_unrel, k)
        hts.append(ht_score_unrel.annotate(pop=pop, related=False))

    # merged table now contains PCs, pop, and relatedness information
    ht = hts[0].union(*hts[1:])
    ht.write(get_custom_pc_path(COVAR_PATH, iteration), overwrite=True)

    # filtering iteration 0
    # EUR: 34576 variants
    # AMR: 48101 variants
    # AFR: 65611 variants
    # EAS: 27044 variants
    # SAS: 46434 variants