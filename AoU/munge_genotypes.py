#!/usr/bin/env python3
# This file includes functions to process AoU genomes.
# All pathing to new files is provided by SaigeImporters.py.
# All pathing to AoU files is provided by paths.py.

# Need to determine sample QC and variant QC.
import argparse
import hail as hl
import pandas as pd

from AoU.paths import *
from utils.SaigeImporters import *
from AoU.covariates import get_all_demographics
from copy import deepcopy
from typing import Union


def get_adj_expr(
    gt_expr: hl.expr.CallExpression,
    gq_expr: Union[hl.expr.Int32Expression, hl.expr.Int64Expression],
    ad_expr: hl.expr.ArrayNumericExpression,
    adj_gq: int = 30,
    adj_ab: float = 0.2,
) -> hl.expr.BooleanExpression:
    """
    Get adj genotype annotation.

    Defaults correspond to gnomAD values.
    """
    return (
        (gq_expr >= adj_gq)
        & (
            hl.case()
            .when(~gt_expr.is_het(), True)
            .when(gt_expr.is_het_ref(), ad_expr[gt_expr[1]] / hl.sum(ad_expr) >= adj_ab)
            .default(
                (ad_expr[gt_expr[0]] / hl.sum(ad_expr) >= adj_ab)
                & (ad_expr[gt_expr[1]] / hl.sum(ad_expr) >= adj_ab)
            )
        )
    )


def annotate_adj(
        mt: hl.MatrixTable,
        adj_gq: int = 30,
        adj_ab: float = 0.2,
) -> hl.MatrixTable:
    """
    Annotate genotypes with adj criteria (assumes diploid).

    Defaults correspond to gnomAD values.
    """
    if "GT" not in mt.entry and "LGT" in mt.entry:
        print("No GT field found, using LGT instead.")
        gt_expr = mt.LGT
    else:
        gt_expr = mt.GT

    if "AD" not in mt.entry and "LAD" in mt.entry:
        print("No AD field found, using LAD instead.")
        ad_expr = mt.LAD
    else:
        ad_expr = mt.AD

    return mt.annotate_entries(
        adj=get_adj_expr(
            gt_expr, mt.GQ, ad_expr, adj_gq, adj_ab
        )
    )


def get_n_samples_vec(analysis_type, sample_qc, use_array_for_variant, use_drc_ancestry_data=False):
    mt = get_filtered_genotype_mt(analysis_type=analysis_type,
                                  pop='all',
                                  filter_samples=sample_qc,
                                  filter_variants=False,
                                  use_array_for_variant=use_array_for_variant,
                                  use_drc_ancestry_data=use_drc_ancestry_data)
    ht = mt.cols()
    ht_ct = ht.group_by(ht.pop).aggregate(hl.agg.count())
    df_ct = ht_ct.to_pandas()


def get_filtered_genotype_mt(analysis_type, pop,
                             filter_samples=True, filter_variants=True,
                             use_array_for_variant=False,
                             use_drc_ancestry_data=False):
    if analysis_type == 'gene':
        mt_path = EXOME_PATH
    elif analysis_type == 'variant':
        mt_path = ARRAY_PATH if use_array_for_variant else WGS_PATH
    else:
        raise ValueError('ERROR: analysis_type can only be gene or variant.')
    
    mt = hl.read_matrix_table(mt_path)
    if (analysis_type == 'gene') or (analysis_type == 'variant' and not use_array_for_variant):
        mt = mt.filter_entries(hl.is_missing(mt.FT) | (mt.FT == 'PASS'))
        mt = mt.drop('variant_qc')
    else:
        mt = mt.repartition(9000)
    
    meta_ht = get_all_demographics(use_drc_ancestry_data=use_drc_ancestry_data)
    mt = mt.annotate_cols(**meta_ht[mt.col_key])
    mt = mt.annotate_cols(pop = mt.ancestry.pop)

    if pop is not None and pop != 'all':
        mt = mt.filter_cols(mt.pop == pop)

    if filter_samples:
        mt = mt.filter_cols(mt.pass_qc)
    
    if filter_variants:
        mt = mt.filter_rows(mt.info.AC[0] > 0)

        if (analysis_type == 'gene') or (analysis_type == 'variant' and not use_array_for_variant):
            mt = annotate_adj(mt)
            mt = mt.filter_entries(mt.adj)
        
        mt = mt.filter_rows(hl.agg.any(mt.GT.n_alt_alleles() > 0))
    
    return mt


def generate_call_stats_ht(sample_qc, analysis_type, overwrite,
                           use_drc_ancestry_data=False, use_array_for_variant=False):
    # TODO add or remove VEP as an option
    pops = deepcopy(POPS)
    pops.append('all')
    for pop in pops:
        path = get_call_stats_ht_path(GENO_PATH, pop=pop, sample_qc=sample_qc, 
                                      analysis_type=analysis_type, 
                                      use_drc_ancestry_data=use_drc_ancestry_data,
                                      use_array_for_variant=use_array_for_variant)
        if not hl.hadoop_exists(os.path.join(path, '_SUCCESS')) or overwrite:
            mt = get_filtered_genotype_mt(analysis_type=analysis_type, pop=pop, 
                                          filter_variants=True, filter_samples=sample_qc,
                                          use_array_for_variant=use_array_for_variant,
                                          use_drc_ancestry_data=use_drc_ancestry_data)
            call_stats_ht = mt.annotate_rows(call_stats=hl.agg.call_stats(mt.GT, mt.alleles)).rows()
            call_stats_ht = call_stats_ht.naive_coalesce(1000).checkpoint(path, overwrite=overwrite)
    
    # if with_vep:
    #     full_path = get_call_stats_ht_path(GENO_PATH, pop='full', sample_qc=sample_qc, analysis_type=analysis_type)
    #     if not hl.hadoop_exists(f'{full_path}/_SUCCESS') or overwrite:            
    #         if not hl.hadoop_exists(f'{get_call_stats_ht_path(GENO_PATH, pop="all", sample_qc=(not sample_qc), analysis_type=analysis_type)}/_SUCCESS'):
    #             _ = generate_call_stats_ht(~sample_qc, analysis_type, overwrite, with_vep=False)
            
    #         ht = hl.read_table(get_aou_util_path(name="vep_full"))
    #         global_ht = hl.read_table(get_call_stats_ht_path(GENO_PATH, pop='all', sample_qc=False, analysis_type=analysis_type))
    #         ht = ht.filter(hl.is_defined(global_ht[ht.key]))
    #         ht = ht.annotate(
    #             freq = hl.struct(**{f'{pop.upper()}': hl.read_table(get_call_stats_ht_path(GENO_PATH, pop=pop, sample_qc=True, analysis_type=analysis_type))[ht.key].call_stats for pop in pops}),
    #             freq_raw = hl.struct(**{f'{pop.upper()}': hl.read_table(get_call_stats_ht_path(GENO_PATH, pop=pop, sample_qc=False, analysis_type=analysis_type))[ht.key].call_stats for pop in pops})
    #         )

    #         ht.describe()
    #         ht.naive_coalesce(1000).checkpoint(full_path, overwrite=args.overwrite)


def get_call_stats_ht(pop, sample_qc, analysis_type, use_drc_ancestry_data, use_array_for_variant, overwrite):
    path = get_call_stats_ht_path(GENO_PATH, pop=pop, sample_qc=sample_qc, 
                                  analysis_type=analysis_type, 
                                  use_drc_ancestry_data=use_drc_ancestry_data, 
                                  use_array_for_variant=use_array_for_variant)
    if not hl.hadoop_exists(os.path.join(path, '_SUCCESS')) or overwrite:
        _ = generate_call_stats_ht(sample_qc=sample_qc, analysis_type=analysis_type, 
                                   use_drc_ancestry_data=use_drc_ancestry_data,
                                   use_array_for_variant=use_array_for_variant)
    ht = hl.read_table(path)
    return ht


def mac_category_case_builder(call_stats_ac_expr, call_stats_af_expr, min_maf_common_variants: float = 0.01):
    return (
        hl.case()
        .when(call_stats_ac_expr <= 5, call_stats_ac_expr)
        .when(call_stats_ac_expr <= 10, 10)
        .when(call_stats_ac_expr <= 20, 20)
        .when(call_stats_af_expr <= 0.001, 0.001)
        .when(call_stats_af_expr <= min_maf_common_variants, min_maf_common_variants)
        .default(0.99)
    )


def filter_variants_for_grm(pop, analysis_type, use_array_for_variant, sample_qc, 
                            use_drc_ancestry_data=False, overwrite=False,
                            n_common_variants_to_keep=50000, # 100000 for per pop
                            min_call_rate=CALLRATE_CUTOFF, min_maf_common_variants=0.01, 
                            variants_per_mac_category=2000, variants_per_maf_category=10000):
    
    ht_sites_path = get_sites_for_grm_path(GENO_PATH, extension='ht',
                                           pop=pop, analysis_type=analysis_type, sample_qc=sample_qc,
                                           use_drc_ancestry_data=use_drc_ancestry_data,
                                           use_array_for_variant=use_array_for_variant,
                                           ld_pruned=False,
                                           n_common=n_common_variants_to_keep, 
                                           n_maf=variants_per_maf_category,
                                           n_mac=variants_per_mac_category)
    
    if overwrite or not hl.hadoop_exists(os.path.join(ht_sites_path, '_SUCCESS')):
        print(f'Number of common variants to sample: {n_common_variants_to_keep}')
        n_samples = get_n_samples_vec(analysis_type, sample_qc, use_array_for_variant=use_array_for_variant,
                                    use_drc_ancestry_data=use_drc_ancestry_data)
        ht = get_call_stats_ht(pop=pop, sample_qc=sample_qc, analysis_type=analysis_type,
                               use_drc_ancestry_data=use_drc_ancestry_data, 
                               use_array_for_variant=use_array_for_variant)
        ht = ht.filter(
            (ht.locus.in_autosome())
            & (ht.call_stats.AN >= (n_samples[pop] * 2 * min_call_rate))
            & (ht.call_stats.AC[1] > 0)
        )

        ht = ht.annotate(
            mac_category=mac_category_case_builder(ht.call_stats.AC[1], ht.call_stats.AF[1], min_maf_common_variants)
        )

        # From: https://hail.zulipchat.com/#narrow/stream/123010-Hail-Query-0.2E2-support/topic/.E2.9C.94.20randomly.20sample.20table/near/388162012
        bins = ht.aggregate(hl.agg.collect_as_set(ht.mac_category))
        ac_bins = [bin for bin in bins if (bin >= 1) and (bin<=5)]
        af_bins = [bin for bin in bins if (bin < 0.99) or (bin > 5)]

        sampled_common_variants = ht.aggregate(
            hl.agg.filter(
                    ht.call_stats.AF[1] > min_maf_common_variants,
                    hl.agg._reservoir_sample(ht.key, n_common_variants_to_keep),
                ),
        )
        print('Finished sampling common variants...')
        common_variants = [variant for variant in sampled_common_variants]

        binned_variants_af = ht.aggregate(
            hl.agg.array_agg(
                lambda x: hl.agg.filter(
                    ht.mac_category == x,
                    hl.agg._reservoir_sample(ht.key, variants_per_maf_category),
                ),
                hl.literal(af_bins),
            ),
        )
        print('Finished sampling rare variants...')

        binned_variants_ac = ht.aggregate(
            hl.agg.array_agg(
                lambda x: hl.agg.filter(
                    ht.mac_category == x,
                    hl.agg._reservoir_sample(ht.key, variants_per_mac_category),
                ),
                hl.literal(ac_bins),
            )
        )
        print('Finished sampling ultra-rare variants...')

        binned_rare_variants = binned_variants_ac + binned_variants_af
        rare_variants = [variant for bin in binned_rare_variants for variant in bin]

        print(f"N rare variants sampled: {len(rare_variants)}")
        print(f"N common variants sampled: {len(common_variants)}")
        rare_ht = hl.Table.parallelize(rare_variants).key_by(*ht.key.keys())
        common_ht = hl.Table.parallelize(common_variants).key_by(*ht.key.keys())
        ht = rare_ht.union(common_ht)
        ht = ht.checkpoint(ht_sites_path)
    else:
        ht = hl.read_table(ht_sites_path)
    
    ht.describe()

    return ht


def filter_mt_for_grm(pop, analysis_type, use_array_for_variant, sample_qc, 
                      use_drc_ancestry_data=False, overwrite=False,
                      n_common_variants_to_keep=50000, # 100000 for per pop
                      min_call_rate=CALLRATE_CUTOFF, min_maf_common_variants=0.01, 
                      variants_per_mac_category=2000, variants_per_maf_category=10000):
    
    mt_sites_path = get_sites_for_grm_path(GENO_PATH, extension='mt',
                                           pop=pop, analysis_type=analysis_type, sample_qc=sample_qc,
                                           use_drc_ancestry_data=use_drc_ancestry_data,
                                           use_array_for_variant=use_array_for_variant,
                                           ld_pruned=False,
                                           n_common=n_common_variants_to_keep, 
                                           n_maf=variants_per_maf_category,
                                           n_mac=variants_per_mac_category)
    
    if overwrite or not hl.hadoop_exists(os.path.join(mt_sites_path, '_SUCCESS')):
        mt = get_filtered_genotype_mt(analysis_type=analysis_type, pop=pop, filter_samples=sample_qc, filter_variants=True,
                                      use_array_for_variant=use_array_for_variant, use_drc_ancestry_data=use_drc_ancestry_data)
        filtered_ht = filter_variants_for_grm(pop=pop, analysis_type=analysis_type,
                                              use_array_for_variant=use_array_for_variant,
                                              use_drc_ancestry_data=use_drc_ancestry_data,
                                              overwrite=overwrite, 
                                              n_common_variants_to_keep=n_common_variants_to_keep,
                                              min_call_rate=min_call_rate, 
                                              min_maf_common_variants=min_maf_common_variants,
                                              variants_per_mac_category=variants_per_mac_category,
                                              variants_per_maf_category=variants_per_maf_category)
        
        print(f'Number of variants sampled for {pop.upper()}: {filtered_ht.count()}')
        mt = mt.filter_rows(hl.is_defined(filtered_ht[mt.row_key]))

        print("Removing HLA and inversion...")
        # Common inversion taken from Table S4 of https://www.ncbi.nlm.nih.gov/pubmed/27472961
        # (converted to GRCh38 by: https://liftover.broadinstitute.org/#input=chr8%3A8055789-11980649&hg=hg19-to-hg38 )
        # Also removing HLA, from https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38
        mt = mt.filter_rows(
            ~hl.parse_locus_interval(
                "chr8:8198267-12123140", reference_genome="GRCh38"
            ).contains(mt.locus)
            & ~hl.parse_locus_interval(
                "chr6:28510120-33480577", reference_genome="GRCh38"
            ).contains(mt.locus)
        )

        mt = mt.naive_coalesce(1000).checkpoint(mt_sites_path)
    else:
        mt = hl.read_matrix_table(mt_sites_path)

    print(mt.count())
    return mt


def main():
    a = 'a'


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument()

    args = parser.parse_args()
    main(**args)