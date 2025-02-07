#!/usr/bin/env python3
# This file includes functions to process AoU genomes.
# All pathing to new files is provided by SaigeImporters.py.
# All pathing to AoU files is provided by paths.py.

# Need to determine sample QC and variant QC.
import argparse
import json
import hail as hl
import pandas as pd
import subprocess

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


def get_n_samples_per_pop_vec(analysis_type, sample_qc, use_array_for_variant, use_drc_pop=False, use_hail=False):
    vec_path = get_n_samples_per_pop_path(GENO_PATH, analysis_type=analysis_type, sample_qc=sample_qc, 
                                          use_array_for_variant=use_array_for_variant,
                                          use_drc_pop=use_drc_pop)
    if not hl.hadoop_exists(vec_path):
        mt = get_filtered_genotype_mt(analysis_type=analysis_type,
                                      pop='all',
                                      filter_samples=sample_qc,
                                      filter_variants=False,
                                      use_array_for_variant=use_array_for_variant,
                                      use_drc_pop=use_drc_pop)
        ht = mt.cols()
        ht_ct = ht.group_by(ht.pop).aggregate(N = hl.agg.count())
        df_ct = ht_ct.to_pandas()
        df_ct.to_csv(vec_path, sep='\t', index=False)
    else:
        if use_hail:
            df_ct = hl.import_table(vec_path, impute=True).to_pandas()
        else:
            df_ct = pd.read_csv(vec_path, sep='\t')
    
    return {x['pop']: x.N for _, x in df_ct.iterrows()}


def get_filtered_genotype_mt(analysis_type, pop,
                             filter_samples=True, filter_variants=True,
                             use_array_for_variant=False, use_custom_pcs=None,
                             use_drc_pop=True, remove_related=False):
    if analysis_type == 'gene':
        mt_path = EXOME_PATH
    elif analysis_type == 'variant':
        mt_path = ARRAY_PATH if use_array_for_variant else WGS_PATH
    else:
        raise ValueError('ERROR: analysis_type can only be gene or variant.')
    
    if (analysis_type == 'gene') or (analysis_type == 'variant' and not use_array_for_variant):
        mt = hl.read_matrix_table(mt_path)
        mt = mt.filter_entries(hl.is_missing(mt.FT) | (mt.FT == 'PASS'))
        mt = mt.drop('variant_qc')
    else:
        mt = hl.read_matrix_table(mt_path, _n_partitions=15000)
    
    meta_ht = get_all_demographics(use_drc_pop=use_drc_pop, use_custom_pcs=use_custom_pcs)
    mt = mt.annotate_cols(**meta_ht[mt.col_key])
    mt = mt.annotate_cols(pop = mt.ancestry.pop)

    if pop is not None and pop != 'all':
        mt = mt.filter_cols(mt.pop == pop)

    if filter_samples:
        mt = mt.filter_cols(mt.pass_qc)
        if remove_related:
            related_samples = hl.import_table(get_aou_related_samples(), key='sample_id')
            mt = mt.anti_join_cols(related_samples)

    if filter_variants:
        mt = mt.filter_rows(mt.info.AC[0] > 0)

        if (analysis_type == 'gene') or (analysis_type == 'variant' and not use_array_for_variant):
            mt = annotate_adj(mt)
            mt = mt.filter_entries(mt.adj)
        
        mt = mt.filter_rows(hl.agg.any(mt.GT.n_alt_alleles() > 0))
    
    return mt


def generate_call_stats_ht(sample_qc, analysis_type, overwrite, checkpoint=False,
                           use_drc_pop=False, use_array_for_variant=False):
    # TODO add or remove VEP as an option

    n_partitions = 1000 if analysis_type=='variant' and use_array_for_variant else 10000
    
    pops = deepcopy(POPS)
    pops.append('all')
    
    for pop in pops:
        path = get_call_stats_ht_path(GENO_PATH, pop=pop, sample_qc=sample_qc, 
                                      analysis_type=analysis_type, 
                                      use_drc_pop=use_drc_pop,
                                      use_array_for_variant=use_array_for_variant)
        if not hl.hadoop_exists(os.path.join(path, '_SUCCESS')) or overwrite:
            print(f'Generating call stats for {analysis_type}, pop {pop}...')
            mt = get_filtered_genotype_mt(analysis_type=analysis_type, pop=pop, 
                                          filter_variants=True, filter_samples=sample_qc,
                                          use_array_for_variant=use_array_for_variant,
                                          use_drc_pop=use_drc_pop)
            if checkpoint:
                # NOTE: this may be a giant MatrixTable!
                # NOTE: you probably never want to enable this!
                mt_staging_path = f'{TEMP_PATH}/tmp_mt_genotypes_{analysis_type}_for_callstats_{pop}.mt'
                mt = mt.checkpoint(mt_staging_path)

            call_stats_ht = mt.annotate_rows(call_stats=hl.agg.call_stats(mt.GT, mt.alleles)).rows()
            call_stats_ht = call_stats_ht.naive_coalesce(n_partitions).checkpoint(path, overwrite=overwrite)
    
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


def get_call_stats_ht(pop, sample_qc, analysis_type, use_drc_pop, use_array_for_variant, overwrite):
    path = get_call_stats_ht_path(GENO_PATH, pop=pop, sample_qc=sample_qc, 
                                  analysis_type=analysis_type, 
                                  use_drc_pop=use_drc_pop, 
                                  use_array_for_variant=use_array_for_variant)
    if not hl.hadoop_exists(os.path.join(path, '_SUCCESS')) or overwrite:
        _ = generate_call_stats_ht(sample_qc=sample_qc, analysis_type=analysis_type, 
                                   use_drc_pop=use_drc_pop,
                                   use_array_for_variant=use_array_for_variant,
                                   overwrite=overwrite)
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


def get_call_rate_filtered_variants(pop, analysis_type, sample_qc, use_array_for_variant, use_drc_pop, 
                                    min_call_rate=CALLRATE_CUTOFF, only_autosomes=False, overwrite=False, ac_filter_override=0,
                                    use_hail_nsamp=False):
    n_samples = get_n_samples_per_pop_vec(analysis_type, sample_qc, use_array_for_variant=use_array_for_variant,
                                          use_drc_pop=use_drc_pop, use_hail=use_hail_nsamp)
    ht = get_call_stats_ht(pop=pop, sample_qc=sample_qc, analysis_type=analysis_type,
                            use_drc_pop=use_drc_pop, 
                            use_array_for_variant=use_array_for_variant, overwrite=overwrite)
    ht = ht.filter(
        (ht.call_stats.AN >= (n_samples[pop] * 2 * min_call_rate))
        & (ht.call_stats.AC[1] > ac_filter_override)
    )

    if only_autosomes:
        ht = ht.filter(ht.locus.in_autosome())
    
    return ht


def filter_common_variants_for_null(pop, analysis_type, use_array_for_variant, sample_qc, 
                                    use_drc_pop=False, overwrite=False,
                                    n_common_variants_to_keep=50000, # 100000 for per pop
                                    min_call_rate=CALLRATE_CUTOFF, min_maf_common_variants=0.01):
    
    ht_sites_path = get_sites_for_null_path(GENO_PATH, extension='common.ht',
                                            pop=pop, analysis_type=analysis_type, sample_qc=sample_qc,
                                            use_drc_pop=use_drc_pop,
                                            use_array_for_variant=use_array_for_variant,
                                            ld_pruned=False,
                                            n_common=n_common_variants_to_keep, 
                                            n_maf=0,
                                            n_mac=0)
    
    if overwrite or not hl.hadoop_exists(os.path.join(ht_sites_path, '_SUCCESS')):
        print(f'Number of common variants to sample: {n_common_variants_to_keep}')
        ht = get_call_rate_filtered_variants(pop=pop, analysis_type=analysis_type, sample_qc=sample_qc,
                                             use_array_for_variant=use_array_for_variant, use_drc_pop=use_drc_pop,
                                             min_call_rate=min_call_rate, only_autosomes=True)

        sampled_common_variants = ht.aggregate(
            hl.agg.filter(
                    ht.call_stats.AF[1] > min_maf_common_variants,
                    hl.agg._reservoir_sample(ht.key, n_common_variants_to_keep),
                ),
        )
        print('Finished sampling common variants...')
        common_variants = [variant for variant in sampled_common_variants]
        print(f"N common variants sampled: {len(common_variants)}")
        
        common_ht = hl.Table.parallelize(common_variants).key_by(*ht.key.keys())
        ht = common_ht
        ht = ht.checkpoint(ht_sites_path)
    else:
        ht = hl.read_table(ht_sites_path)
    
    ht.describe()

    return ht


def filter_rare_variants_for_null(pop, analysis_type, use_array_for_variant, sample_qc, 
                                  use_drc_pop=False, overwrite=False,
                                  min_call_rate=CALLRATE_CUTOFF, min_maf_common_variants = 0.01,
                                  variants_per_mac_category=2000, variants_per_maf_category=10000):
    
    ht_sites_path = get_sites_for_null_path(GENO_PATH, extension='rare.ht',
                                            pop=pop, analysis_type=analysis_type, sample_qc=sample_qc,
                                            use_drc_pop=use_drc_pop,
                                            use_array_for_variant=use_array_for_variant,
                                            ld_pruned=False,
                                            n_common=0,
                                            n_maf=variants_per_maf_category,
                                            n_mac=variants_per_mac_category)
    
    if overwrite or not hl.hadoop_exists(os.path.join(ht_sites_path, '_SUCCESS')):

        ht = get_call_rate_filtered_variants(pop=pop, analysis_type=analysis_type, sample_qc=sample_qc,
                                             use_array_for_variant=use_array_for_variant, use_drc_pop=use_drc_pop,
                                             min_call_rate=min_call_rate, only_autosomes=True)

        ht = ht.annotate(
            mac_category=mac_category_case_builder(ht.call_stats.AC[1], ht.call_stats.AF[1], min_maf_common_variants)
        )

        # From: https://hail.zulipchat.com/#narrow/stream/123010-Hail-Query-0.2E2-support/topic/.E2.9C.94.20randomly.20sample.20table/near/388162012
        bins = ht.aggregate(hl.agg.collect_as_set(ht.mac_category))
        ac_bins = [bin for bin in bins if (bin >= 1) and (bin<=5)]
        af_bins = [bin for bin in bins if (bin < 0.99) or (bin > 5)]

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
        rare_ht = hl.Table.parallelize(rare_variants).key_by(*ht.key.keys())
        ht = rare_ht
        ht = ht.checkpoint(ht_sites_path)
    else:
        ht = hl.read_table(ht_sites_path)
    
    ht.describe()

    return ht


def generate_plink_files_for_null(pop, sample_qc, use_drc_pop, overwrite,
                                  n_common_variants_to_keep=50000, min_maf_common_variants=0.01,
                                  variants_per_mac_category=2000, variants_per_maf_category=10000,
                                  min_call_rate=CALLRATE_CUTOFF):
    
    # first produce filtered variants
    ht_common = filter_common_variants_for_null(pop=pop, analysis_type='variant',
                                                use_array_for_variant=False,
                                                sample_qc=sample_qc, 
                                                use_drc_pop=use_drc_pop,
                                                overwrite=overwrite, 
                                                n_common_variants_to_keep=n_common_variants_to_keep,
                                                min_call_rate=min_call_rate, 
                                                min_maf_common_variants=min_maf_common_variants)
    ht_rare = filter_rare_variants_for_null(pop=pop, analysis_type='gene',
                                            use_array_for_variant=False,
                                            sample_qc=sample_qc, 
                                            use_drc_pop=use_drc_pop,
                                            overwrite=overwrite, 
                                            min_call_rate=min_call_rate, 
                                            min_maf_common_variants=min_maf_common_variants,
                                            variants_per_mac_category=variants_per_mac_category,
                                            variants_per_maf_category=variants_per_maf_category)

    # then produce filtered mt
    mt_sites_path = get_sites_for_null_path(GENO_PATH, extension='mt',
                                            pop=pop, analysis_type='both',
                                            sample_qc=sample_qc,
                                            use_drc_pop=use_drc_pop,
                                            use_array_for_variant=False,
                                            ld_pruned=False,
                                            n_common=n_common_variants_to_keep,
                                            n_maf=variants_per_maf_category,
                                            n_mac=variants_per_mac_category)
    
    if overwrite or not hl.hadoop_exists(os.path.join(mt_sites_path, '_SUCCESS')):
        print(f'Exporting subsampled genotype MatrixTables for {pop}...')
        mt_variant = get_filtered_genotype_mt(analysis_type='variant', pop=pop, 
                                              filter_samples=sample_qc, 
                                              filter_variants=True,
                                              use_array_for_variant=False, 
                                              use_drc_pop=use_drc_pop)
        mt_variant = mt_variant.semi_join_rows(ht_common)
        mt_variant = mt_variant.naive_coalesce(30000).checkpoint(f'{TEMP_PATH}/tmp_mt_genotypes_wgs_subsample_{pop}.mt', overwrite=True)

        mt_gene = get_filtered_genotype_mt(analysis_type='gene', pop=pop, 
                                           filter_samples=sample_qc, 
                                           filter_variants=True,
                                           use_array_for_variant=False, 
                                           use_drc_pop=use_drc_pop)
        mt_gene = mt_gene.semi_join_rows(ht_rare)
        mt_gene = mt_gene.naive_coalesce(30000).checkpoint(f'{TEMP_PATH}/tmp_mt_genotypes_exome_subsample_{pop}.mt', overwrite=True)

        mt = mt_variant.union_rows(mt_gene)
        mt = mt.filter_rows(~hl.parse_locus_interval(INVERSION_LOCUS, reference_genome="GRCh38").contains(mt.locus) & \
                            ~hl.parse_locus_interval(HLA_LOCUS, reference_genome="GRCh38").contains(mt.locus))
        
        mt = mt.naive_coalesce(2500).checkpoint(mt_sites_path)
    else:
        mt = hl.read_matrix_table(mt_sites_path)
        
    # then perform LD pruning
    ht_pruned_path = get_sites_for_null_path(GENO_PATH, extension='ht',
                                             pop=pop, analysis_type='both',
                                             sample_qc=sample_qc,
                                             use_drc_pop=use_drc_pop,
                                             use_array_for_variant=False,
                                             ld_pruned=True,
                                             n_common=n_common_variants_to_keep,
                                             n_maf=variants_per_maf_category,
                                             n_mac=variants_per_mac_category)
    
    if overwrite or not hl.hadoop_exists(os.path.join(ht_pruned_path, '_SUCCESS')):
        print(f'Performing LD pruning for pop {pop} for subsampled variants...')
        mt = mt.unfilter_entries()
        ht_pruned_sites = hl.ld_prune(mt.GT,
                                      r2=0.1,
                                      bp_window_size=int(1e7),
                                      memory_per_core=1024)
        ht_pruned_sites = ht_pruned_sites.checkpoint(ht_pruned_path)
    else:
        ht_pruned_sites = hl.read_table(ht_pruned_path)

    # then save as PLINK files
    bed_path, _, _ = get_plink_for_null_path(GENO_PATH, pop=pop,
                                             sample_qc=sample_qc,
                                             analysis_type='both',
                                             ld_pruned=True,
                                             use_drc_pop=use_drc_pop,
                                             use_array_for_variant=False,
                                             n_common=n_common_variants_to_keep,
                                             n_maf=variants_per_maf_category,
                                             n_mac=variants_per_mac_category)
    if overwrite or not hl.hadoop_exists(bed_path):
        print(f'Exporting subsampled plink files for pop {pop}...')
        mt = hl.read_matrix_table(mt_sites_path, _n_partitions=250)
        mt_for_export = mt.semi_join_rows(ht_pruned_sites)
        hl.export_plink(mt_for_export, os.path.splitext(bed_path)[0])


def get_filtered_array_mt_for_pruning(pop, sample_qc, min_af=0.01,
                                      min_call_rate=CALLRATE_CUTOFF,
                                      use_drc_pop=False, overwrite=False):
    n_samples = get_n_samples_per_pop_vec(analysis_type='variant', sample_qc=sample_qc, 
                                            use_array_for_variant=True,
                                            use_drc_pop=use_drc_pop)
    ht = get_call_stats_ht(pop=pop, sample_qc=sample_qc, analysis_type='variant',
                            use_drc_pop=use_drc_pop, 
                            use_array_for_variant=True,
                            overwrite=overwrite)
    
    # prior to LD pruning, ensure that autosomes only are found
    ht = ht.filter(
        (ht.locus.in_autosome())
        & (ht.call_stats.AN >= (n_samples[pop] * 2 * min_call_rate))
        & (ht.call_stats.AC[1] > 0)
        & (ht.call_stats.AC[0] > 0)
    )

    ht = ht.filter(~hl.parse_locus_interval(INVERSION_LOCUS, reference_genome="GRCh38").contains(ht.locus) & \
                    ~hl.parse_locus_interval(HLA_LOCUS, reference_genome="GRCh38").contains(ht.locus))

    # now filtering based on MAF
    ht = ht.annotate(maf = hl.min(ht.call_stats.AF))
    ht = ht.filter(ht.maf >= min_af)
    ht = ht.checkpoint(f'{TEMP_PATH}/tmp_ht_call_stats_for_pruning_{pop}.ht', overwrite=True)

    geno_mt = get_filtered_genotype_mt(analysis_type='variant', pop=pop, filter_samples=sample_qc, filter_variants=True,
                                        use_array_for_variant=True, use_drc_pop=use_drc_pop)
    geno_mt = geno_mt.semi_join_rows(ht)
    geno_mt = geno_mt.unfilter_entries()
    return geno_mt


def plink_ld_pruned_mt(sample_qc, saige_importers_path, wdl_path, min_af=0.01,
                       min_call_rate=CALLRATE_CUTOFF, 
                       use_drc_pop=False, 
                       overwrite=False):

    from cromwell.classes import CromwellManager

    baseline = {'ld_prune.step_size': 1,
                'ld_prune.window_size': 10000,
                'ld_prune.r2': 0.1,
                'ld_prune.SaigeImporters': saige_importers_path,
                'ld_prune.n_cpu': 2}

    plink_roots = []
    output_roots = []

    for pop in POPS:
        ld_pruned_ht_path = get_ld_pruned_array_data_path(GENO_PATH, pop=pop, sample_qc=sample_qc,
                                                          use_drc_pop=use_drc_pop,
                                                          af_cutoff=min_af, extension='ht', 
                                                          window='1e7',
                                                          use_plink=True)

        if overwrite or not hl.hadoop_exists(os.path.join(ld_pruned_ht_path, '_SUCCESS')):

            geno_mt = get_filtered_array_mt_for_pruning(pop=pop, sample_qc=sample_qc, min_af=min_af,
                                                        min_call_rate=min_call_rate, 
                                                        use_drc_pop=use_drc_pop, 
                                                        overwrite=overwrite)
            
            for chr in AUTOSOMES:
                plink_path = get_plink_inputs_ld_prune(GENO_PATH, pop=pop, chr=chr, sample_qc=sample_qc,
                                                       use_drc_pop=use_drc_pop,
                                                       af_cutoff=min_af, pruned=None, extension='bed')
                plink_root = os.path.splitext(plink_path)[0]

                plink_out = get_plink_inputs_ld_prune(GENO_PATH, pop=pop, chr=chr, sample_qc=sample_qc,
                                                     use_drc_pop=use_drc_pop,
                                                     af_cutoff=min_af, pruned='1e7', extension='txt')
                output_root = os.path.splitext(plink_out)[0]
                
                if overwrite or not hl.hadoop_exists(plink_out):
                    # output file from pruning not found, so add to pruning queue
                    plink_roots.append(plink_root)
                    output_roots.append(output_root)

                    if overwrite or not hl.hadoop_exists(plink_path):
                        print(f'Exporting plink files for pruning for {chr} in pop {pop}...')

                        if overwrite or not hl.hadoop_exists(plink_path):
                            this_chr_geno = geno_mt.filter_rows(geno_mt.locus.contig == chr)
                            hl.export_plink(this_chr_geno, plink_root)

    this_run = {'bedfile': [x + '.bed' for x in plink_roots],
                'bimfile': [x + '.bim' for x in plink_roots],
                'famfile': [x + '.fam' for x in plink_roots],
                'gs_output_path_root': output_roots}
    df = pd.DataFrame(this_run)
    
    if df.shape[0] > 0:
        with open(os.path.abspath('./saige_template.json'), 'w') as j:
            json.dump(baseline, j)

        # run LD pruning using Cromwell
        manager = CromwellManager(run_name='ld_prune',
                                  inputs_file=df,
                                  json_template_path=os.path.abspath('./saige_template.json'),
                                  wdl_path=wdl_path,
                                  batch=None, limit=199, n_parallel_workflows=199, 
                                  add_requester_pays_parameter=False,
                                  restart=False, batches_precomputed=False, 
                                  submission_sleep=0, check_freq=60, quiet=True)
        manager.run_pipeline(submission_retries=0, cromwell_timeout=60, skip_waiting=False)

    mt_dict = {}
    for pop in POPS:

        ld_pruned_ht_path = get_ld_pruned_array_data_path(GENO_PATH, pop=pop, sample_qc=sample_qc,
                                                          use_drc_pop=use_drc_pop,
                                                          af_cutoff=min_af, extension='ht', 
                                                          window='1e7',
                                                          use_plink=True)

        if overwrite or not hl.hadoop_exists(os.path.join(ld_pruned_ht_path, '_SUCCESS')):

            ht_list = []

            for chr in AUTOSOMES:
                plink_out = get_plink_inputs_ld_prune(GENO_PATH, pop=pop, chr=chr, sample_qc=sample_qc,
                                                      use_drc_pop=use_drc_pop,
                                                      af_cutoff=min_af, pruned='1e7', extension='txt')
                ht = hl.import_table(plink_out, no_header=True, types={'f0':'str'})
                ht = ht.annotate(locus = hl.locus(contig = ht.f0.split(':')[0],
                                                    pos = hl.int(ht.f0.split(':')[1]),
                                                    reference_genome='GRCh38'),
                                    alleles = ht.f0.split(':')[2:4])
                ht = ht.key_by('locus','alleles').select()
                ht_list.append(ht)
            
            ht_prune = hl.Table.union(*ht_list).repartition(50).checkpoint(ld_pruned_ht_path)
        else:
            ht_prune = hl.read_table(ld_pruned_ht_path)
        
        print(f'After pruning using PLINK for {pop}, there are {str(ht_prune.count())} variants.')

        mt = get_filtered_genotype_mt(analysis_type='variant', pop=pop, filter_samples=sample_qc, filter_variants=True,
                                      use_array_for_variant=True, use_drc_pop=use_drc_pop)
        mt = mt.semi_join_rows(ht_prune)
        mt_dict.update({pop: mt})
    
    return mt_dict


def hail_ld_pruned_mt(pop, sample_qc, min_af=0.01,
                      min_call_rate=CALLRATE_CUTOFF, 
                      use_drc_pop=False, 
                      overwrite=False):
    ld_pruned_ht_path = get_ld_pruned_array_data_path(GENO_PATH, pop=pop, sample_qc=sample_qc,
                                                      use_drc_pop=use_drc_pop,
                                                      af_cutoff=min_af, window='1e7',
                                                      extension='ht', use_plink=False)
    
    if overwrite or not hl.hadoop_exists(os.path.join(ld_pruned_ht_path, '_SUCCESS')):
        geno_mt = get_filtered_array_mt_for_pruning(pop=pop, sample_qc=sample_qc, min_af=min_af,
                                                    min_call_rate=min_call_rate, 
                                                    use_drc_pop=use_drc_pop, 
                                                    overwrite=overwrite)
        mt_staging_path = f'{TEMP_PATH}/tmp_mt_genotypes_for_pruning_{pop}.mt'
        geno_mt.write(mt_staging_path, overwrite=True)
        geno_mt = hl.read_matrix_table(mt_staging_path)

        ht_prune = hl.ld_prune(geno_mt.GT,
                               r2=0.1,
                               bp_window_size=int(1e7),
                               memory_per_core=2048)
        ht_prune = ht_prune.naive_coalesce(1000).checkpoint(ld_pruned_ht_path, overwrite=True)
    else:
        ht_prune = hl.read_table(ld_pruned_ht_path)
    
    print(f'After pruning using Hail for pop {pop}, there are {str(ht_prune.count())} variants.')
    mt = get_filtered_genotype_mt(analysis_type='variant', pop=pop, filter_samples=sample_qc, filter_variants=True,
                                  use_array_for_variant=True, use_drc_pop=use_drc_pop)
    mt = mt.semi_join_rows(ht_prune)
    return mt


def generate_plink_files_for_grm(mt_dict, sample_qc, min_af=0.01,
                                 use_drc_pop=False, 
                                 overwrite=False, use_plink=False):
    for pop, mt in mt_dict.items():
        plink_path = get_ld_pruned_array_data_path(GENO_PATH, pop=pop, sample_qc=sample_qc,
                                                   use_drc_pop=use_drc_pop,
                                                   af_cutoff=min_af, extension='bed', 
                                                   window='1e7',
                                                   use_plink=use_plink)
        if overwrite or not hl.hadoop_exists(plink_path):
            mt_staging_path = f'{TEMP_PATH}/tmp_mt_genotypes_for_pruning_prior_to_plinkexport_{pop}.mt'
            mt.write(mt_staging_path, overwrite=True)
            mt = hl.read_matrix_table(mt_staging_path, _n_partitions=250)
            hl.export_plink(mt, os.path.splitext(plink_path)[0])
            print(f'PLINK files generated for sparse GRM construction for {pop}.')
        else:
            print(f'PLINK files already generated for sparse GRM construction for {pop}.')


def generate_sparse_grm_distributed(pops, sample_qc, af_cutoff,
                                    saige_importers_path,
                                    wdl_path,
                                    n_markers=2000,
                                    relatedness=0.125,
                                    no_wait=False,
                                    n_cpu_sparse=32,
                                    use_drc_pop=False,
                                    overwrite=False,
                                    use_plink=True,
                                    use_array_data=True,
                                    n_common=50000,
                                    n_maf=10000,
                                    n_mac=2000):
    """
    Uses cromwell to distribute across pops.
    """

    from cromwell.classes import CromwellManager

    pops_to_queue = []

    for pop in pops:
        if use_array_data:
            mtx, _ = get_sparse_grm_path(GENO_PATH, pop=pop,
                                         n_markers=n_markers,
                                         relatedness=relatedness,
                                         sample_qc=sample_qc,
                                         af_cutoff=af_cutoff,
                                         use_drc_pop=use_drc_pop,
                                         use_plink=use_plink,
                                         use_array_data='')
        else:
            mtx, _ = get_sparse_grm_path(GENO_PATH, pop=pop,
                                         n_markers=n_markers,
                                         relatedness=relatedness,
                                         sample_qc=sample_qc,
                                         af_cutoff=af_cutoff,
                                         use_drc_pop=use_drc_pop,
                                         use_plink=use_plink,
                                         use_array_data=f'{str(n_common)}_{str(n_maf)}_{str(n_mac)}')
        if overwrite or not hl.hadoop_exists(mtx):
            pops_to_queue.append(pop)

    if len(pops_to_queue) > 0:
        # make json
        def remove_bucket(path):
            return re.sub('^'+BUCKET, '', path)

        this_run = {'pop': pops_to_queue}
        df = pd.DataFrame(this_run)
        df.to_csv(os.path.abspath('./this_run.tsv'), index=False, sep='\t') # ADD THIS

        baseline = {'saige_sparse_grm.relatednessCutoff': relatedness,
                    'saige_sparse_grm.min_af': af_cutoff,
                    'saige_sparse_grm.num_markers': n_markers,
                    'saige_sparse_grm.use_array': use_array_data,
                    'saige_sparse_grm.n_common': n_common,
                    'saige_sparse_grm.n_maf': n_maf,
                    'saige_sparse_grm.n_mac': n_mac,
                    'saige_sparse_grm.gs_bucket': BUCKET,
                    'saige_sparse_grm.gs_genotype_path': remove_bucket(GENO_PATH),
                    'saige_sparse_grm.SaigeImporters': saige_importers_path,
                    'saige_sparse_grm.use_drc_pop': use_drc_pop,
                    'saige_sparse_grm.use_plink': use_plink,
                    'saige_sparse_grm.sample_qc': sample_qc,
                    'saige_sparse_grm.n_cpu': n_cpu_sparse}
        with open(os.path.abspath('./saige_template.json'), 'w') as j:
            json.dump(baseline, j)

        # run sparse GRM analysis
        print('NOW COMMENCING GENERATION OF SPARSE GRMs.')
        print('This stage will use Cromwell.')
        manager = CromwellManager(run_name='saige_sparse_grm_multipop_aou',
                                  inputs_file=df,
                                  json_template_path=os.path.abspath('./saige_template.json'),
                                  wdl_path=wdl_path,
                                  batch=len(pops), limit=len(pops)+1, n_parallel_workflows=len(pops)+1, 
                                  add_requester_pays_parameter=False,
                                  restart=False, batches_precomputed=False, continue_when_possible=True,
                                  submission_sleep=0, check_freq=120, quiet=False)
        manager.run_pipeline(submission_retries=0, cromwell_timeout=60, skip_waiting=no_wait)
    
    print('Generation of all sparse GRMs completed.')


def get_variant_intervals(pop, overwrite):
    this_chunk_size = int(CHUNK_SIZE[pop])
    variant_interval_path = get_saige_interval_path(GENO_PATH, analysis_type='variant', pop=pop)
    if not hl.hadoop_exists(variant_interval_path) or overwrite:
        print(f'Generating interval file for SAIGE (chunk size: {this_chunk_size})...')
        chr = []
        start = []
        end = []
        for chrom in CHROMOSOMES:
            CHROMOSOME_LEN = hl.get_reference('GRCh38').lengths
            chrom_length = CHROMOSOME_LEN[chrom]
            for start_pos in range(1, chrom_length, this_chunk_size):
                end_pos = (
                    chrom_length
                    if start_pos + this_chunk_size > chrom_length
                    else (start_pos + this_chunk_size)
                )
                chr.append(chrom)
                start.append(start_pos)
                end.append(end_pos)

        df = pd.DataFrame({'chrom': chr, 'start': start, 'end': end})
        df.to_csv(variant_interval_path, sep='\t', index=False)
    
    df = pd.read_csv(variant_interval_path, sep='\t')
    print(pop)
    print(df.head(5))
    print(df.shape)
    return(df)


def get_gene_based_intervals(pop, overwrite):

    gene_interval_path = get_saige_interval_path(GENO_PATH, analysis_type='gene', pop=pop)
    if not hl.hadoop_exists(gene_interval_path) or overwrite:
        # load gene boundaries
        gene_boundaries = generate_gene_boundries_via_vat(overwrite=False)
        gene_boundaries = gene_boundaries.annotate(contig = gene_boundaries.contig.find(lambda x: True))
        gene_boundaries = gene_boundaries.annotate(interval = hl.locus_interval(gene_boundaries.contig, gene_boundaries.min_pos, gene_boundaries.max_pos, reference_genome='GRCh38'))
        gene_boundaries = gene_boundaries.key_by('interval')

        ht_intervals = hl.Table.from_pandas(read_variant_intervals(GENO_PATH, pop, 'variant'))
        ht_intervals = ht_intervals.annotate(locus_start = hl.locus(ht_intervals.chrom, ht_intervals.start, reference_genome='GRCh38'),
                                            locus_end = hl.locus(ht_intervals.chrom, ht_intervals.end, reference_genome='GRCh38'))

        # find where there are interval starts that are overlapping genes
        ht_intervals_start = ht_intervals.rename({'locus_start': 'locus_start_orig'}).key_by('locus_start_orig')
        ht_intervals_start = ht_intervals_start.annotate(overlap_gene = gene_boundaries[ht_intervals_start.key])
        ht_intervals_start_mod = ht_intervals_start
        ht_intervals_start_mod.show()

        still_overlapping = True
        this_iter = 1
        print(f'Prior to loop, there were {str(ht_intervals_start_mod.aggregate(hl.agg.count_where(hl.is_defined(ht_intervals_start_mod.overlap_gene))))} overlapping start sites with genes.')
        while still_overlapping:
            ht_intervals_start_mod = ht_intervals_start_mod.persist()
            ht_intervals_start_mod = ht_intervals_start_mod.annotate(locus_start_pre = hl.locus(ht_intervals_start_mod.chrom, hl.if_else(hl.is_missing(ht_intervals_start_mod.overlap_gene), 
                                                                                                                                        ht_intervals_start_mod.key[0].position, 
                                                                                                                                        ht_intervals_start_mod.overlap_gene.min_pos-50), reference_genome='GRCh38'),
                                                                    locus_start_post = hl.locus(ht_intervals_start_mod.chrom, hl.if_else(hl.is_missing(ht_intervals_start_mod.overlap_gene), 
                                                                                                                                        ht_intervals_start_mod.key[0].position, 
                                                                                                                                        ht_intervals_start_mod.overlap_gene.max_pos+50), reference_genome='GRCh38'))
            ht_intervals_start_mod = ht_intervals_start_mod.key_by('locus_start_pre')
            ht_intervals_start_mod = ht_intervals_start_mod.annotate(overlap_gene_pre = gene_boundaries[ht_intervals_start_mod.key])
            ht_intervals_start_mod = ht_intervals_start_mod.key_by('locus_start_post')
            ht_intervals_start_mod = ht_intervals_start_mod.annotate(overlap_gene_post = gene_boundaries[ht_intervals_start_mod.key])
            ht_intervals_start_mod = ht_intervals_start_mod.annotate(locus_start = hl.if_else(hl.is_missing(ht_intervals_start_mod.overlap_gene_pre), 
                                                                                            ht_intervals_start_mod.locus_start_pre, 
                                                                                            hl.if_else(hl.is_missing(ht_intervals_start_mod.overlap_gene_post),
                                                                                                        ht_intervals_start_mod.locus_start_post,
                                                                                                        ht_intervals_start_mod.locus_start_pre)),
                                                                    overlap_gene = hl.if_else(hl.is_missing(ht_intervals_start_mod.overlap_gene_pre), 
                                                                                            ht_intervals_start_mod.overlap_gene_pre, 
                                                                                            hl.if_else(hl.is_missing(ht_intervals_start_mod.overlap_gene_post),
                                                                                                        ht_intervals_start_mod.overlap_gene_post,
                                                                                                        ht_intervals_start_mod.overlap_gene_pre)))
            ht_intervals_start_mod = ht_intervals_start_mod.key_by('locus_start')
            ht_intervals_start_mod.filter(~hl.is_missing(ht_intervals_start_mod.overlap_gene)).show()
            
            noverlap = ht_intervals_start_mod.aggregate(hl.agg.count_where(hl.is_defined(ht_intervals_start_mod.overlap_gene)))
            print(f'After loop {str(this_iter)}, there were {str(noverlap)} overlapping start sites with genes.')
            if noverlap == 0:
                still_overlapping = False
                print('Completed.')
            else:
                this_iter += 1
        
        # now produce a fixed interval table by shifting around interval ending
        ht_intervals_start_mod_for_end = ht_intervals_start_mod.key_by()
        ht_intervals_start_mod_for_end = ht_intervals_start_mod_for_end.select(chrom=ht_intervals_start_mod_for_end.chrom, 
                                                                            start=ht_intervals_start_mod_for_end.locus_start.position, 
                                                                            end=ht_intervals_start_mod_for_end.locus_end.position)
        df_orig = ht_intervals_start_mod_for_end.to_pandas()

        chr = []
        start = []
        end = []
        CHROMOSOME_LEN = hl.get_reference('GRCh38').lengths

        prev_chr = None
        prev_end = None
        for _, row in df_orig.iterrows():
            this_chr = row['chrom']
            chrom_length = CHROMOSOME_LEN[this_chr]
            this_start = max(1, row['start'])
            this_end = min(row['end'], chrom_length)
            if prev_chr is not None and row['chrom'] == prev_chr:
                if this_start != prev_end:
                    end[-1] = this_start

            chr.append(this_chr)
            start.append(this_start)
            end.append(this_end)

            prev_chr = this_chr
            prev_end = this_end
                
        df = pd.DataFrame({'chrom': chr, 'start': start, 'end': end})
        df.to_csv(gene_interval_path, sep='\t', index=False)
    
    df = pd.read_csv(gene_interval_path, sep='\t')
    print(pop)
    print(df.head(5))
    print(df.shape)
    return(df)


def create_bgen_split_intervals(pop, git_path, wdl_path, callrate_filter, min_ac, analysis_type='variant',
                                mean_impute_missing=True, use_drc_pop=True, encoding='additive', 
                                limit=5000, n_cpu=8):

    from cromwell.classes import CromwellManager

    subprocess.run(['tar', '-czf', './saige_wdl.tar.gz', git_path])
    repo_tarball = os.path.join(BUCKET,'scripts/saige_wdl.tar.gz')
    _ = subprocess.run(['gsutil','cp','./saige_wdl.tar.gz',repo_tarball])

    baseline = {'split_bgen_intervals.sample_qc': True,
                'split_bgen_intervals.variant_qc': True,
                'split_bgen_intervals.use_drc_pop': use_drc_pop,
                'split_bgen_intervals.mean_impute_missing': mean_impute_missing,
                'split_bgen_intervals.call_rate_filter': callrate_filter,
                'split_bgen_intervals.min_ac': min_ac,
                'split_bgen_intervals.analysis_type': analysis_type,
                'split_bgen_intervals.encoding': encoding,
                'split_bgen_intervals.repo_tarball': repo_tarball,
                'split_bgen_intervals.tar_folder_path': git_path.lstrip('/'),
                'split_bgen_intervals.n_cpu': n_cpu}

    if type(pop) == list:
        this_pop_list = ','.join(pop)

        dct = {'chr': [],
               'start': [],
               'end': [],
               'bgen_prefix': [],
               'pop': []}

        for this_pop in pop:
            if analysis_type == 'variant':
                interval_list = get_variant_intervals(pop=this_pop, overwrite=False)
            elif analysis_type == 'gene':
                interval_list = get_gene_based_intervals(pop=this_pop, overwrite=False)
            bgen_prefix = get_wildcard_path_intervals_bgen(GENO_PATH, pop=this_pop, use_drc_pop=use_drc_pop, encoding=encoding, analysis_type=analysis_type)
        
            for _, row in interval_list.iterrows():
                this_bgen_prefix = bgen_prefix.replace('@', row['chrom']).replace('#', str(row['start'])).replace('?', str(row['end']))
                if not hl.hadoop_exists(this_bgen_prefix + '.bgen') or not hl.hadoop_exists(this_bgen_prefix + '.bgen.bgi'):
                    dct['chr'].append(row['chrom'])
                    dct['start'].append(row['start'])
                    dct['end'].append(row['end'])
                    dct['bgen_prefix'].append(this_bgen_prefix)
                    dct['pop'].append(this_pop)

    else:
        this_pop_list = pop
        baseline.update({'split_bgen_intervals.pop': pop})
    
        if analysis_type == 'variant':
            interval_list = get_variant_intervals(pop=pop, overwrite=False)
        elif analysis_type == 'gene':
            interval_list = get_gene_based_intervals(pop=pop, overwrite=False)
        bgen_prefix = get_wildcard_path_intervals_bgen(GENO_PATH, pop=pop, use_drc_pop=use_drc_pop, encoding=encoding, analysis_type=analysis_type)
        
        dct = {'chr': [],
               'start': [],
               'end': [],
               'bgen_prefix': []}
        
        for _, row in interval_list.iterrows():
            this_bgen_prefix = bgen_prefix.replace('@', row['chrom']).replace('#', str(row['start'])).replace('?', str(row['end']))
            if not hl.hadoop_exists(this_bgen_prefix + '.bgen') or not hl.hadoop_exists(this_bgen_prefix + '.bgen.bgi'):
                dct['chr'].append(row['chrom'])
                dct['start'].append(row['start'])
                dct['end'].append(row['end'])
                dct['bgen_prefix'].append(this_bgen_prefix)


    df = pd.DataFrame(dct)
    df.to_csv(os.path.abspath(f'./this_{this_pop_list}_run.tsv'), index=False, sep='\t')    

    with open(os.path.abspath(f'./saige_template_{this_pop_list}.json'), 'w') as j:
        json.dump(baseline, j)

    # run sparse GRM analysis
    print('NOW COMMENCING GENERATION OF BGENs.')
    print('This stage will use Cromwell.')
    manager = CromwellManager(run_name=f'saige_aou_split_{analysis_type}_bgen_{this_pop_list}',
                              inputs_file=df,
                              json_template_path=os.path.abspath(f'./saige_template_{this_pop_list}.json'),
                              wdl_path=wdl_path,
                              limit=limit, n_parallel_workflows=limit, 
                              add_requester_pays_parameter=True,
                              restart=False, batches_precomputed=False, 
                              submission_sleep=0, check_freq=120, quiet=False)
    manager.run_pipeline(submission_retries=0, cromwell_timeout=60, skip_waiting=True)
    return manager


def gt_to_gp(mt, location: str = "GP"):
    return mt.annotate_entries(
        **{
            location: hl.or_missing(
                hl.is_defined(mt.GT),
                hl.map(
                    lambda i: hl.if_else(
                        mt.GT.unphased_diploid_gt_index() == i, 1.0, 0.0
                    ),
                    hl.range(0, hl.triangle(hl.len(mt.alleles))),
                ), # makes calls into [0/1, 0/1, 0/1] vector where [1, 0, 0] is hom ref, etc
            )
        }
    )


def impute_missing_gp(mt, location: str = "GP", mean_impute: bool = True):
    mt = mt.annotate_entries(_gp=mt[location])
    if mean_impute:
        mt = mt.annotate_rows(
            _mean_gp=hl.agg.array_agg(lambda x: hl.agg.mean(x), mt._gp)
        )
        gp_expr = mt._mean_gp
    else:
        gp_expr = [1.0, 0.0, 0.0]
    return mt.annotate_entries(**{location: hl.or_else(mt._gp, gp_expr)}).drop("_gp")


def generate_vat_ht(overwrite=False):
    # import giant table and save as ht
    ht = hl.import_table(os.path.join(AUX_PATH, 'vat/vat_complete_v7.1.bgz.tsv.gz'), force_bgz=True, min_partitions=20000)
    ht = ht.checkpoint(os.path.join(TEMP_PATH, 'vat_prelim.ht'), _read_if_exists=not overwrite).key_by('contig', 'position', 'ref_allele', 'alt_allele')
    ht = ht.filter(hl.is_defined(ht.gene_symbol))
    ht = ht.transmute(consequences = ht.consequence.split(', ')).drop('vid')
    ht = ht.transmute(transcript_consequences = hl.struct(consequence_terms = ht.consequences,
                                                        gene_symbol = ht.gene_symbol,
                                                        is_canonical_transcript = ht.is_canonical_transcript,
                                                        transcript = ht.transcript))
    ht = ht.drop(*[x for x in ht.row if re.search('^gvs', x)])
    ht = ht.drop(*[x for x in ht.row if re.search('^gnomad', x) and x != 'gnomad_failed_filter'])
    ht = ht.checkpoint(os.path.join(TEMP_PATH, 'vat_truly_final.ht'), _read_if_exists=not overwrite)

    # set up processing functions
    csqs = hl.literal(CSQ_ORDER)
    PENALIZE_LOFTEE = False
    csq_dict = hl.literal(dict(zip(CSQ_ORDER, range(len(CSQ_ORDER)))))

    def _get_most_severe_consequence_expr(
        csq_expr: hl.expr.ArrayExpression,
        csq_order = None,
    ) -> hl.expr.StringExpression:
        if csq_order is None:
            csq_order = CSQ_ORDER
        csqs = hl.literal(csq_order)
        return csqs.find(lambda c: csq_expr.contains(c))


    def _add_most_severe_consequence_to_consequence(
        tc,
        csq_order = None,
        most_severe_csq_field: str = "most_severe_consequence",
    ):
        csq = lambda x: _get_most_severe_consequence_expr(x.consequence_terms, csq_order)
        if isinstance(tc, hl.expr.StructExpression):
            return tc.annotate(**{most_severe_csq_field: csq(tc)})
        else:
            return tc.map(lambda x: x.annotate(**{most_severe_csq_field: csq(x)}))


    def _find_worst_transcript_consequence(
        tcl: hl.expr.ArrayExpression,
        include_loftee=False, include_polyphen=False
    ) -> hl.expr.StructExpression:
        """
        Find the worst transcript consequence in an array of transcript consequences.

        :param tcl: Array of transcript consequences.
        :return: Worst transcript consequence.
        """
        flag = 500
        no_flag = flag * (1 + PENALIZE_LOFTEE)

        # Score each consequence based on the order in csq_order.
        score_expr = tcl.map(
            lambda tc: csq_dict[csqs.find(lambda x: x == tc.most_severe_consequence)]
        )

        if include_loftee:
            # Determine the score adjustment based on the consequence's LOF and LOF flags.
            sub_expr = tcl.map(
                lambda tc: (
                    hl.case(missing_false=True)
                    .when((tc.lof == "HC") & hl.or_else(tc.lof_flags == "", True), no_flag)
                    .when((tc.lof == "HC") & (tc.lof_flags != ""), flag)
                    .when(tc.lof == "OS", 20)
                    .when(tc.lof == "LC", 10)
                    .default(0)
                )
            )

        if include_polyphen:
            # If requested, determine the score adjustment based on the consequence's
            # PolyPhen prediction.
                polyphen_sub_expr = tcl.map(
                    lambda tc: (
                        hl.case(missing_false=True)
                        .when(tc.polyphen_prediction == "probably_damaging", 0.5)
                        .when(tc.polyphen_prediction == "possibly_damaging", 0.25)
                        .when(tc.polyphen_prediction == "benign", 0.1)
                        .default(0)
                    )
                )
                sub_expr = hl.map(lambda s, ps: s + ps, sub_expr, polyphen_sub_expr)

        # Calculate the final consequence score.
        if include_loftee:
            tcl = hl.map(
                lambda tc, s, ss: tc.annotate(csq_score=s - ss), tcl, score_expr, sub_expr
            )
        else:
            tcl = hl.map(
                lambda tc, s: tc.annotate(csq_score=s), tcl, score_expr
            )

        # Return the worst consequence based on the calculated score.
        return hl.or_missing(hl.len(tcl) > 0, hl.sorted(tcl, lambda x: x.csq_score)[0])


    # process VAT table
    # Annotate each transcript consequence with the 'most_severe_consequence'.
    ht_f = ht.filter(ht.transcript_consequences.gene_symbol != '')
    ht_f = ht_f.collect_by_key()
    ht_f = ht_f.annotate(vep = hl.struct())

    transcript_csqs = _add_most_severe_consequence_to_consequence(
        ht_f.values.transcript_consequences, CSQ_ORDER
    )

    # Group transcript consequences by gene and find the worst consequence for each.
    gene_dict = transcript_csqs.group_by(lambda tc: tc.gene_symbol)
    worst_csq_gene = gene_dict.map_values(_find_worst_transcript_consequence).values()
    sorted_scores = hl.sorted(worst_csq_gene, key=lambda tc: tc.csq_score)

    # Filter transcript consequences to only include canonical transcripts and find the
    # worst consequence for each gene.
    canonical = transcript_csqs.filter(lambda csq: csq.is_canonical_transcript == 'true')
    gene_canonical_dict = canonical.group_by(lambda tc: tc.gene_symbol)
    worst_csq_gene_canonical = gene_canonical_dict.map_values(
        _find_worst_transcript_consequence
    ).values()
    sorted_canonical_scores = hl.sorted(
        worst_csq_gene_canonical, key=lambda tc: tc.csq_score
    )

    vep_data = ht_f['vep'].annotate(
        transcript_consequences=transcript_csqs,
        worst_consequence_term=csqs.find(
            lambda c: transcript_csqs.map(
                lambda csq: csq.most_severe_consequence
            ).contains(c)
        ),
        worst_csq_by_gene=sorted_scores,
        worst_csq_by_gene_canonical=sorted_canonical_scores
    )
    ht_vep = ht_f.annotate(vep = vep_data)
    ht_vep = ht_vep.checkpoint(os.path.join(TEMP_PATH, 'vat_post_processed_vep_style_pre.ht'), _read_if_exists=True)
    ht_vep = ht_vep.annotate(vep = ht_vep.vep.annotate(
        worst_csq_for_variant=hl.or_missing(
            hl.len(ht_vep.vep.worst_csq_by_gene) > 0, ht_vep.vep.worst_csq_by_gene[0]
        ),
        worst_csq_for_variant_canonical=hl.or_missing(
            hl.len(ht_vep.vep.worst_csq_by_gene_canonical) > 0, ht_vep.vep.worst_csq_by_gene_canonical[0]
        ),
    ))
    ht_vep = ht_vep.checkpoint(get_processed_vat_path(ANNOT_PATH), _read_if_exists=True)
    
    return ht_vep


def generate_gene_boundries_via_vat(remove_cross_contig_genes=True, overwrite=False):
    ht = generate_vat_ht(overwrite=overwrite)
    ht_vep_e = ht.explode('values')
    ht_gene_position = ht_vep_e.group_by(ht_vep_e.values.transcript_consequences.gene_symbol
                            ).aggregate(contig = hl.agg.collect_as_set(ht_vep_e.contig),
                                        min_pos = hl.agg.min(hl.int(ht_vep_e.position)),
                                        max_pos = hl.agg.max(hl.int(ht_vep_e.position)))
    ht_gene_position = ht_gene_position.naive_coalesce(500).checkpoint(os.path.join(ANNOT_PATH, 'gene_boundry_using_vat_variants.ht'), _read_if_exists=True)

    if remove_cross_contig_genes:
        ht_gene_position = ht_gene_position.filter(hl.len(ht_gene_position.contig) == 1)
        ht_gene_position = ht_gene_position.filter(~ht_gene_position.gene_symbol.startswith('SNORA') & \
                                                   ~ht_gene_position.gene_symbol.startswith('SCARNA') & \
                                                    ~ht_gene_position.gene_symbol.startswith('SNORD') & \
                                                        (ht_gene_position.gene_symbol != 'U6atac'))
        
        # filter using UCSC table
        ucsc_table = get_ucsc_gene_table(BUCKET, ANNOT_PATH, overwrite=overwrite)
        annot_table = ht_gene_position.annotate(ucsc_table = ucsc_table[ht_gene_position.gene_symbol])
        annot_table = annot_table.annotate(size = annot_table.max_pos - annot_table.min_pos)
        annot_table_multi = annot_table.filter(annot_table.ucsc_table.has_multiple)
        annot_table_multi = annot_table_multi.annotate(size_ratio = (annot_table_multi.size - annot_table_multi.ucsc_table.size)/annot_table_multi.ucsc_table.size)
        genes_too_long = annot_table_multi.aggregate(hl.agg.filter((annot_table_multi.size_ratio >= 2.5) & (annot_table_multi.size >= 30000), hl.agg.collect(annot_table_multi.gene_symbol)))

        annot_table_single = annot_table.filter(~annot_table.ucsc_table.has_multiple)
        annot_table_single = annot_table_single.annotate(delta = (annot_table_single.size - annot_table_single.ucsc_table.size))
        annot_table_single = annot_table_single.filter(annot_table_single.size > 30000)
        annot_table_single = annot_table_single.annotate(size_ratio = annot_table_single.delta / annot_table_single.ucsc_table.size)
        annot_table_single = annot_table_single.filter((annot_table_single.size_ratio > 30))
        genes_too_long += annot_table_single.aggregate(hl.agg.collect(annot_table_single.gene_symbol))
        
        ht_gene_position = ht_gene_position.filter(~hl.literal(genes_too_long).contains(ht_gene_position.gene_symbol))
        return ht_gene_position
    else: 
        return ht_gene_position


def get_overlapping_genes(pop, analysis_type, overwrite=False):
    ht_gene_position = generate_gene_boundries_via_vat(overwrite=False)

    ht_gene_position_annot = ht_gene_position.repartition(20000)
    ht_gene_position_annot = ht_gene_position_annot.annotate(all_positions = hl.range(ht_gene_position_annot.min_pos, ht_gene_position_annot.max_pos))
    ht_gene_position_annot = ht_gene_position_annot.annotate(contig = ht_gene_position_annot.contig.find(lambda x: True))
    ht_gene_position_annot = ht_gene_position_annot.explode('all_positions')
    ht_gene_position_annot = ht_gene_position_annot.annotate(locus = hl.locus(ht_gene_position_annot.contig, ht_gene_position_annot.all_positions, reference_genome='GRCh38'))
    ht_gene_position_annot = ht_gene_position_annot.key_by('locus')
    ht_gene_position_annot = ht_gene_position_annot.checkpoint(os.path.join(TEMP_PATH, 'ht_gene_position.ht'), _read_if_exists=True)

    ht_intervals = hl.Table.from_pandas(read_variant_intervals(GENO_PATH, pop, analysis_type))
    ht_intervals = ht_intervals.annotate(interval = hl.locus_interval(ht_intervals.chrom, ht_intervals.start, ht_intervals.end, reference_genome='GRCh38'))
    ht_intervals = ht_intervals.annotate(interval_id = ht_intervals.chrom + ':' + hl.str(ht_intervals.start) + '-' + hl.str(ht_intervals.end)).key_by('interval')

    ht_gene_position_annot_interval = ht_gene_position_annot.annotate(interval_id = ht_intervals[ht_gene_position_annot.key].interval_id)
    intervals_by_gene = ht_gene_position_annot_interval.group_by(ht_gene_position_annot_interval.gene_symbol
                                                      ).aggregate(gene_start = hl.agg.collect_as_set(ht_gene_position_annot_interval.min_pos),
                                                                  gene_end = hl.agg.collect_as_set(ht_gene_position_annot_interval.max_pos),
                                                                  overlapping_intervals = hl.agg.collect_as_set(ht_gene_position_annot_interval.interval_id)
                                                      ).naive_coalesce(500)
    intervals_by_gene = intervals_by_gene.checkpoint(os.path.join(ANNOT_PATH, f'gene_to_{analysis_type}_interval_map_{pop}.ht'), overwrite=overwrite, _read_if_exists=not overwrite)
    return intervals_by_gene, ht_gene_position_annot


def generate_gene_group_files(pop, overwrite=False, use_canonical=False, min_cr=CALLRATE_CUTOFF):
    ht = generate_vat_ht(overwrite=False)
    ht = ht.annotate(var_id = ht.contig.replace('chr','') + ':' + ht.position + ':' + ht.ref_allele + ':' + ht.alt_allele)

    # load call rate table and remove variants with CR < 90%
    print(f'Prior to CR filtation, we have {str(ht.count())} records in pop {pop}.')
    ht_call_rate = get_call_rate_filtered_variants(pop=pop, analysis_type='gene',
                                                   sample_qc=True, use_array_for_variant=False,
                                                   use_drc_pop=True, min_call_rate=min_cr)
    ht_call_rate = ht_call_rate.key_by(chrom = ht_call_rate.locus.contig, 
                                       pos = hl.str(ht_call_rate.locus.position), 
                                       ref = ht_call_rate.alleles[0], 
                                       alt = ht_call_rate.alleles[1])
    ht = ht.semi_join(ht_call_rate)
    print(f'After CR filtation, we have {str(ht.count())} records in pop {pop}.')

    if use_canonical:
        ht = ht.annotate(csq_struct = ht.vep.worst_csq_by_gene_canonical)
    else:
        ht = ht.annotate(csq_struct = ht.vep.worst_csq_by_gene)
    hte = ht.explode('csq_struct')
    all_csq = PLOF_CSQS + MISSENSE_CSQS + SYNONYMOUS_CSQS
    lookup_csq = {}
    lookup_csq.update({x: 'lof' for x in PLOF_CSQS})
    lookup_csq.update({x: 'missense' for x in MISSENSE_CSQS})
    lookup_csq.update({x: 'synonymous' for x in SYNONYMOUS_CSQS})
    hte_f = hte.filter(hl.literal(all_csq).contains(hte.csq_struct.most_severe_consequence))
    hte_f = hte_f.annotate(csq_class = hl.literal(lookup_csq)[hte_f.csq_struct.most_severe_consequence])
    htgene = hte_f.group_by(hte_f.csq_struct.gene_symbol
                ).aggregate(contig = hl.agg.take(hte_f.contig, 1)[0], 
                            joint = hl.agg.collect([hte_f.var_id, hte_f.csq_class]))
    htgene = htgene.annotate(var = hl.str(' ').join(hl.map(lambda x: x[0], htgene.joint)),
                            anno = hl.str(' ').join(hl.map(lambda x: x[1], htgene.joint))).drop('joint')
    htgene = htgene.annotate(jt = ['var ' + htgene.var, 'anno ' + htgene.anno]).explode('jt').key_by()
    htgene = htgene.transmute(final_column = htgene.gene_symbol + ' ' + htgene.jt).drop('var', 'anno')

    for chrom in CHROMOSOMES:
        annotation_path = get_gene_annotation_path(ANNOT_PATH, chrom, pop)
        if not hl.hadoop_exists(annotation_path) or overwrite:
            htgene_this = htgene.filter(htgene.contig == chrom).drop('contig')
            htgene_this.export(annotation_path, header=False)
    
    return None


def main():
    sample_qc = True
    min_af = 0.01
    use_drc_pop=True
    overwrite=False
    no_wait=False

    n_markers = 2000
    relatedness = 0.125

    saige_importers_path = os.path.join(BUCKET,'scripts/SaigeImporters.py')
    plink_wdl_path = ''
    sparse_wdl_path = ''
    use_plink = True

    n_common = 50000
    min_maf_common = 0.01
    n_mac = 2000
    n_maf = 10000

    use_array_data_for_sparse = True
    
    if use_plink:
        mt_dict = plink_ld_pruned_mt(sample_qc=sample_qc,
                                     saige_importers_path=saige_importers_path,
                                     wdl_path=plink_wdl_path,
                                     min_af=min_af,
                                     use_drc_pop=use_drc_pop,
                                     overwrite=overwrite)

    else:
        mt_dict = {}
        for pop in POPS:
            # export GRM sample list
            mt = hail_ld_pruned_mt(pop=pop,
                                   sample_qc=sample_qc,
                                   min_af=min_af,
                                   use_drc_pop=use_drc_pop,
                                   overwrite=overwrite)
            mt_dict.update({pop: mt})

    for pop in POPS:
        with hl.hadoop_open(get_aou_samples_file_path(GENO_PATH, pop, sample_qc=sample_qc, use_plink=use_plink,
                                                      use_drc_pop=use_drc_pop), 'w') as f:
            f.write('\n'.join(mt.s.collect()) + '\n')

        # produce downsampled plink files
        # note this will include a bunch of markers from the rare spectrum using WGS for common variants and WES for rare variants
        _ = generate_plink_files_for_null(pop, sample_qc=sample_qc,
                                          use_drc_pop=use_drc_pop,
                                          overwrite=overwrite,
                                          n_common_variants_to_keep=n_common,
                                          min_maf_common_variants=min_maf_common,
                                          variants_per_mac_category=n_mac,
                                          variants_per_maf_category=n_maf)

    # export plink files for GRM construction
    _ = generate_plink_files_for_grm(mt_dict,
                                     sample_qc=sample_qc,
                                     min_af=min_af,
                                     use_drc_pop=use_drc_pop,
                                     overwrite=overwrite,
                                     use_plink=use_plink)
    
    # produce sparse GRMs
    generate_sparse_grm_distributed(pops=POPS, sample_qc=sample_qc,
                                    af_cutoff=min_af,
                                    n_markers=n_markers,
                                    relatedness=relatedness,
                                    no_wait=no_wait,
                                    saige_importers_path=saige_importers_path,
                                    wdl_path=sparse_wdl_path,
                                    use_drc_pop=use_drc_pop,
                                    overwrite=overwrite,
                                    use_plink=use_plink,
                                    use_array_data=use_array_data_for_sparse,
                                    n_common=n_common,
                                    n_maf=n_maf,
                                    n_mac=n_mac)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument()

    args = parser.parse_args()
    main(**args)