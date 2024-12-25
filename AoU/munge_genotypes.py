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

from cromwell.classes import CromwellManager

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


def get_n_samples_per_pop_vec(analysis_type, sample_qc, use_array_for_variant, use_drc_pop=False):
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
                                    min_call_rate=CALLRATE_CUTOFF, only_autosomes=False, overwrite=False, ac_filter_override=0):
    n_samples = get_n_samples_per_pop_vec(analysis_type, sample_qc, use_array_for_variant=use_array_for_variant,
                                            use_drc_pop=use_drc_pop)
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
                                    use_plink=True):
    """
    Uses cromwell to distribute across pops.
    """
    pops_to_queue = []

    for pop in pops:
        mtx, _ = get_sparse_grm_path(GENO_PATH, pop=pop,
                                     n_markers=n_markers,
                                     relatedness=relatedness,
                                     sample_qc=sample_qc,
                                     af_cutoff=af_cutoff,
                                     use_drc_pop=use_drc_pop,
                                     use_plink=use_plink)
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
                                restart=False, batches_precomputed=False, 
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
                    chrom_length + 1
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


def create_variant_bgen_split_intervals(pop, git_path, wdl_path, callrate_filter, min_ac,
                                        mean_impute_missing=True, use_drc_pop=True, encoding='additive', 
                                        limit=5000, n_cpu=8):
    interval_list = get_variant_intervals(pop=pop, overwrite=False)
    bgen_prefix = get_wildcard_path_intervals_bgen(GENO_PATH, pop=pop, use_drc_pop=use_drc_pop, encoding=encoding)
    
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
    df.to_csv(os.path.abspath(f'./this_{pop}_run.tsv'), index=False, sep='\t')

    subprocess.run(['tar', '-czf', './saige_wdl.tar.gz', git_path])
    repo_tarball = os.path.join(BUCKET,'scripts/saige_wdl.tar.gz')
    _ = subprocess.run(['gsutil','cp','./saige_wdl.tar.gz',repo_tarball])


    baseline = {'split_bgen_intervals.pop': pop,
                'split_bgen_intervals.sample_qc': True,
                'split_bgen_intervals.variant_qc': True,
                'split_bgen_intervals.use_drc_pop': use_drc_pop,
                'split_bgen_intervals.mean_impute_missing': mean_impute_missing,
                'split_bgen_intervals.call_rate_filter': callrate_filter,
                'split_bgen_intervals.min_ac': min_ac,
                'split_bgen_intervals.analysis_type': 'variant',
                'split_bgen_intervals.encoding': encoding,
                'split_bgen_intervals.repo_tarball': repo_tarball,
                'split_bgen_intervals.n_cpu': n_cpu}
    with open(os.path.abspath(f'./saige_template_{pop}.json'), 'w') as j:
        json.dump(baseline, j)

    # run sparse GRM analysis
    print('NOW COMMENCING GENERATION OF BGENs.')
    print('This stage will use Cromwell.')
    manager = CromwellManager(run_name=f'saige_aou_split_bgen_{pop}',
                              inputs_file=df,
                              json_template_path=os.path.abspath(f'./saige_template_{pop}.json'),
                              wdl_path=wdl_path,
                              limit=limit, n_parallel_workflows=500, 
                              add_requester_pays_parameter=False,
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
                                    use_plink=use_plink)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument()

    args = parser.parse_args()
    main(**args)