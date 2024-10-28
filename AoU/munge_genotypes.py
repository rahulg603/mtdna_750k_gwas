#!/usr/bin/env python3
# This file includes functions to process AoU genomes.
# All pathing to new files is provided by SaigeImporters.py.
# All pathing to AoU files is provided by paths.py.

# Need to determine sample QC and variant QC.
import argparse
import json
import hail as hl
import pandas as pd

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


def get_n_samples_per_pop_vec(analysis_type, sample_qc, use_array_for_variant, use_drc_ancestry_data=False):
    vec_path = get_n_samples_per_pop_path(GENO_PATH, analysis_type=analysis_type, sample_qc=sample_qc, 
                                          use_array_for_variant=use_array_for_variant,
                                          use_drc_ancestry_data=use_drc_ancestry_data)
    if not hl.hadoop_exists(vec_path):
        mt = get_filtered_genotype_mt(analysis_type=analysis_type,
                                      pop='all',
                                      filter_samples=sample_qc,
                                      filter_variants=False,
                                      use_array_for_variant=use_array_for_variant,
                                      use_drc_ancestry_data=use_drc_ancestry_data)
        ht = mt.cols()
        ht_ct = ht.group_by(ht.pop).aggregate(N = hl.agg.count())
        df_ct = ht_ct.to_pandas()
        df_ct.to_csv(vec_path, sep='\t', index=False)
    else:
        df_ct = pd.read_csv(vec_path, sep='\t')
    
    return {x['pop']: x.N for _, x in df_ct.iterrows()}


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
    
    if (analysis_type == 'gene') or (analysis_type == 'variant' and not use_array_for_variant):
        mt = hl.read_matrix_table(mt_path)
        mt = mt.filter_entries(hl.is_missing(mt.FT) | (mt.FT == 'PASS'))
        mt = mt.drop('variant_qc')
    else:
        mt = hl.read_matrix_table(mt_path, _n_partitions=15000)
    
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


def filter_variants_for_null(pop, analysis_type, use_array_for_variant, sample_qc, 
                             use_drc_ancestry_data=False, overwrite=False,
                             n_common_variants_to_keep=50000, # 100000 for per pop
                             min_call_rate=CALLRATE_CUTOFF, min_maf_common_variants=0.01, 
                             variants_per_mac_category=2000, variants_per_maf_category=10000):
    
    ht_sites_path = get_sites_for_null_path(GENO_PATH, extension='ht',
                                            pop=pop, analysis_type=analysis_type, sample_qc=sample_qc,
                                            use_drc_ancestry_data=use_drc_ancestry_data,
                                            use_array_for_variant=use_array_for_variant,
                                            ld_pruned=False,
                                            n_common=n_common_variants_to_keep, 
                                            n_maf=variants_per_maf_category,
                                            n_mac=variants_per_mac_category)
    
    if overwrite or not hl.hadoop_exists(os.path.join(ht_sites_path, '_SUCCESS')):
        print(f'Number of common variants to sample: {n_common_variants_to_keep}')
        n_samples = get_n_samples_per_pop_vec(analysis_type, sample_qc, use_array_for_variant=use_array_for_variant,
                                              use_drc_ancestry_data=use_drc_ancestry_data)
        ht = get_call_stats_ht(pop=pop, sample_qc=sample_qc, analysis_type=analysis_type,
                               use_drc_ancestry_data=use_drc_ancestry_data, 
                               use_array_for_variant=use_array_for_variant, overwrite=overwrite)
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


def filter_mt_for_null(pop, analysis_type, use_array_for_variant, sample_qc, 
                       use_drc_ancestry_data=False, overwrite=False,
                       n_common_variants_to_keep=50000, # 100000 for per pop
                       min_call_rate=CALLRATE_CUTOFF, min_maf_common_variants=0.01, 
                       variants_per_mac_category=2000, variants_per_maf_category=10000):
    
    mt_sites_path = get_sites_for_null_path(GENO_PATH, extension='mt',
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
        filtered_ht = filter_variants_for_null(pop=pop, analysis_type=analysis_type,
                                               sample_qc=sample_qc,
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
        mt = mt.filter_rows(~hl.parse_locus_interval(INVERSION_LOCUS, reference_genome="GRCh38").contains(mt.locus) & \
                            ~hl.parse_locus_interval(HLA_LOCUS, reference_genome="GRCh38").contains(mt.locus))

        mt = mt.naive_coalesce(1000).checkpoint(mt_sites_path)
    else:
        mt = hl.read_matrix_table(mt_sites_path)

    print(mt.count())
    return mt


def get_filtered_array_mt_for_pruning(pop, sample_qc, min_af=0.01,
                                      min_call_rate=CALLRATE_CUTOFF,
                                      use_drc_ancestry_data=False, overwrite=False):
    n_samples = get_n_samples_per_pop_vec(analysis_type='variant', sample_qc=sample_qc, 
                                            use_array_for_variant=True,
                                            use_drc_ancestry_data=use_drc_ancestry_data)
    ht = get_call_stats_ht(pop=pop, sample_qc=sample_qc, analysis_type='variant',
                            use_drc_ancestry_data=use_drc_ancestry_data, 
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
                                        use_array_for_variant=True, use_drc_ancestry_data=use_drc_ancestry_data)
    geno_mt = geno_mt.semi_join_rows(ht)
    geno_mt = geno_mt.unfilter_entries()
    return geno_mt


def plink_ld_pruned_mt(sample_qc, saige_importers_path, wdl_path, min_af=0.01,
                       min_call_rate=CALLRATE_CUTOFF, 
                       use_drc_ancestry_data=False, 
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
                                                          use_drc_ancestry_data=use_drc_ancestry_data,
                                                          af_cutoff=min_af, extension='ht', 
                                                          window='1e7',
                                                          use_plink=True)

        if overwrite or not hl.hadoop_exists(os.path.join(ld_pruned_ht_path, '_SUCCESS')):

            geno_mt = get_filtered_array_mt_for_pruning(pop=pop, sample_qc=sample_qc, min_af=min_af,
                                                        min_call_rate=min_call_rate, 
                                                        use_drc_ancestry_data=use_drc_ancestry_data, 
                                                        overwrite=overwrite)
            
            for chr in AUTOSOMES:
                plink_path = get_plink_inputs_ld_prune(GENO_PATH, pop=pop, chr=chr, sample_qc=sample_qc,
                                                       use_drc_ancestry_data=use_drc_ancestry_data,
                                                       af_cutoff=min_af, pruned=None, extension='bed')
                plink_root = os.path.splitext(plink_path)[0]

                plink_out = get_plink_inputs_ld_prune(GENO_PATH, pop=pop, chr=chr, sample_qc=sample_qc,
                                                     use_drc_ancestry_data=use_drc_ancestry_data,
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
                                                            use_drc_ancestry_data=use_drc_ancestry_data,
                                                            af_cutoff=min_af, extension='ht', 
                                                            window='1e7',
                                                            use_plink=True)

        if overwrite or not hl.hadoop_exists(os.path.join(ld_pruned_ht_path, '_SUCCESS')):

            ht_list = []

            for chr in AUTOSOMES:
                plink_out = get_plink_inputs_ld_prune(GENO_PATH, pop=pop, chr=chr, sample_qc=sample_qc,
                                                    use_drc_ancestry_data=use_drc_ancestry_data,
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
                                      use_array_for_variant=True, use_drc_ancestry_data=use_drc_ancestry_data)
        mt = mt.semi_join_rows(ht_prune)
        mt_dict.update({pop: mt})
    
    return mt_dict


def hail_ld_pruned_mt(pop, sample_qc, min_af=0.01,
                      min_call_rate=CALLRATE_CUTOFF, 
                      use_drc_ancestry_data=False, 
                      overwrite=False):
    ld_pruned_ht_path = get_ld_pruned_array_data_path(GENO_PATH, pop=pop, sample_qc=sample_qc,
                                                      use_drc_ancestry_data=use_drc_ancestry_data,
                                                      af_cutoff=min_af, window='1e7',
                                                      extension='ht', use_plink=False)
    
    if overwrite or not hl.hadoop_exists(os.path.join(ld_pruned_ht_path, '_SUCCESS')):
        geno_mt = get_filtered_array_mt_for_pruning(pop=pop, sample_qc=sample_qc, min_af=min_af,
                                                    min_call_rate=min_call_rate, 
                                                    use_drc_ancestry_data=use_drc_ancestry_data, 
                                                    overwrite=overwrite)
        mt_staging_path = f'{TEMP_PATH}/tmp_mt_genotypes_for_pruning_{pop}.mt'
        geno_mt.write(mt_staging_path, overwrite=True)
        geno_mt = hl.read_matrix_table(mt_staging_path)

        ht_prune = hl.ld_prune(geno_mt.GT,
                               r2=0.1,
                               bp_window_size=int(1e7),
                               memory_per_core=2048)
        ht_prune = ht_prune.checkpoint(ld_pruned_ht_path, overwrite=True)
    else:
        ht_prune = hl.read_table(ld_pruned_ht_path)
    
    print(f'After pruning using Hail for pop {pop}, there are {str(ht_prune.count())} variants.')
    mt = get_filtered_genotype_mt(analysis_type='variant', pop=pop, filter_samples=sample_qc, filter_variants=True,
                                  use_array_for_variant=True, use_drc_ancestry_data=use_drc_ancestry_data)
    mt = mt.semi_join_rows(ht_prune)
    return mt


def generate_plink_files_for_grm(mt_dict, sample_qc, min_af=0.01,
                                 use_drc_ancestry_data=False, 
                                 overwrite=False, use_plink=False):
    for pop, mt in mt_dict.items():
        plink_path = get_ld_pruned_array_data_path(GENO_PATH, pop=pop, sample_qc=sample_qc,
                                                   use_drc_ancestry_data=use_drc_ancestry_data,
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
                                    use_drc_ancestry_data=False,
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
                                     use_drc_ancestry_data=use_drc_ancestry_data,
                                     use_plink=use_plink)
        if overwrite or not hl.hadoop_exists(mtx):
            pops_to_queue.append(pop)
    
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
                'saige_sparse_grm.use_drc_ancestry_data': use_drc_ancestry_data,
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


def main():
    sample_qc = True
    min_af = 0.01
    use_drc_ancestry_data=True
    overwrite=False
    no_wait=False

    n_markers = 2000
    relatedness = 0.125

    saige_importers_path = os.path.join(BUCKET,'scripts/SaigeImporters.py')
    plink_wdl_path = ''
    sparse_wdl_path = ''
    use_plink = True
    
    if use_plink:
        mt_dict = plink_ld_pruned_mt(sample_qc=sample_qc,
                                     saige_importers_path=saige_importers_path,
                                     wdl_path=plink_wdl_path,
                                     min_af=min_af,
                                     use_drc_ancestry_data=use_drc_ancestry_data,
                                     overwrite=overwrite)

    else:
        mt_dict = {}
        for pop in POPS:
            # export GRM sample list
            mt = hail_ld_pruned_mt(pop=pop,
                                   sample_qc=sample_qc,
                                   min_af=min_af,
                                   use_drc_ancestry_data=use_drc_ancestry_data,
                                   overwrite=overwrite)
            mt_dict.update({pop: mt})

    for pop in POPS:
        with hl.hadoop_open(get_aou_samples_file_path(GENO_PATH, pop, sample_qc=sample_qc, use_plink=use_plink,
                                                      use_drc_ancestry_data=use_drc_ancestry_data), 'w') as f:
            f.write('\n'.join(mt.s.collect()) + '\n')

        # produce downsampled plink files
        # TODO finish this

    # export plink files for GRM construction
    _ = generate_plink_files_for_grm(mt_dict,
                                     sample_qc=sample_qc,
                                     min_af=min_af,
                                     use_drc_ancestry_data=use_drc_ancestry_data,
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
                                    use_drc_ancestry_data=use_drc_ancestry_data,
                                    overwrite=overwrite,
                                    use_plink=use_plink)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument()

    args = parser.parse_args()
    main(**args)