import os, re
import json
import hail as hl
import pandas as pd
from AoU.paths import *
from utils.SaigeImporters import *
from cromwell.classes import CromwellManager


def distributed_export(wdl_path, saige_importers, suffix, encoding, gene_analysis,
                       cross_biobank_meta=False, legacy_exponentiate_p=True, use_drc_pop=True, use_custom_pcs='custom', 
                       specific_pop=None, n_cpu=32, 
                       overwrite=False, skip_waiting=True, disable_cache=True):
    """
    This function exports sumstats with meta-analyses as flat files in the pan ancestry format.
    In its current form, summary statistics without meta-analyses will NOT be exported via this method.
    TODO update to export any summary statistic.
    """
    suffix_updated = update_suffix(suffix, use_drc_pop, use_custom_pcs)
    
    if cross_biobank_meta:
        meta_mt = hl.read_matrix_table(get_saige_cross_biobank_meta_mt_path(GWAS_PATH, suffix_updated, encoding, gene_analysis=gene_analysis, pop=specific_pop))
    else:
        if specific_pop is not None:
            raise NotImplementedError('ERROR: specific_pop can only be used if cross_biobank_meta is enabled.')
        meta_mt = hl.read_matrix_table(get_saige_meta_mt_path(GWAS_PATH, suffix_updated, encoding, gene_analysis=gene_analysis))
    
    ht = meta_mt.cols()
    phenotype_list = ht.annotate(phenotype_id = hl.str('-').join([ht[x] for x in PHENO_KEY_FIELDS])).phenotype_id.collect()
    sumstat_files = [format_pheno_dir(x) + '.tsv.bgz' for x in phenotype_list]

    if cross_biobank_meta:
        path_to_meta_sumstats = get_saige_cross_biobank_meta_sumstats_tsv_folder(GWAS_PATH, suffix_updated, encoding, gene_analysis, pop=specific_pop)
    else:
        path_to_meta_sumstats = get_saige_sumstats_tsv_folder(GWAS_PATH, suffix_updated, encoding, gene_analysis)
    
    if hl.hadoop_exists(path_to_meta_sumstats) and not overwrite:
        existing_files = [os.path.basename(x['path']) for x in hl.hadoop_ls(path_to_meta_sumstats) if not x['is_dir']]
        phenos_to_run = []
        pheno_output_locations = []
        for pheno, file in zip(phenotype_list, sumstat_files):
            if file not in existing_files:
                phenos_to_run.append(pheno)
                pheno_output_locations.append(os.path.join(path_to_meta_sumstats, file))

    else:
        phenos_to_run = phenotype_list
        pheno_output_locations = [os.path.join(path_to_meta_sumstats, x) for x in sumstat_files]

    this_run = {'phenotype_id': phenos_to_run, 'output_path': pheno_output_locations}
    df = pd.DataFrame(this_run)
    df.to_csv(os.path.abspath('./this_run.tsv'), index=False, sep='\t') # ADD THIS

    baseline = {'export_single_saige_sumstats.encoding': encoding,
                'export_single_saige_sumstats.suffix': suffix,
                'export_single_saige_sumstats.gene_analysis': gene_analysis,
                'export_single_saige_sumstats.use_drc_pop': use_drc_pop,
                'export_single_saige_sumstats.use_custom_pcs': use_custom_pcs,
                'export_single_saige_sumstats.legacy_exponentiate_p': legacy_exponentiate_p,
                'export_single_saige_sumstats.remove_low_quality_sites': True,
                'export_single_saige_sumstats.gs_bucket': BUCKET.rstrip('/'),
                'export_single_saige_sumstats.gs_gwas_path': remove_bucket(GWAS_PATH),
                'export_single_saige_sumstats.SaigeImporters': saige_importers,
                'export_single_saige_sumstats.n_cpu': n_cpu,
                'export_single_saige_sumstats.mem': 40}
    
    if specific_pop is not None:
        baseline.update({'export_single_saige_sumstats.specific_pop': specific_pop})

    with open(os.path.abspath('./saige_template.json'), 'w') as j:
        json.dump(baseline, j)

    # run manhattan plotting
    print('Starting flat file export...')
    print('This stage will use Cromwell.')
    manager = CromwellManager(run_name='export_saige_sumstats',
                              inputs_file=df,
                              json_template_path=os.path.abspath('./saige_template.json'),
                              wdl_path=wdl_path,
                              batch=None, limit=None, n_parallel_workflows=999, 
                              add_requester_pays_parameter=False,
                              restart=True, batches_precomputed=False, 
                              submission_sleep=0, check_freq=120, quiet=False,
                              disable_cache=disable_cache)
    manager.run_pipeline(submission_retries=0, cromwell_timeout=60, skip_waiting=skip_waiting)
    
    return manager


def list_all_saige_sumstats(suffix, encoding, gene_analysis, cross_biobank_meta=False, use_drc_pop=True, use_custom_pcs='custom'):
    suffix_updated = update_suffix(suffix, use_drc_pop, use_custom_pcs)
    if cross_biobank_meta:
        meta_mt = hl.read_matrix_table(get_saige_cross_biobank_meta_mt_path(GWAS_PATH, suffix_updated, encoding, gene_analysis=gene_analysis))
    else:
        meta_mt = hl.read_matrix_table(get_saige_meta_mt_path(GWAS_PATH, suffix_updated, encoding, gene_analysis=gene_analysis))
    ht = meta_mt.cols()
    phenotype_list = ht.annotate(phenotype_id = hl.str('-').join([ht[x] for x in PHENO_KEY_FIELDS])).phenotype_id.collect()
    sumstat_files = [format_pheno_dir(x) + '.tsv.bgz' for x in phenotype_list]

    if cross_biobank_meta:
        path_to_meta_sumstats = get_saige_cross_biobank_meta_sumstats_tsv_folder(GWAS_PATH, suffix_updated, encoding, gene_analysis)
    else:
        path_to_meta_sumstats = get_saige_sumstats_tsv_folder(GWAS_PATH, suffix_updated, encoding, gene_analysis)
    
    phenos_to_plot = []
    files_to_plot = []

    if hl.hadoop_exists(path_to_meta_sumstats):
        existing_files = [os.path.basename(x['path']) for x in hl.hadoop_ls(path_to_meta_sumstats) if not x['is_dir'] and re.search('.tsv.bgz$', x['path'])]
        for pheno, file in zip(phenotype_list, sumstat_files):
            if file in existing_files:
                phenos_to_plot.append(pheno)
                files_to_plot.append(os.path.join(path_to_meta_sumstats, file))

    return phenos_to_plot, files_to_plot


def export_rvas_meta(export_dir, suffix, use_drc_pop, use_custom_pcs, encoding='additive'):
    mt = hl.read_matrix_table(get_saige_cross_biobank_meta_mt_path(GWAS_PATH, update_suffix(suffix, use_drc_pop, use_custom_pcs), encoding=encoding, gene_analysis=True))
    mt = mt.annotate_cols(pheno_id = hl.str('-').join([mt[x] for x in PHENO_KEY_FIELDS]))
    all_phenos = list(mt.aggregate_cols(hl.agg.collect_as_set(mt.pheno_id)))
    for pheno in all_phenos:
        this_file_id = format_pheno_dir(pheno) + '_merged.tsv.bgz'
        mt_this = mt.select_entries(**mt.meta_analysis[0])
        mt_this = mt_this.annotate_cols(meta_analysis_data = mt_this.meta_analysis_data[0])
        mt_this = mt_this.filter_cols(mt_this.pheno_id == pheno)
        ht_this = mt_this.entries()
        ht_this = ht_this.rename({'gene_symbol': 'Region',
                                  'group': 'Group',
                                  'Pvalue_SKATO': 'Pvalue',
                                  'Pvalue_Burden': 'Pvalue_Burden_Stouffer'})
        ht_this = ht_this.rename({'BETA_IV_Burden': 'Beta_Burden',
                                  'Pvalue_IV_Burden': 'Pvalue_Burden',
                                  'SE_IV_Burden': 'SE_Burden'})
        ht_this = ht_this.key_by('Region','Group','max_MAF').drop('pheno_data', 'meta_analysis_data', *PHENO_KEY_FIELDS, 'pheno_id')
        ht_this.export(os.path.join(export_dir, this_file_id))


def export_cross_biobank_sample_size(suffix, use_drc_pop, use_custom_pcs, encoding='additive', gene_analysis=False):
    mt = hl.read_matrix_table(get_saige_cross_biobank_meta_mt_path(GWAS_PATH, update_suffix(suffix, use_drc_pop, use_custom_pcs), encoding=encoding, gene_analysis=gene_analysis))
    ht = mt.cols()
    ht = ht.select(ht.pheno_data, **ht.meta_analysis_data[0])
    ht = ht.annotate(zipped = hl.zip(ht.cohort, ht.pop))
    ht = ht.annotate(pops = hl.str('|').join(ht.zipped.map(lambda x: x[0] + ':' + hl.str(',').join(x[1]))),
                     ukb_n_cases = ht.pheno_data.filter(lambda x: x.cohort == 'ukb').n_cases.first(),
                     aou_n_cases = ht.pheno_data.filter(lambda x: x.cohort == 'aou').n_cases.first())
    ht = ht.drop('zipped', 'cohort', 'pop', 'n_controls', 'pheno_data')

    return ht