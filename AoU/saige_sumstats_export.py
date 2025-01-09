import os
import json
import hail as hl
import pandas as pd
from AoU.paths import *
from utils.SaigeImporters import *
from cromwell.classes import CromwellManager


def distributed_export_meta(wdl_path, saige_importers, suffix, encoding, gene_analysis, legacy_exponentiate_p=True, use_drc_pop=True, use_custom_pcs='custom', n_cpu=32):
    """
    This function exports sumstats with meta-analyses as flat files in the pan ancestry format.
    In its current form, summary statistics without meta-analyses will NOT be exported via this method.
    TODO update to export any summary statistic.
    """
    suffix_updated = update_suffix(suffix, use_drc_pop, use_custom_pcs)
    meta_mt = hl.read_matrix_table(get_saige_meta_mt_path(GWAS_PATH, suffix_updated, encoding, gene_analysis=gene_analysis))
    ht = meta_mt.cols()
    phenotype_list = ht.annotate(phenotype_id = hl.str('-').join([ht[x] for x in PHENO_KEY_FIELDS])).phenotype_id.collect()
    sumstat_files = [format_pheno_dir(x) + '.tsv.bgz' for x in phenotype_list]

    path_to_meta_sumstats = get_saige_sumstats_tsv_folder(GWAS_PATH, suffix_updated, encoding, gene_analysis)
    
    if hl.hadoop_exists(path_to_meta_sumstats):
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

    with open(os.path.abspath('./saige_template.json'), 'w') as j:
        json.dump(baseline, j)

    # run manhattan plotting
    print('Starting flat file export...')
    print('This stage will use Cromwell.')
    manager = CromwellManager(run_name='export_saige_sumstats',
                              inputs_file=df,
                              json_template_path=os.path.abspath('./saige_template.json'),
                              wdl_path=wdl_path,
                              batch=None, limit=None, n_parallel_workflows=99, 
                              add_requester_pays_parameter=False,
                              restart=True, batches_precomputed=False, 
                              submission_sleep=0, check_freq=120, quiet=False)
    #manager.run_pipeline(submission_retries=0, cromwell_timeout=60, skip_waiting=True)
    
    return manager