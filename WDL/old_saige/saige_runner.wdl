version 1.0

import "https://personal.broadinstitute.org/rahul/saige/saige_sparse_grm.wdl" as saige_tools
import "https://personal.broadinstitute.org/rahul/saige/saige_tests.wdl" as saige_tests

workflow saige_multi {

    input {
        Array[String] pheno
        String suffix
        String pop

        Array[String] export_pheno
        Array[Array[String]] null_model
        Array[Array[Array[String]]] tests
        Array[Array[String]] hail_merge

        File sample_ids

        Float test_min_mac = 0.5
        Float test_min_maf = 0

        File? group_file
        File? groups
        Float? max_maf_for_group

        File bedfile_vr_markers
        File bimfile_vr_markers
        File famfile_vr_markers
        File? sparse_grm
        File? sparse_grm_ids

        String bgen_prefix

        Float rel_cutoff

        String covariates
        File? additional_covariates

        # path to gs folders
        String gs_bucket
        String gs_genotype_path
        String gs_phenotype_path
        String gs_covariate_path
        String gs_output_path
        String? google_project_req_pays

        # options
        Boolean rvas_mode
        Boolean use_drc_ancestry_data
        Boolean always_use_sparse_grm
        Boolean force_inverse_normalize
        Boolean disable_loco

        # runtime
        Int n_cpu_null
        Int n_cpu_test
        Int n_cpu_merge = 16

        # helper functions
        File SaigeImporters

        # docker images
        String HailDocker
        String SaigeDocker
    }

    String analysis_type = if rvas_mode then 'gene' else 'variant'

    scatter (per_pheno_data in zip(zip(zip(zip(hail_merge, tests), null_model), export_pheno), pheno)) {

        if (per_pheno_data.left.right == "") {

            call export_phenotype_files {
                # this function will read in a single phenotype flat file, munge them into a correct format, and output the phenotypes to process
                input:
                    phenotype_id = per_pheno_data.right,
                    pop = pop,
                    suffix = suffix,
                    additional_covariates = additional_covariates,
                    use_drc_ancestry_data = use_drc_ancestry_data,
                    
                    gs_bucket = gs_bucket,
                    gs_phenotype_path = gs_phenotype_path,
                    gs_covariate_path = gs_covariate_path,

                    SaigeImporters = SaigeImporters,
                    HailDocker = HailDocker
            }

        }

        File pheno_file_comb = select_first([export_phenotype_files.pheno_file, per_pheno_data.left.right])

        if (per_pheno_data.left.left.right[0] == '') {

            call null {
                input:
                    # runs the null model
                    phenotype_id = per_pheno_data.right,
                    phenotype_file = pheno_file_comb,

                    pop = pop,
                    suffix = suffix,

                    covariates = covariates,

                    bedfile_vr_markers = bedfile_vr_markers,
                    bimfile_vr_markers = bimfile_vr_markers,
                    famfile_vr_markers = famfile_vr_markers,
                    sparse_grm = sparse_grm,
                    sparse_grm_ids = sparse_grm_ids,

                    gs_bucket = gs_bucket,
                    gs_output_path = gs_output_path,

                    rel_cutoff = rel_cutoff,
                    n_cpu_null = n_cpu_null,

                    analysis_type = analysis_type,
                    force_inverse_normalize = force_inverse_normalize,
                    disable_loco = disable_loco,

                    SaigeImporters = SaigeImporters,
                    SaigeDocker = SaigeDocker
            }

            call saige_tools.upload as u1 {
                input:
                    paths = [null.null_rda_path, null.null_var_path, null.log_path],
                    files = [null.null_rda, null.null_var_ratio, null.log],
                    HailDocker = HailDocker
            }

        }
        File null_rda_comb = select_first([null.null_rda, per_pheno_data.left.left.right[0]])
        File null_var_ratio_comb = select_first([null.null_var_ratio, per_pheno_data.left.left.right[1]])
        File null_log_comb = select_first([null.log, per_pheno_data.left.left.right[2]])

        call saige_tests.saige_tests as test_runner {
            input:
                pheno = per_pheno_data.right,
                suffix = suffix,
                pop = pop,

                null_rda = null_rda_comb,
                null_var_ratio = null_var_ratio_comb,
                sample_list = sample_ids,

                sparse_grm = sparse_grm,
                sparse_grm_ids = sparse_grm_ids,

                bgen_prefix = bgen_prefix,

                tests = per_pheno_data.left.left.left.right,

                group_file = group_file,
                groups = groups,
                max_maf_for_group = max_maf_for_group,

                min_mac = test_min_mac,
                min_maf = test_min_maf,

                gs_bucket = gs_bucket, 
                gs_genotype_path = gs_genotype_path, 
                gs_output_path = gs_output_path, 
                google_project_req_pays = google_project_req_pays, 

                rvas_mode = rvas_mode,
                always_use_sparse_grm = always_use_sparse_grm,
                disable_loco = disable_loco,

                n_cpu_test = n_cpu_test,
                SaigeImporters = SaigeImporters,
                HailDocker = HailDocker,
                SaigeDocker = SaigeDocker
        }

        if (per_pheno_data.left.left.left.left[0] == '') {

            call merge {
                input:
                    phenotype_id = per_pheno_data.right,
                    suffix = suffix,
                    pop = pop,
                    analysis_type = analysis_type,

                    null_log = null_log_comb,
                    test_logs = test_runner.test_logs,
                    single_test = test_runner.single_variant,
                    gene_test = test_runner.gene_test,

                    gs_bucket = gs_bucket, 
                    gs_output_path = gs_output_path, 
                    
                    SaigeImporters = SaigeImporters,
                    HailDocker = HailDocker,

                    n_cpu_merge = n_cpu_merge
            }

            call saige_tools.upload as u2 {
                input:
                    paths = [merge.single_variant_flat_path],
                    files = [merge.single_variant_flat_file],
                    HailDocker = HailDocker
            }
        }
        
        File sumstats = select_first([merge.single_variant_flat_file, per_pheno_data.left.left.left.left[1]])
        
        # call manhattan {
        # }

    }

    output {

    }

}


task export_phenotype_files {
    input {
        
        String phenotype_id
        String pop
        String suffix

        String gs_bucket
        String gs_phenotype_path
        String gs_covariate_path
        File? additional_covariates
        Boolean use_drc_ancestry_data

        File SaigeImporters
        String HailDocker
    }

    String addl_cov_file = select_first([additional_covariates, ''])
    String drc = if use_drc_ancestry_data then 'drc' else 'custom'

    command <<<
        set -e

        python3.8 <<CODE
    import importlib
    import os, sys
    import json
    import hail as hl

    this_temp_path = '/cromwell_root/tmp/'
    hl.init(log='log.log', tmp_dir=this_temp_path, default_reference='GRCh38')

    flpath = os.path.dirname('~{SaigeImporters}')
    scriptname = os.path.basename('~{SaigeImporters}')
    sys.path.append(flpath)
    load_module = importlib.import_module(os.path.splitext(scriptname)[0])
    globals().update(vars(load_module))

    gs_prefix = parse_bucket('~{gs_bucket}')
    gs_phenotype_path = os.path.join(gs_prefix, '~{gs_phenotype_path}'.lstrip('/'))
    gs_covariate_path = os.path.join(gs_prefix, '~{gs_covariate_path}'.lstrip('/'))

    pheno_export_dir = get_pheno_export_dir(gs_phenotype_path, '~{suffix}', '~{pop}')
    pheno_dct = pheno_str_to_dict('~{phenotype_id}')
    pheno_export_path = get_pheno_output_path(pheno_export_dir, pheno_dct)

    addl_cov = None if '~{addl_cov_file}' == '' else '~{addl_cov_file}'
    drc_tf = '~{drc}' == 'drc'

    binary_trait = SAIGE_PHENO_TYPES[pheno_dct['trait_type']] != 'quantitative'

    suffix = '~{suffix}'
    mt = get_custom_ukb_pheno_mt(gs_phenotype_path, gs_covariate_path, addl_cov, suffix, "~{pop}", drc=drc_tf)
    mt = mt.filter_cols(hl.all(lambda x: x, [mt[k] == pheno_dct[k] for k in PHENO_KEY_FIELDS if k != 'pheno_sex']))
    pheno_sex_mt = mt.filter_cols(mt.pheno_sex == pheno_dct['pheno_sex'])
    if pheno_sex_mt.count_cols() == 1:
        mt = pheno_sex_mt
    else:
        mt = mt.filter_cols(mt.pheno_sex == 'both_sexes')
    mt = mt.select_entries(value=mt[pheno_dct['pheno_sex']])
    if binary_trait:
        mt = mt.select_entries(value=hl.int(mt.value))
    
    ht = mt.key_cols_by().select_cols().entries()
    ht.export(pheno_export_path)

    with open('export_path.txt', 'w') as f:
        f.write(pheno_export_path)
    
    CODE
    >>>

    runtime {
        docker: HailDocker
        memory: '4 GB'
    }

    output {
        String pheno_file = read_string('export_path.txt')
    }
}


task null {

    input {
        String phenotype_id
        File phenotype_file

        String pop
        String suffix

        String covariates
        #?qCovarColList # NOTE NEED TO FIGURE THIS OUT
        #isCovariateOffset # NOTE NEED TO FIGURE THIS OUT

        File bedfile_vr_markers
        File bimfile_vr_markers
        File famfile_vr_markers
        File? sparse_grm
        File? sparse_grm_ids

        String gs_bucket
        String gs_output_path

        Boolean force_inverse_normalize
        Boolean disable_loco

        String analysis_type

        Float rel_cutoff

        Int n_cpu_null

        File SaigeImporters
        String SaigeDocker
    }

    String tf_defined_spGRM = if defined(sparse_grm) then "defined" else "not"
    String loco = if disable_loco then "noloco" else "loco"
    String invnorm = if force_inverse_normalize then "inv_normal" else "no_normal"
    String output_prefix = phenotype_id + "." + analysis_type + "." + pop + "." + suffix

    command <<<
        set -e
        set -o pipefail

        python3.8 <<CODE
    import importlib
    import os, sys
    import json
    import subprocess
    import pandas as pd

    flpath = os.path.dirname('~{SaigeImporters}')
    scriptname = os.path.basename('~{SaigeImporters}')
    sys.path.append(flpath)
    load_module = importlib.import_module(os.path.splitext(scriptname)[0])
    globals().update(vars(load_module))

    pheno_dct = pheno_str_to_dict('~{phenotype_id}')
    trait_type = SAIGE_PHENO_TYPES[pheno_dct['trait_type']]
    
    gs_prefix = parse_bucket('~{gs_bucket}')
    gs_output_path = os.path.join(gs_prefix, '~{gs_output_path}'.lstrip('/'))
    null_model_dir = get_null_model_path(gs_output_path, '~{suffix}', '~{pop}')
    rda_path, var_ratio_path, log_path = get_null_model_file_paths(null_model_dir, pheno_dct, '~{analysis_type}')

    # we need to recode the bimfile to have numeric chromosome names
    df = pd.read_csv('~{bimfile_vr_markers}', header=None, sep='\t')
    l1 = [x.replace('chr', '') for x in CHROMOSOMES]
    mapper = {y: '23' if x == 'X' else ('24' if x == 'Y' else x) for x, y in zip(l1, CHROMOSOMES)}
    if df[0].isin(mapper.keys()).any():
        df[0] = df[0].map(mapper)
        df.to_csv('updated_bim.bim', header=None, index=None, sep='\t')
        path_to_bim = 'updated_bim.bim'
    else:
        path_to_bim = '~{bimfile_vr_markers}'

    saige_step_1 = ['Rscript', '/usr/local/bin/step1_fitNULLGLMM.R',
                    '--bedFile=~{bedfile_vr_markers}',
                    f'--bimFile={path_to_bim}',
                    '--famFile=~{famfile_vr_markers}',
                    '--phenoFile=~{phenotype_file}',
                    '--outputPrefix=~{output_prefix}',
                    '--outputPrefix_varRatio=~{output_prefix}',
                    '--covarColList=~{covariates}',
                    '--phenoCol=value',
                    f'--sampleIDColinphenoFile={PHENO_SAMPLE_ID}',
                    f'--traitType={trait_type}',
                    '--minCovariateCount=1',
                    '--nThreads=~{n_cpu_null}']

    if "~{tf_defined_spGRM}" == 'defined':
        saige_step_1 = saige_step_1 + ['--relatednessCutoff=~{rel_cutoff}',
                                       '--sparseGRMFile=~{sparse_grm}',
                                       '--sparseGRMSampleIDFile=~{sparse_grm_ids}',
                                       '--useSparseGRMtoFitNULL=TRUE',
                                       '--useSparseGRMforVarRatio=TRUE',
                                       '--LOCO=FALSE']
    else:
        if '~{loco}' == 'loco':
            saige_step_1 = saige_step_1 + ['--LOCO=TRUE']
        else:
            saige_step_1 = saige_step_1 + ['--LOCO=FALSE']
    
    if "~{analysis_type}" == 'gene':
        saige_step_1 = saige_step_1 + ['--isCateVarianceRatio=TRUE',
                                       '--cateVarRatioMinMACVecExclude=0.5,1.5,2.5,3.5,4.5,5.5,10.5,15.5,20.5',
                                       '--cateVarRatioMaxMACVecInclude=1.5,2.5,3.5,4.5,5.5,10.5,15.5,20.5']

    if "~{invnorm}" == "inv_normal":
        saige_step_1 = saige_step_1 + ['--invNormalize=TRUE']
    else:
        saige_step_1 = saige_step_1 + ['--invNormalize=FALSE']

    with open('~{output_prefix}' + ".log", "wb") as f:
        process = subprocess.Popen(saige_step_1, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for c in iter(lambda: process.stdout.read(1), b""):
            sys.stdout.buffer.write(c)
            f.write(c)

    with open('rda_path.txt', 'w') as f:
        f.write(rda_path)

    with open('null_var_ratio_path.txt', 'w') as f:
        f.write(var_ratio_path)

    with open('null_log_path.txt', 'w') as f:
        f.write(log_path)
    
    CODE

    ls -lh

    >>>

    runtime {
        docker: SaigeDocker
        memory: '64 GB'
        cpu: n_cpu_null
    }

    output {
        File null_rda = output_prefix + '.rda'
        File null_var_ratio = output_prefix + '.varianceRatio.txt'
        File log = output_prefix + ".log"
        String null_rda_path = read_string('rda_path.txt')
        String null_var_path = read_string('null_var_ratio_path.txt')
        String log_path = read_string('null_log_path.txt')
    }
}


task merge {

    input {
        String phenotype_id
        String suffix
        String pop
        String analysis_type

        File null_log
        Array[File] test_logs
        Array[File] single_test
        Array[File?] gene_test
        
        String gs_bucket
        String gs_output_path

        File SaigeImporters
        String HailDocker

        Int n_cpu_merge
    }

    Int disk = ceil((size(single_test, 'G') + size(gene_test, 'G')) * 2)
    String output_prefix = phenotype_id + "." + analysis_type + "." + pop + "." + suffix

    command <<<
        set -e

        python3.8 <<CODE
    import hail as hl
    import importlib
    import os, sys
    from datetime import date

    curdate = date.today().strftime("%y%m%d")

    this_temp_path = '/cromwell_root/tmp/'
    hl.init(master=f'local[str(~{n_cpu_merge})]',
            log='load_results.log', tmp_dir=this_temp_path)

    # import relevant objects
    flpath = os.path.dirname('~{SaigeImporters}')
    scriptname = os.path.basename('~{SaigeImporters}')
    sys.path.append(flpath)
    load_module = importlib.import_module(os.path.splitext(scriptname)[0])
    globals().update(vars(load_module))

    gs_prefix = parse_bucket('~{gs_bucket}')
    gs_output_path = os.path.join(gs_prefix, '~{gs_output_path}'.lstrip('/'))

    pheno_dct = pheno_str_to_dict('~{phenotype_id}')
    trait_type = SAIGE_PHENO_TYPES[pheno_dct['trait_type']]
    result_dir = get_result_path(gs_output_path, '~{suffix}', '~{pop}')
    pheno_results_dir = get_pheno_output_path(result_dir, pheno_dct, '')
    results_prefix = get_results_prefix(pheno_results_dir, pheno_dct, chr)
    results_files = get_results_files(results_prefix, '~{analysis_type}')

    variant_ht = get_merged_ht_path(gs_output_path, "~{suffix}", "~{pop}", pheno_dct)

    ht = load_variant_data(output_ht_path=variant_ht,
                           paths='~{sep="," single_test}'.split(','),
                           extension='',
                           trait_type=trait_type,
                           pheno_dict=pheno_dct,
                           null_log='~{null_log}',
                           test_logs='~{sep="," test_logs}'.split(','))


    gene_ht = get_merged_ht_path(gs_output_path, "~{suffix}", "~{pop}", pheno_dct, gene_analysis=True)
    if "~{analysis_type}" == "gene":
        load_gene_data(output_ht_path=gene_ht,
                       paths='~{sep="," single_test}'.split(','),
                       extension='',
                       trait_type=trait_type,
                       pheno_dict=pheno_dct,
                       null_log='~{null_log}',
                       test_logs='~{sep="," test_logs}'.split(','))

    ht.export('~{output_prefix + ".tsv.bgz"}')

    single_flat = get_merged_flat_path(gs_output_path, "~{suffix}", "~{pop}", pheno_dct)

    with open('flat_path.txt', 'w') as f:
        f.write(single_flat)

    CODE

    >>>

    runtime {
        docker: HailDocker
        memory: '16 GB'
        cpu: n_cpu_merge
        disks: 'local-disk ' + disk + ' SSD'
    }

    output {
        File single_variant_flat_file = output_prefix + ".tsv.bgz"
        String single_variant_flat_path = read_string('flat_path.txt')
    }

}