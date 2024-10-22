version 1.0

workflow saige_multi {

    input {
        Array[String] pheno
        String suffix
        String pop

        Array[String] export_pheno
        Array[Array[String]] null_model
        Array[Array[Array[String]]] tests
        Array[String] hail_merge

        File bedfile_vr_markers
        File bimfile_vr_markers
        File famfile_vr_markers
        File? sparse_grm
        File? sparse_grm_ids

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

        File pheno_file = select_first([export_phenotype_files.pheno_file, per_pheno_data.left.right])

        if (per_pheno_data.left.left.right[0] == '') {

            call null {
                # runs the null model
                phenotype_id = per_pheno_data.right,
                pheno_file = pheno_file,
                covariates = process_phenotype_table.covariates,

                bedfile_vr_markers = bedfile_vr_markers,
                bimfile_vr_markers = bimfile_vr_markers,
                famfile_vr_markers = famfile_vr_markers,
                sparse_grm = sparse_grm,
                sparse_grm_ids = sparse_grm_ids,

                analysis_type = analysis_type,
                force_inverse_normalize = force_inverse_normalize,
                disable_loco = disable_loco,

                SaigeImporters = SaigeImporters,
                SaigeDocker = SaigeDocker
            }

        }
        File null_rda = select_first([null.null_rda, per_pheno_data.left.left.right])
        File null_var_ratio = select_first([null.null_var_ratio, per_pheno_data.left.left.right])

        # scatter (chr_run in per_pheno_data.left.left.left.right) {
        #     if (chr_run[0] != '') {
        #         call run_tests {

        #         }
        #     }
        #     if rvas_mode {
        #         File gene_file = select_first([run_tests.gene_file, chr_run[1]])
        #     }
        #     File marker_file = select_first([run_tests.marker_file, chr_run[0]])
        # }

        # if (per_pheno_data.left.left.left.left == '') {
        #     call merge {

        #     }
        # }
        # String ht_path = select_first([merge.ht_path, per_pheno_data.left.left.left.left])

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

        String covariates
        ?qCovarColList # NOTE NEED TO FIGURE THIS OUT
        isCovariateOffset # NOTE NEED TO FIGURE THIS OUT

        File bedfile_vr_markers
        File bimfile_vr_markers
        File famfile_vr_markers
        File? sparse_grm
        File? sparse_grm_ids

        Boolean force_inverse_normalize
        Boolean disable_loco

        String analysis_type

        File SaigeImporters
        String SaigeDocker
    }
    String tf_defined_spGRM = if defined(sparse_grm) then "defined" else "not"
    String loco = if disable_loco then "noloco" else "loco"
    String invnorm = if force_inverse_normalize then "inv_normal" else "no_normal"

    command <<<
        set -e
        set -o pipefail

        python3.8 <<CODE
    import importlib
    import os, sys
    import json
    import subprocess

    flpath = os.path.dirname('~{SaigeImporters}')
    scriptname = os.path.basename('~{SaigeImporters}')
    sys.path.append(flpath)
    load_module = importlib.import_module(os.path.splitext(scriptname)[0])
    globals().update(vars(load_module))

    gs_prefix = parse_bucket('~{gs_bucket}')
    gs_phenotype_path = os.path.join(gs_prefix, '~{gs_phenotype_path}'.lstrip('/'))
    gs_covariate_path = os.path.join(gs_prefix, '~{gs_covariate_path}'.lstrip('/'))

    pheno_dct = pheno_str_to_dict('~{phenotype_id}')
    rda, var_ratio = get_null_model_file_paths(null_model_dir, pheno_dct, '~{analysis_type}')
    trait_type = SAIGE_PHENO_TYPES[pheno_dct['trait_type']]

    saige_step_1 = ['Rscript', '/usr/local/bin/step1_fitNULLGLMM.R',
                    '--bedFile', '~{bedfile_vr_markers}',
                    '--bimFile', '~{bimfile_vr_markers}',
                    '--famFile', '~{famfile_vr_markers}',
                    '--phenoFile', '~{phenotype_file}',
                    '--outputPrefix', '~{phenotype_id}',
                    '--outputPrefix_varRatio', '~{phenotype_id}' + '.' + '~{analysis_type}',
                    '--covarColList', '~{covariates}',
                    '--phenoCol', 'value',
                    '--sampleIDColinphenoFile', PHENO_SAMPLE_ID,
                    '--traitType', trait_type,
                    '--minCovariateCount', 1]

    if "~{tf_defined_spGRM}" == 'defined':
        saige_step_1 = saige_step_1 + ['--relatednessCutoff', 0.05,
                                       '--sparseGRMFile', '~{sparse_grm}',
                                       '--sparseGRMSampleIDFile', '~{sparse_grm_ids}',
                                       '--useSparseGRMtoFitNULL', "TRUE",
                                       '--LOCO', "FALSE"]
    else:
        saige_step_1 = saige_step_1 + ['--nThreads', 16]
        if '~{loco}' == 'loco':
            saige_step_1 = saige_step_1 + ['--LOCO', "TRUE"]
        else:
            saige_step_1 = saige_step_1 + ['--LOCO', "FALSE"]
    
    if "~{analysis_type}" == 'gene':
        saige_step_1 = saige_step_1 + ['--isCateVarianceRatio', "TRUE"]

    if "~{invnorm}" == "inv_normal":
        saige_step_1 = saige_step_1 + ['--invNormalize', "TRUE"]
    else:
        saige_step_1 = saige_step_1 + ['--invNormalize', "FALSE"]

    with open('~{phenotype_id}' + ".log", "wb") as f:
        process = subprocess.Popen(saige_step_1, stdout=subprocess.PIPE)
        for c in iter(lambda: process.stdout.read(1), b""):
            sys.stdout.buffer.write(c)
            f.buffer.write(c)

    CODE


    >>>

    runtime {
        docker: SaigeDocker
        memory: '64 GB'
        cpu: '16'
    }

    output {
        String pheno_file = read_string('export_path.txt')
    }
}