version 1.0

import "https://personal.broadinstitute.org/rahul/saige/saige_runner.wdl" as saige_runner

workflow saige_manager {

    input {

        String suffix
        String trait_type
        String sex = 'both_sexes'
        String modifier

        File phenotype_flat_file
        String sample_col = 's'
        File? additional_covariates

        Array[String] pops = ['eur']

        # path to gs folders
        String gs_bucket
        String gs_genotype_path
        String gs_phenotype_path
        String gs_covariate_path
        String gs_output_path
        String? google_project_req_pays

        # helper functions
        File SaigeImporters

        # docker images
        String SaigeDocker = 'us-docker.pkg.dev/mito-wgs/mito-wgs-docker-repo/saige:1.3.6'
        String HailDocker = 'us-docker.pkg.dev/mito-wgs/mito-wgs-docker-repo/rgupta-hail-utils:0.2.119'

        # options
        Boolean rvas_mode
        Boolean always_use_sparse_grm # Will use sparse GRM for both common- and rare-variant analyses.
                                      # Sparse GRM is always used for rare-variant analyses.
        Boolean force_inverse_normalize # force inverse normalization. Intended for continuous traits only
        Boolean disable_loco # disables leave-one-chromosome-out

        Boolean include_base_covariates
        Boolean use_drc_ancestry_data = true
        String sex_stratified = '' # '' = none, 'all' = all, 'only' = only individual sexes
        String? specific_phenos

        Int num_pcs
        Int min_cases = 50

        Boolean append_pheno = false
        Boolean overwrite_pheno = false
        Boolean overwrite_pheno_export = false
        Boolean overwrite_null = false
        Boolean overwrite_tests = false
        Boolean overwrite_hail_results = false

    }

    String suffix_this = suffix + if use_drc_ancestry_data then '_drccovar' else ''

    call process_phenotype_table {
        # this function will read in a single phenotype flat file, munge it into a correct format, and output the phenotypes to process
        input:
            phenotype_flat_file = phenotype_flat_file,
            suffix = suffix_this,
            trait_type = trait_type,
            modifier = modifier,
            sample_col = sample_col,
            
            additional_covariates = additional_covariates,
            gs_bucket = gs_bucket,
            gs_phenotype_path = gs_phenotype_path,
            gs_covariate_path = gs_covariate_path,
            use_drc_ancestry_data = use_drc_ancestry_data,
            include_base_covariates = include_base_covariates,
            num_pcs = num_pcs,
            append_pheno = append_pheno,
            overwrite_pheno = overwrite_pheno,

            SaigeImporters = SaigeImporters,
            HailDocker = HailDocker
    }

    scatter (pop in pops) {

        call get_tasks_to_run as tasks {
            # this function checks to see which of each subtask is done for this population
            # will also output locations of grm files needed for null model construction
            # outputs:
            # phe: [pheno1, pheno2, pheno3, pheno4]
            # export_phe: [File, File, File, File]
            # null: [[File, File], [File, File], [File, File], [File, File]]
            # test: [[[File, File],...,[File, File]], [[File, File],...,[File, File]],...]
            # merge: [String, String, String, String]
            input:
                suffix = suffix_this,
                pop = pop,

                specific_phenos = specific_phenos,
                min_cases = min_cases,
                sex_stratified = sex_stratified,
                rvas_mode = rvas_mode,
                use_drc_ancestry_data = use_drc_ancestry_data,
                sparse_n_markers = sparse_n_markers,

                gs_bucket = gs_bucket,
                gs_phenotype_path = gs_phenotype_path,
                gs_genotype_path = gs_genotype_path,
                gs_output_path = gs_output_path,

                overwrite_pheno_export = overwrite_pheno_export,
                overwrite_null = overwrite_null,
                overwrite_tests = overwrite_tests,
                overwrite_hail_results = overwrite_hail_results,

                SaigeImporters = SaigeImporters,
                HailDocker = HailDocker,

                wait_for_pheno_mt = process_phenotype_table.task_complete
        }

        
        Array[String] pheno = read_json(tasks.phe)
        Array[String] export_pheno = read_json(tasks.export_phe)
        Array[Array[String]] null_model = read_json(tasks.null)
        Array[Array[Array[String]]] tests = read_json(tasks.test)
        Array[String] hail_merge = read_json(tasks.merge)


        if (always_use_sparse_grm) {
            File bed_subset = get_tasks_to_run.bed_subset
            File bim_subset = get_tasks_to_run.bim_subset
            File fam_subset = get_tasks_to_run.fam_subset
            File sparse_grm = get_tasks_to_run.sparse_grm
            File sparse_grm_ids = get_tasks_to_run.sparse_grm_ids
        }
        if (!always_use_sparse_grm) {
            File bed_all = get_tasks_to_run.bed
            File bim_all = get_tasks_to_run.bim
            File fam_all = get_tasks_to_run.fam
        }

        File bedfile_vr_markers = select_first([bed_subset, bed_all])
        File bimfile_vr_markers = select_first([bim_subset, bim_all])
        File famfile_vr_markers = select_first([fam_subset, fam_all])


        call saige_runner.saige_multi as saige {
            
            input:
                pheno = pheno,
                suffix = suffix_this,
                pop = pop,

                covariates = process_phenotype_table.covariate_list,
                additional_covariates = additional_covariates,

                export_pheno = export_pheno,
                null_model = null_model,
                tests = tests,
                hail_merge = hail_merge,

                bedfile_vr_markers = bedfile_vr_markers,
                bimfile_vr_markers = bimfile_vr_markers,
                famfile_vr_markers = famfile_vr_markers,
                sparse_grm = sparse_grm,
                sparse_grm_ids = sparse_grm_ids,
                
                gs_bucket = gs_bucket,
                gs_genotype_path = gs_genotype_path,
                gs_phenotype_path = gs_phenotype_path,
                gs_covariate_path = gs_covariate_path,
                gs_output_path = gs_output_path,
                google_project_req_pays = google_project_req_pays,
                
                rvas_mode = rvas_mode,
                use_drc_ancestry_data = use_drc_ancestry_data,
                force_inverse_normalize = force_inverse_normalize,
                disable_loco = disable_loco,

                SaigeImporters = SaigeImporters,
                HailDocker = HailDocker,
                SaigeDocker = SaigeDocker
 
        }

    }

    output {
        #Array[Array[String]] pheno_files = pheno_file
    }

}


task process_phenotype_table {

    input {
        File phenotype_flat_file
        File? additional_covariates

        String sample_col
        String suffix
        String trait_type
        String modifier

        String gs_bucket
        String gs_phenotype_path
        String gs_covariate_path

        Boolean include_base_covariates
        Boolean use_drc_ancestry_data
        Int num_pcs
        Boolean append_pheno
        Boolean overwrite_pheno

        File SaigeImporters
        String HailDocker
    }

    String addl_cov_file = select_first([additional_covariates, ''])
    String include_base_covar = if include_base_covariates then 'base' else 'no_base'
    String append = if append_pheno then 'append' else 'no_append'
    String overwrite = if overwrite_pheno then 'overwrite' else 'no_overwrite'
    String drc = if use_drc_ancestry_data then 'drc' else 'custom'

    command <<<
        set -e

        python3.8 <<CODE
    import hail as hl
    import importlib
    import os, sys
    from datetime import date

    curdate = date.today().strftime("%y%m%d")

    this_temp_path = '/cromwell_root/tmp/'
    hl.init(log='log.log', tmp_dir=this_temp_path)

    # import relevant objects
    flpath = os.path.dirname('~{SaigeImporters}')
    scriptname = os.path.basename('~{SaigeImporters}')
    sys.path.append(flpath)
    load_module = importlib.import_module(os.path.splitext(scriptname)[0])
    globals().update(vars(load_module))

    gs_prefix = parse_bucket('~{gs_bucket}')
    gs_phenotype_path = os.path.join(gs_prefix, '~{gs_phenotype_path}'.lstrip('/'))
    gs_covariate_path = os.path.join(gs_prefix, '~{gs_covariate_path}'.lstrip('/'))

    addl_cov = None if '~{addl_cov_file}' == '' else '~{addl_cov_file}'
    drc_tf = '~{drc}' == 'drc'

    kwargs = {'data_path': '~{phenotype_flat_file}',
              'trait_type': '~{trait_type}',
              'modifier': '~{modifier}',
              'cov_folder': gs_covariate_path,
              'custom': addl_cov,
              'sample_col': '~{sample_col}',
              'drc': drc_tf}
    mt, cust_covar_list = load_custom_pheno_with_covariates(**kwargs)

    if '~{include_base_covar}' == 'base':
        basic_covars = BASE_NONPC_COVARS
    else:
        basic_covars = []
    covariates = ','.join(basic_covars + cust_covar_list + [f'PC{x}' for x in range(1, ~{num_pcs} + 1)])

    # write string listing the covariates
    with open('this_covar.txt', 'w') as f:
        f.write(str(covariates))

    suffix = '~{suffix}'
    mt_path = get_custom_ukb_pheno_mt_path(gs_phenotype_path, suffix, drc_tf)
    overwrite_tf = '~{overwrite}' == 'overwrite'
    append_tf = '~{append}' == 'append'

    if not hl.hadoop_exists(f'{mt_path}/_SUCCESS') or (overwrite_tf or append_tf):
        mt_this = mt.group_rows_by('pop').aggregate(
            n_cases=hl.agg.count_where(mt.both_sexes == 1.0),
            n_controls=hl.agg.count_where(mt.both_sexes == 0.0),
            n_defined=hl.agg.count_where(hl.is_defined(mt.both_sexes))
        ).entries()
        mt_this.drop(*[x for x in PHENO_COLUMN_FIELDS if x != 'description' and x in mt_this.row]).show(100, width=180)
        
        # save table
        if append_tf and hl.hadoop_exists(f'{mt_path}/_SUCCESS'):
            original_mt = hl.read_matrix_table(mt_path)
            original_mt = original_mt.checkpoint(get_custom_ukb_pheno_mt_path(gs_phenotype_path, f'{suffix}_before_{curdate}'), overwrite=overwrite_tf)
            original_mt.cols().export(get_custom_phenotype_summary_backup_path(gs_phenotype_path, suffix, curdate))
            original_mt.union_cols(mt, row_join_type='outer').write(mt_path, overwrite=overwrite_tf)
        else:
            mt.write(mt_path, overwrite=overwrite_tf)
        
        summarize_data(gs_phenotype_path, suffix, overwrite=overwrite_tf)

    CODE

    >>>
    
    runtime {
        docker: HailDocker
        memory: '12 GB'
        cpu: '4'
    }

    output {
        String covariate_list = read_string("this_covar.txt")
        Boolean task_complete = true
    }
}


task get_tasks_to_run {
    
    input {
        
        String pop
        String suffix
        String? specific_phenos
        
        String sex_stratified
        Int min_cases

        String gs_bucket
        String gs_phenotype_path
        String gs_genotype_path
        String gs_output_path

        Boolean overwrite_pheno_export
        Boolean overwrite_null
        Boolean overwrite_tests
        Boolean overwrite_hail_results

        Boolean rvas_mode
        Int sparse_n_markers

        File SaigeImporters
        String HailDocker

        Boolean wait_for_pheno_mt
    }

    String overwrite_p = if overwrite_pheno_export then 'overwrite' else 'no_overwrite'
    String overwrite_n = if overwrite_null then 'overwrite' else 'no_overwrite'
    String overwrite_t = if overwrite_tests then 'overwrite' else 'no_overwrite'
    String overwrite_h = if overwrite_hail_results then 'overwrite' else 'no_overwrite'
    String specific_phenos_sel = select_first([specific_phenos, ''])
    String analysis_type = if rvas_mode then 'gene' else 'variant'

    command <<<
        set -e

        python3.8 <<CODE
    import hail as hl
    import importlib
    import os, sys
    import json

    this_temp_path = '/cromwell_root/tmp/'
    hl.init(log='log.log', tmp_dir=this_temp_path)

    flpath = os.path.dirname('~{SaigeImporters}')
    scriptname = os.path.basename('~{SaigeImporters}')
    sys.path.append(flpath)
    load_module = importlib.import_module(os.path.splitext(scriptname)[0])
    globals().update(vars(load_module))

    gs_prefix = parse_bucket('~{gs_bucket}')
    gs_phenotype_path = os.path.join(gs_prefix, '~{gs_phenotype_path}'.lstrip('/'))
    gs_genotype_path = os.path.join(gs_prefix, '~{gs_genotype_path}'.lstrip('/'))
    gs_output_path = os.path.join(gs_prefix, '~{gs_output_path}'.lstrip('/'))

    ht = hl.read_table(get_custom_phenotype_summary_path(gs_phenotype_path, '~{suffix}'))
    ht = ht.filter(ht.pop == '~{pop}')

    criteria = True
    criteria &= (ht.n_cases_by_pop >= ~{min_cases})

    ht = ht.filter(criteria).key_by()

    if len('~{sex_stratified}') > 0:
        ht_sex_specific = ht.annotate(pheno_sex='males').union(ht.annotate(pheno_sex='females'))
        if '~{sex_stratified}' == 'all':
            ht = ht.union(ht_sex_specific)
        else:
            ht = ht_sex_specific

    out = set([tuple(x[field] for field in PHENO_KEY_FIELDS) for x in ht.select(*PHENO_KEY_FIELDS).collect()])
    if len('~{specific_phenos_sel}') > 0:
        specific_phenos = specific_phenos.split(',')
        out = [x for x in out if all(map(lambda y: y is not None, x)) and any([re.match(pcd, '-'.join(x)) for pcd in specific_phenos])]

    pheno_key_dict = [dict(zip(PHENO_KEY_FIELDS, x)) for x in out]
    pheno_export_dir = get_pheno_export_dir(gs_phenotype_path, '~{suffix}', '~{pop}')

    if not ('~{overwrite_p}' == 'overwrite') and hl.hadoop_exists(pheno_export_dir):
        phenos_already_exported = {x['path'] for x in hl.hadoop_ls(pheno_export_dir)}
    else:
        phenos_already_exported = {}

    phenotypes = []
    export_phenotype = []
    null_model = []
    run_tests = []
    run_hail_merge = []
    null_model_dir = get_null_model_path(gs_output_path, '~{suffix}', '~{pop}')
    result_dir = get_result_path(gs_output_path, '~{suffix}', '~{pop}')
    for pheno_dct in pheno_key_dict:
        pheno_name = pheno_dict_to_str(pheno_dct)
        phenotypes.append(pheno_name)
        
        # phenotype export
        pheno_export_path = get_pheno_output_path(pheno_export_dir, pheno_dct)
        if ('~{overwrite_p}' == 'overwrite') or (pheno_export_path not in phenos_already_exported):
            export_phenotype.append('')
        else:
            export_phenotype.append(pheno_export_path)

        # null model
        rda, var_ratio = get_null_model_file_paths(null_model_dir, pheno_dct, '~{analysis_type}')
        overwrite_null_tf = ('~{overwrite_n}' == 'overwrite')
        if not overwrite_null_tf and hl.hadoop_exists(null_model_dir):
            null_models_existing = {x['path'] for x in hl.hadoop_ls(null_model_dir)}
        else:
            null_models_existing = {}
        files_found = rda in null_models_existing and var_ratio in null_models_existing
        if overwrite_null_tf or (not files_found):
            null_model.append(['',''])
        else:
            null_model.append([rda, var_ratio])

        # tests
        this_pheno_result_holder = []
        overwrite_test_tf = ('~{overwrite_t}' == 'overwrite')
        pheno_results_dir = get_pheno_output_path(result_dir, pheno_dct, '')
        if not overwrite_test_tf and hl.hadoop_exists(pheno_results_dir):
            results_already_created = {x['path'] for x in hl.hadoop_ls(pheno_results_dir)}
        else:
            results_already_created = {}
        
        for chr in CHROMOSOMES:
            results_prefix = get_results_prefix(pheno_results_dir, pheno_key_dict, chr)
            results_files = get_results_files(results_prefix, '~{analysis_type}')
            if '~{analysis_type}' == 'variant':
                res_found = results_files[0] in results_already_created
                if overwrite_test_tf or not res_found:
                    this_pheno_result_holder.append([chr, '', ''])
                else:
                    this_pheno_result_holder.append([chr, results_files[0], ''])
            else:
                res_found = (results_files[0] in results_already_created) and (results_files[1] in results_already_created)
                if overwrite_test_tf or not res_found:
                    this_pheno_result_holder.append([chr, '', ''])
                else:
                    this_pheno_result_holder.append([chr, results_files[0], results_files[1]])
        run_tests.append(this_pheno_result_holder)
        
        # merged hail table
        merged_ht_path = get_merged_ht_path(gs_output_path, '~{suffix}', '~{pop}', pheno_dct)
        overwrite_hail_tf = ('~{overwrite_h}' == 'overwrite')
        if overwrite_test_tf or overwrite_hail_tf or \
                merged_ht_path not in results_already_created or \
                not hl.hadoop_exists(f'{merged_ht_path}/_SUCCESS'):
            run_hail_merge.append('')
        else:
            run_hail_merge.append(merged_ht_path)

    # phenotype names
    with open('pheno.json', 'w') as f:
        json.dump(phenotypes, f)
    
    # tf phenotype export
    with open('pheno_export.json', 'w') as f:
        json.dump(export_phenotype, f)
    
    # tf null model
    with open('null_mod.json', 'w') as f:
        json.dump(null_model, f)
    
    # tf tests
    with open('run_tests.json', 'w') as f:
        json.dump(run_tests, f)
    
    # tf hail merge
    with open('merge.json', 'w') as f:
        json.dump(run_hail_merge, f)

    #### NOW generate paths for null model construction
    bed, bim, fam = get_plink_for_grm_path(gs_genotype_path, '~{pop}', '~{analysis_type}')
    bed_subset, bim_subset, fam_subset = get_plink_for_grm_path(gs_genotype_path, '~{pop}', '~{analysis_type}', ~{sparse_n_markers})
    mtx, ix = get_sparse_grm_path(gs_genotype_path, '~{pop}', '~{analysis_type}', ~{sparse_n_markers})

    with open('bed.txt', 'w') as f:
        f.write(bed)

    with open('bim.txt', 'w') as f:
        f.write(bim)

    with open('fam.txt', 'w') as f:
        f.write(fam)

    with open('bed_subset.txt', 'w') as f:
        f.write(bed_subset)

    with open('bim_subset.txt', 'w') as f:
        f.write(bim_subset)

    with open('fam_subset.txt', 'w') as f:
        f.write(fam_subset)

    with open('mtx.txt', 'w') as f:
        f.write(mtx)

    with open('idx.txt', 'w') as f:
        f.write(idx)

    CODE

    >>>
    
    runtime {
        docker: HailDocker
        memory: '4 GB'
    }

    output {
        File phe = "pheno.json"
        File export_phe = "pheno_export.json"
        File null = "null_mod.json"
        File test = "run_tests.json"
        File merge = "merge.json"

        String bed = read_string('bed.txt')
        String bim = read_string('bim.txt')
        String fam = read_string('fam.txt')

        String bed_subset = read_string('bed_subset.txt')
        String bim_subset = read_string('bim_subset.txt')
        String fam_subset = read_string('fam_subset.txt')

        String mtx = read_string('mtx.txt')
        String ix = read_string('ix.txt')
    }
}

