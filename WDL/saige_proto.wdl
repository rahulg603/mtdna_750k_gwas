version 1.0

workflow saige {

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
        String gs_phenotype_path
        String gs_covariate_path
        String gs_output_path

        # helper functions
        File SaigeImporters

        # options
        Boolean rvas_mode
        Boolean include_base_covariates
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

    call process_phenotype_table {
        # this function will read in a single phenotype flat file, munge it into a correct format, and output the phenotypes to process
        input:
            phenotype_flat_file = phenotype_flat_file,
            suffix = suffix,
            trait_type = trait_type,
            modifier = modifier,
            sample_col = sample_col,
            
            additional_covariates = additional_covariates,
            gs_bucket = gs_bucket,
            gs_phenotype_path = gs_phenotype_path,
            gs_covariate_path = gs_covariate_path,
            include_base_covariates = include_base_covariates,
            num_pcs = num_pcs,
            append_pheno = append_pheno,
            overwrite_pheno = overwrite_pheno,

            SaigeImporters = SaigeImporters
    }

    scatter (pop in pops) {

        call get_tasks_to_run as tasks {
            # this function checks to see which of each subtask is done for this population
            # outputs:
            # phe: [pheno1, pheno2, pheno3, pheno4]
            # export_phe: [File, File, File, File]
            # null: [[File, File], [File, File], [File, File], [File, File]]
            # test: [[[File, File],...,[File, File]], [[File, File],...,[File, File]],...]
            # merge: [String, String, String, String]
            input:
                suffix = suffix,
                pop = pop,

                specific_phenos = specific_phenos,
                min_cases = min_cases,
                sex_stratified = sex_stratified,
                rvas_mode = rvas_mode,

                gs_bucket = gs_bucket,
                gs_phenotype_path = gs_phenotype_path,
                gs_output_path = gs_output_path,

                overwrite_pheno_export = overwrite_pheno_export,
                overwrite_null = overwrite_null,
                overwrite_tests = overwrite_tests,
                overwrite_hail_results = overwrite_hail_results,

                SaigeImporters = SaigeImporters,

                wait_for_pheno_mt = process_phenotype_table.task_complete
        }

        scatter (per_pheno_data in zip(zip(zip(zip(tasks.merge, tasks.test), tasks.null), tasks.export_phe), tasks.phe)) {

            if (per_pheno_data.left.right == '') {
                call export_phenotype_files {
                    # this function will read in a single phenotype flat file, munge them into a correct format, and output the phenotypes to process
                    input:
                        phenotype_id = per_pheno_data.right,
                        pop = pop,
                        suffix = suffix,
                        additional_covariates = additional_covariates,
                        
                        gs_bucket = gs_bucket,
                        gs_phenotype_path = gs_phenotype_path,
                        SaigeImporters = SaigeImporters
                }
            }
            String pheno_file = select_first([export_phenotype_files.pheno_file, per_pheno_data.left.right])
        }

    }

    output {
        Array[Array[String]] pheno_files = pheno_file
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
        Int num_pcs
        Boolean append_pheno
        Boolean overwrite_pheno

        File SaigeImporters
    }

    String addl_cov_file = select_first([additional_covariates, ''])
    String include_base_covar = if include_base_covariates then 'base' else 'no_base'
    String append = if append_pheno then 'append' else 'no_append'
    String overwrite = if overwrite_pheno then 'overwrite' else 'no_overwrite'

    command <<<
        set -e

        mkdir tmp
        export _JAVA_OPTIONS="-Djava.io.tmpdir=$(pwd)/tmp/"

        python3.8 <<CODE
    import hail as hl
    import importlib
    import os, sys
    from datetime import date

    curdate = date.today().strftime("%y%m%d")

    this_temp_path = os.path.abspath('./tmp/')
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

    kwargs = {'data_path': '~{phenotype_flat_file}',
              'trait_type': '~{trait_type}',
              'modifier': '~{modifier}',
              'cov_folder': gs_covariate_path,
              'custom': addl_cov,
              'sample_col': '~{sample_col}'}
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
    mt_path = get_custom_ukb_pheno_mt_path(gs_phenotype_path, suffix)
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
        docker: 'us-docker.pkg.dev/mito-wgs/mito-wgs-docker-repo/rgupta-hail-utils:0.2.119'
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
        String gs_output_path

        Boolean overwrite_pheno_export
        Boolean overwrite_null
        Boolean overwrite_tests
        Boolean overwrite_hail_results

        Boolean rvas_mode

        File SaigeImporters

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

        mkdir tmp
        export _JAVA_OPTIONS="-Djava.io.tmpdir=$(pwd)/tmp/"

        python3.8 <<CODE
    import hail as hl
    import importlib
    import os, sys
    import json

    this_temp_path = os.path.abspath('./tmp/')
    hl.init(log='log.log', tmp_dir=this_temp_path)

    flpath = os.path.dirname('~{SaigeImporters}')
    scriptname = os.path.basename('~{SaigeImporters}')
    sys.path.append(flpath)
    load_module = importlib.import_module(os.path.splitext(scriptname)[0])
    globals().update(vars(load_module))

    gs_prefix = parse_bucket('~{gs_bucket}')
    gs_phenotype_path = os.path.join(gs_prefix, '~{gs_phenotype_path}'.lstrip('/'))
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
                    this_pheno_result_holder.append(['', ''])
                else:
                    this_pheno_result_holder.append([results_files[0], ''])
            else:
                res_found = (results_files[0] in results_already_created) and (results_files[1] in results_already_created)
                if overwrite_test_tf or not res_found:
                    this_pheno_result_holder.append(['', ''])
                else:
                    this_pheno_result_holder.append([results_files[0], results_files[1]])
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

    CODE

    >>>
    
    runtime {
        docker: 'us-docker.pkg.dev/mito-wgs/mito-wgs-docker-repo/rgupta-hail-utils:0.2.119'
        memory: '4 GB'
    }

    output {
        Array[String] phe = read_json("pheno.json")
        Array[String] export_phe = read_json("pheno_export.json")
        Array[Array[String]] null = read_json("null_mod.json")
        Array[Array[Array[String]]] test = read_json("run_tests.json")
        Array[String] merge = read_json("merge.json")
    }
}

task export_phenotype_files {
    input {
        
        String phenotype_id
        String pop
        String suffix

        String gs_bucket
        String gs_phenotype_path
        File? additional_covariates

        File SaigeImporters
    }

    command <<<
        set -e

        python3.8 <<CODE
    import importlib
    import os, sys
    import json

    flpath = os.path.dirname('~{SaigeImporters}')
    scriptname = os.path.basename('~{SaigeImporters}')
    sys.path.append(flpath)
    load_module = importlib.import_module(os.path.splitext(scriptname)[0])
    globals().update(vars(load_module))

    gs_prefix = parse_bucket('~{gs_bucket}')
    gs_phenotype_path = os.path.join(gs_prefix, '~{gs_phenotype_path}'.lstrip('/'))

    pheno_export_dir = get_pheno_export_dir(gs_phenotype_path, '~{suffix}', '~{pop}')
    pheno_dct = pheno_str_to_dict('~{phenotype_id}')
    pheno_export_path = get_pheno_output_path(pheno_export_dir, pheno_dct)

    with open('a.txt', 'w') as f:
        f.write(pheno_export_path)
    >>>

    runtime {
        docker: 'us-docker.pkg.dev/mito-wgs/mito-wgs-docker-repo/rgupta-hail-utils:0.2.119'
        memory: '2 GB'
    }

    output {
        String pheno_file = read_string('a.txt')
    }
}
