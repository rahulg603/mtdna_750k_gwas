version 1.0

import "https://personal.broadinstitute.org/rahul/saige/saige_sparse_grm.wdl" as saige_tools
import "https://personal.broadinstitute.org/rahul/saige/saige_tests.wdl" as saige_tests
import "https://personal.broadinstitute.org/rahul/saige/ManhattanPlotter.wdl" as ManhattanPlotter

workflow saige_manager {

    input {

        String pop = 'eur'
        String suffix
        String trait_type
        String sex = 'both_sexes'
        String modifier
        String encoding = 'additive'

        File phenotype_flat_file
        String sample_col = 's'
        File? additional_covariates
        String covariate_list

        # files and parameters for test
        Float test_min_mac = 0.5
        Float test_min_maf = 0
        Int test_markers_per_chunk = 10000
        File? group_file
        File? groups
        Float? max_maf_for_group

        # QC parameters (for low_confidence filter)
        Float min_call_rate

        # constants for pathing
        Int sparse_n_markers
        Float sparse_min_af
        Float sparse_relatedness_cutoff
        
        Int n_markers_common = 50000
        Int n_markers_rare_maf = 10000
        Int n_markers_rare_mac = 2000

        # path to gs folders
        String gs_bucket
        String gs_genotype_path
        String gs_phenotype_path
        String gs_covariate_path
        String gs_temp_path
        String gs_output_path
        String? google_project_req_pays

        # helper functions
        File SaigeImporters

        # docker images
        String SaigeDocker = 'us-docker.pkg.dev/mito-wgs/mito-wgs-docker-repo/saige:1.3.6'
        String HailDocker = 'us-docker.pkg.dev/mito-wgs/mito-wgs-docker-repo/rgupta-hail-utils:0.2.119'

        # options
        Boolean use_plink
        Boolean rvas_mode
        Boolean always_use_sparse_grm # Will use sparse GRM for both common- and rare-variant analyses.
                                      # Sparse GRM is always used for rare-variant analyses.
        Boolean force_inverse_normalize # force inverse normalization. Intended for continuous traits only
        Boolean disable_loco # disables leave-one-chromosome-out

        Boolean use_drc_pop = true
        String use_custom_pcs = 'custom' # can be custom, axaou, or none
        Boolean sample_qc = true

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

        Int n_cpu_null
        Int n_cpu_test
        Int n_cpu_merge = 32
        Int mem_merge = 96

    }

    String suffix_this = suffix + (if use_drc_pop then '_drcpop' else '') + (if use_custom_pcs == 'custom' then '_custompcs' else (if use_custom_pcs == 'axaou' then '_axaoupcs' else ''))
    String analysis_type = if rvas_mode then 'gene' else 'variant'

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
            encoding = encoding,

            specific_phenos = specific_phenos,
            min_cases = min_cases,
            sex_stratified = sex_stratified,
            rvas_mode = rvas_mode,
            use_drc_pop = use_drc_pop,
            sample_qc = sample_qc,

            sparse_n_markers = sparse_n_markers,
            min_maf = sparse_min_af,
            relatedness = sparse_relatedness_cutoff,
            n_markers_common = n_markers_common,
            n_markers_rare_maf = n_markers_rare_maf,
            n_markers_rare_mac = n_markers_rare_mac,

            use_plink = use_plink,

            gs_bucket = gs_bucket,
            gs_phenotype_path = gs_phenotype_path,
            gs_genotype_path = gs_genotype_path,
            gs_output_path = gs_output_path,

            overwrite_pheno_export = overwrite_pheno_export,
            overwrite_null = overwrite_null,
            overwrite_tests = overwrite_tests,
            overwrite_hail_results = overwrite_hail_results,

            SaigeImporters = SaigeImporters,
            HailDocker = HailDocker
    }

    
    Array[String] pheno = read_json(tasks.phe)
    Array[String] export_pheno = read_json(tasks.export_phe)
    Array[Array[String]] null_model = read_json(tasks.null)
    Array[Array[Array[String]]] tests = read_json(tasks.test)
    Array[Array[String]] hail_merge = read_json(tasks.merge)


    if (always_use_sparse_grm) {
        File sparse_grm = tasks.mtx
        File sparse_grm_ids = tasks.ix
    }

    File bedfile_vr_markers = tasks.bed
    File bimfile_vr_markers = tasks.bim
    File famfile_vr_markers = tasks.fam
    File sample_ids = tasks.sample_ids


    scatter (per_pheno_data in zip(zip(zip(zip(hail_merge, tests), null_model), export_pheno), pheno)) {

        if (per_pheno_data.left.right == "") {

            call export_phenotype_files {
                # this function will read in a single phenotype flat file, munge them into a correct format, and output the phenotypes to process
                input:
                    phenotype_id = per_pheno_data.right,
                    pop = pop,
                    suffix = suffix_this,
                    additional_covariates = additional_covariates,
                    use_drc_pop = use_drc_pop,
                    use_custom_pcs = use_custom_pcs,
                    
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
                    suffix = suffix_this,

                    covariates = covariate_list,

                    bedfile_vr_markers = bedfile_vr_markers,
                    bimfile_vr_markers = bimfile_vr_markers,
                    famfile_vr_markers = famfile_vr_markers,
                    sparse_grm = sparse_grm,
                    sparse_grm_ids = sparse_grm_ids,

                    gs_bucket = gs_bucket,
                    gs_output_path = gs_output_path,

                    rel_cutoff = sparse_relatedness_cutoff,
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
                suffix = suffix_this,
                pop = pop,
                encoding = encoding,

                null_rda = null_rda_comb,
                null_var_ratio = null_var_ratio_comb,
                sample_list = sample_ids,

                sparse_grm = sparse_grm,
                sparse_grm_ids = sparse_grm_ids,

                tests = per_pheno_data.left.left.left.right,

                group_file = group_file,
                groups = groups,
                max_maf_for_group = max_maf_for_group,

                min_mac = test_min_mac,
                min_maf = test_min_maf,
                markers_per_chunk = test_markers_per_chunk,

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
                    suffix = suffix_this,
                    pop = pop,
                    analysis_type = analysis_type,
                    encoding = encoding,

                    null_log = null_log_comb,
                    test_logs = test_runner.test_logs,
                    single_test = test_runner.single_variant,
                    gene_test = test_runner.gene_test,

                    gs_bucket = gs_bucket, 
                    gs_genotype_path = gs_genotype_path, 
                    gs_output_path = gs_output_path, 
                    gs_temp_path = gs_temp_path,
                    
                    min_call_rate = min_call_rate,
                    use_drc_pop = use_drc_pop,
                    sample_qc = sample_qc,

                    SaigeImporters = SaigeImporters,
                    HailDocker = HailDocker,

                    n_cpu_merge = n_cpu_merge,
                    mem = mem_merge
            }

            call saige_tools.upload as u2 {
                input:
                    paths = [merge.single_variant_flat_path],
                    files = [merge.single_variant_flat_file],
                    HailDocker = HailDocker
            }
        }
        
        File sumstats = select_first([merge.single_variant_flat_file, per_pheno_data.left.left.left.left[1]])

        call ManhattanPlotter.ManhattanPlotter as manhattan {
            input:
                sumstats = sumstats,
                pop = pop,
                pheno = per_pheno_data.right,
                suffix = suffix_this,

                p_col = 'Pvalue',
                af_col = 'AF_Allele2',
                conf_col = 'low_confidence',

                keep_x = true,
                exponentiate_p = false,
                var_as_rsid = true,

                wid = 1300, 
                hei = 640, 
                cex = 1.1, 
                point_size = 18,

                mem = 80
        }

        String results_prefix = per_pheno_data.left.left.left.left[2]

        call saige_tools.upload as u3 {
            input:
                paths = [results_prefix + '.manhattan.png', results_prefix + '.qq.png', results_prefix + '.suggestive.tsv', results_prefix + '.suggestive_genes.tsv'],
                files = [manhattan.manhattan, manhattan.qq, manhattan.sugg, manhattan.sugg_gene],
                HailDocker = HailDocker
        }

    }

    output {
        Array[File] manhattan_plot = manhattan.manhattan
        Array[File] qq_plot = manhattan.qq
        Array[File] sugg_table = manhattan.sugg
        Array[File] sugg_gene_table = manhattan.sugg_gene
    }

}


task get_tasks_to_run {
    
    input {
        
        String pop
        String suffix
        String? specific_phenos
        String encoding
        
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

        Int sparse_n_markers
        Float min_maf
        Float relatedness

        Int n_markers_common
        Int n_markers_rare_maf
        Int n_markers_rare_mac

        Boolean rvas_mode
        Boolean use_drc_pop
        Boolean sample_qc
        Boolean use_plink

        File SaigeImporters
        String HailDocker
    }

    String overwrite_p = if overwrite_pheno_export then 'overwrite' else 'no_overwrite'
    String overwrite_n = if overwrite_null then 'overwrite' else 'no_overwrite'
    String overwrite_t = if overwrite_tests then 'overwrite' else 'no_overwrite'
    String overwrite_h = if overwrite_hail_results then 'overwrite' else 'no_overwrite'
    String specific_phenos_sel = select_first([specific_phenos, ''])
    String analysis_type = if rvas_mode then 'gene' else 'variant'
    String drc = if use_drc_pop then 'drc' else 'custom'
    String qc = if sample_qc then 'qc' else 'no_qc'
    String plink = if use_plink then 'plink' else 'no_plink'

    command <<<
        set -e

        python3.8 <<CODE
    import hail as hl
    import pandas as pd
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
    drc_tf = '~{drc}' == 'drc'
    sample_qc_tf = '~{qc}' == 'qc'
    plink_tf = '~{plink}' == 'plink'

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
    result_dir = get_result_path(gs_output_path, '~{suffix}', '~{pop}', "~{encoding}")
    for pheno_dct in pheno_key_dict:
        running_override = False
        pheno_name = pheno_dict_to_str(pheno_dct)
        phenotypes.append(pheno_name)
        
        # phenotype export
        pheno_export_path = get_pheno_output_path(pheno_export_dir, pheno_dct)
        if ('~{overwrite_p}' == 'overwrite') or (pheno_export_path not in phenos_already_exported) or running_override:
            export_phenotype.append('')
            running_override = True
        else:
            export_phenotype.append(pheno_export_path)

        # null model
        rda, var_ratio, null_log = get_null_model_file_paths(null_model_dir, pheno_dct, '~{analysis_type}')
        overwrite_null_tf = ('~{overwrite_n}' == 'overwrite')
        if not overwrite_null_tf and hl.hadoop_exists(null_model_dir):
            null_models_existing = {x['path'] for x in hl.hadoop_ls(null_model_dir)}
        else:
            null_models_existing = {}
        files_found = rda in null_models_existing and var_ratio in null_models_existing
        if overwrite_null_tf or (not files_found) or running_override:
            null_model.append(['','', ''])
            running_override = True
        else:
            null_model.append([rda, var_ratio, null_log])

        # tests
        this_pheno_result_holder = []
        overwrite_test_tf = ('~{overwrite_t}' == 'overwrite')
        pheno_results_dir = get_pheno_output_path(result_dir, pheno_dct, '')
        if not overwrite_test_tf and hl.hadoop_exists(pheno_results_dir):
            results_already_created = {x['path'] for x in hl.hadoop_ls(pheno_results_dir)}
        else:
            results_already_created = {}

        interval_list = read_variant_intervals(gs_genotype_path, pop="~{pop}", analysis_type="~{analysis_type}")
        bgen_prefix = get_wildcard_path_intervals_bgen(gs_genotype_path, 
                                                       pop="~{pop}", 
                                                       use_drc_pop=drc_tf, 
                                                       encoding="~{encoding}")

        any_test_run = False
        for _, row in interval_list.iterrows():
            this_segment = stringify_interval(row['chrom'], row['start'], row['end'])
            this_bgen_prefix = bgen_prefix.replace('@', row['chrom']).replace('#', str(row['start'])).replace('?', str(row['end']))
            results_prefix = get_results_prefix(pheno_results_dir, pheno_dct, this_segment)
            results_files = get_results_files(results_prefix, '~{analysis_type}')

            if '~{analysis_type}' == 'variant':
                res_found = results_files[0] in results_already_created
                if overwrite_test_tf or not res_found or running_override:
                    this_pheno_result_holder.append([this_segment, '', '', '', row['chrom'], this_bgen_prefix])
                    any_test_run = True
                else:
                    this_pheno_result_holder.append([this_segment, results_files[0], '', results_files[2], row['chrom'], this_bgen_prefix])
            else:
                res_found = (results_files[0] in results_already_created) and (results_files[1] in results_already_created)
                if overwrite_test_tf or not res_found or running_override:
                    this_pheno_result_holder.append([this_segment, '', '', '', row['chrom'], this_bgen_prefix])
                    any_test_run = True
                else:
                    this_pheno_result_holder.append([this_segment, results_files[0], results_files[1], results_files[2], row['chrom'], this_bgen_prefix])

        running_override = True if any_test_run else running_override
        run_tests.append(this_pheno_result_holder)
        
        # merged hail table
        merged_ht_path = get_merged_ht_path(gs_output_path, '~{suffix}', '~{pop}', pheno_dct, '~{encoding}')
        merged_flat_path = get_merged_flat_path(gs_output_path, '~{suffix}', '~{pop}', pheno_dct, '~{encoding}')
        results_path = f'{get_pheno_output_path(result_dir, pheno_dct, "")}/{get_pheno_output_suffix(pheno_dct)}.' + '~{pop}' + '.' + '~{suffix}'

        overwrite_hail_tf = ('~{overwrite_h}' == 'overwrite')
        if overwrite_test_tf or overwrite_hail_tf or \
                merged_ht_path not in results_already_created or \
                not hl.hadoop_exists(f'{merged_ht_path}/_SUCCESS') or running_override:
            run_hail_merge.append(['', '', results_path])
        else:
            run_hail_merge.append([merged_ht_path, merged_flat_path, results_path])

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
    bed, bim, fam = get_plink_for_null_path(geno_folder=gs_genotype_path, 
                                            pop='~{pop}', 
                                            sample_qc=sample_qc_tf, 
                                            analysis_type='both',
                                            ld_pruned=True,
                                            n_common=~{n_markers_common},
                                            n_maf=~{n_markers_rare_maf},
                                            n_mac=~{n_markers_rare_mac},
                                            use_drc_pop=drc_tf, 
                                            use_array_for_variant=False)
    mtx, ix = get_sparse_grm_path(geno_folder=gs_genotype_path, 
                                  pop='~{pop}', 
                                  n_markers=~{sparse_n_markers}, 
                                  relatedness=~{relatedness}, 
                                  sample_qc=sample_qc_tf, 
                                  use_plink=plink_tf,
                                  use_drc_pop=drc_tf, 
                                  af_cutoff=~{min_maf})

    with open('bed.txt', 'w') as f:
        f.write(bed)

    with open('bim.txt', 'w') as f:
        f.write(bim)

    with open('fam.txt', 'w') as f:
        f.write(fam)

    with open('mtx.txt', 'w') as f:
        f.write(mtx)

    with open('ix.txt', 'w') as f:
        f.write(ix)


    #### Get sample IDs file
    sample_id = get_aou_samples_file_path(geno_folder=gs_genotype_path,
                                          pop='~{pop}',
                                          sample_qc=sample_qc_tf,
                                          use_plink=plink_tf,
                                          use_drc_pop=drc_tf)
    
    with open('samp.txt', 'w') as f:
        f.write(sample_id)

    CODE

    >>>
    
    runtime {
        docker: HailDocker
        memory: '4 GB'
        preemptible: 5
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

        String mtx = read_string('mtx.txt')
        String ix = read_string('ix.txt')

        String sample_ids = read_string('samp.txt')
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
        Boolean use_drc_pop
        String use_custom_pcs

        File SaigeImporters
        String HailDocker
    }

    String addl_cov_file = select_first([additional_covariates, ''])
    String use_drc_pop_str = if use_drc_pop then 'drc' else 'aou'

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
    drc_tf = '~{use_drc_pop_str}' == 'drc'

    binary_trait = SAIGE_PHENO_TYPES[pheno_dct['trait_type']] != 'quantitative'

    suffix = '~{suffix}'
    mt = get_custom_ukb_pheno_mt(gs_phenotype_path, gs_covariate_path, addl_cov, suffix, "~{pop}", use_drc_pop=drc_tf, use_custom_pcs="~{use_custom_pcs}")
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
        preemptible: 5
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
        preemptible: 5
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
        String encoding

        File null_log
        Array[File] test_logs
        Array[File] single_test
        Array[File?] gene_test
        
        String gs_bucket
        String gs_genotype_path
        String gs_output_path
        String gs_temp_path

        Float min_call_rate
        Boolean use_drc_pop
        Boolean sample_qc

        File SaigeImporters
        String HailDocker

        Int n_cpu_merge
        Int mem
    }

    Int disk = ceil((size(single_test, 'G') + size(gene_test, 'G')) * 6)
    String output_prefix = phenotype_id + "." + analysis_type + "." + pop + "." + suffix
    String drc = if use_drc_pop then 'drc' else 'custom'
    String qc = if sample_qc then 'qc' else 'no_qc'

    command <<<
        set -e

        python3.8 <<CODE
    import hail as hl
    import importlib
    import os, sys
    from datetime import date

    curdate = date.today().strftime("%y%m%d")

    this_temp_path = '/cromwell_root/tmp/'
    hl.init(master=f'local[~{n_cpu_merge}]',
            log='load_results.log', tmp_dir=this_temp_path,
            spark_conf={'spark.driver.memory': '~{mem}g'})

    # import relevant objects
    flpath = os.path.dirname('~{SaigeImporters}')
    scriptname = os.path.basename('~{SaigeImporters}')
    sys.path.append(flpath)
    load_module = importlib.import_module(os.path.splitext(scriptname)[0])
    globals().update(vars(load_module))

    gs_prefix = parse_bucket('~{gs_bucket}')
    gs_genotype_path = os.path.join(gs_prefix, '~{gs_genotype_path}'.lstrip('/'))
    gs_output_path = os.path.join(gs_prefix, '~{gs_output_path}'.lstrip('/'))
    gs_temp_path = os.path.join(gs_prefix, '~{gs_temp_path}'.lstrip('/'))

    drc_tf = '~{drc}' == 'drc'
    sample_qc_tf = '~{qc}' == 'qc'

    pheno_dct = pheno_str_to_dict('~{phenotype_id}')
    trait_type = SAIGE_PHENO_TYPES[pheno_dct['trait_type']]
    result_dir = get_result_path(gs_output_path, '~{suffix}', '~{pop}', '~{encoding}')
    pheno_results_dir = get_pheno_output_path(result_dir, pheno_dct, '')
    results_prefix = get_results_prefix(pheno_results_dir, pheno_dct, chr)
    results_files = get_results_files(results_prefix, '~{analysis_type}')

    variant_ht = get_merged_ht_path(gs_output_path, "~{suffix}", "~{pop}", pheno_dct, '~{encoding}')
    variant_ht_tmp = get_merged_ht_path(gs_temp_path, "~{suffix}_temp", "~{pop}", pheno_dct, '~{encoding}')

    ht = load_variant_data(output_ht_path=variant_ht,
                           temp_path=variant_ht_tmp,
                           paths='~{sep="," single_test}'.split(','),
                           extension='',
                           trait_type=trait_type,
                           pheno_dict=pheno_dct,
                           null_log='~{null_log}',
                           test_logs='~{sep="," test_logs}'.split(','))

    # now importing call stats
    path_stats = get_call_stats_ht_path(gs_genotype_path, pop="~{pop}", 
                                        analysis_type="~{analysis_type}", 
                                        sample_qc=sample_qc_tf, 
                                        use_drc_pop=drc_tf, 
                                        use_array_for_variant=False)
    ht_stats = hl.read_table(path_stats)
    n_samp_vec_path = get_n_samples_per_pop_path(gs_genotype_path, 
                                                 analysis_type="~{analysis_type}", 
                                                 sample_qc=sample_qc_tf, 
                                                 use_drc_pop=drc_tf,
                                                 use_array_for_variant=False)
    per_pop_N = {x['pop']: x.N for _, x in hl.import_table(n_samp_vec_path, impute=True).to_pandas().iterrows()}
    ht = ht.annotate(AN = ht_stats[ht.row_key].call_stats.AN)

    ht_flat = ht.annotate(variant = ht.locus.contig + ':' + hl.str(ht.locus.position) + ':' + hl.str(':').join(ht.alleles),
                          chr = ht.locus.contig, pos = ht.locus.position, ref = ht.alleles[0], alt = ht.alleles[1],
                          low_confidence = (ht.AC_Allele2 < 20) | ((ht.N - ht.AC_Allele2) < 20) | (ht.AN < (per_pop_N["~{pop}"] * 2 * ~{min_call_rate})))
    ht_flat = ht_flat.key_by('variant').drop('locus', 'alleles', 'trait_type', 'phenocode', 'pheno_sex', 'modifier')
    ht_flat.export('~{output_prefix + ".tsv.bgz"}')

    single_flat = get_merged_flat_path(gs_output_path, "~{suffix}", "~{pop}", pheno_dct, '~{encoding}')
    with open('flat_path.txt', 'w') as f:
        f.write(single_flat)


    gene_ht = get_merged_ht_path(gs_output_path, "~{suffix}", "~{pop}", pheno_dct, '~{encoding}', gene_analysis=True)
    if "~{analysis_type}" == "gene":
        load_gene_data(output_ht_path=gene_ht,
                       paths='~{sep="," single_test}'.split(','),
                       extension='',
                       trait_type=trait_type,
                       pheno_dict=pheno_dct,
                       null_log='~{null_log}',
                       test_logs='~{sep="," test_logs}'.split(','))


    CODE

    >>>

    runtime {
        docker: HailDocker
        memory: mem + ' GB'
        cpu: n_cpu_merge
        disks: 'local-disk ' + disk + ' HDD'
        preemptible: 5
    }

    output {
        File single_variant_flat_file = "/cromwell_root/" + output_prefix + ".tsv.bgz"
        String single_variant_flat_path = read_string('flat_path.txt')
    }

}