version 1.0

workflow export_single_saige_sumstats {

    input {

        String phenotype_id
        String output_path

        String encoding = 'additive'
        String suffix

        Boolean gene_analysis
        Boolean use_drc_pop
        String use_custom_pcs
        Boolean legacy_exponentiate_p = true
        Boolean remove_low_quality_sites = true
        
        String gs_bucket
        String gs_gwas_path

        # helper functions
        File SaigeImporters

        # docker images
        String HailDocker = 'us-docker.pkg.dev/mito-wgs/mito-wgs-docker-repo/rgupta-hail-utils:0.2.119'

        Int n_cpu = 32
        Int mem = 30

    }

    call run_export {
        input:
            encoding = encoding,
            suffix = suffix,
            phenotype_id = phenotype_id,
            output_path = output_path,

            gene_analysis = gene_analysis,
            use_drc_pop = use_drc_pop,
            use_custom_pcs = use_custom_pcs,
            legacy_exponentiate_p = legacy_exponentiate_p,
            remove_low_quality_sites = remove_low_quality_sites,

            gs_bucket = gs_bucket,
            gs_gwas_path = gs_gwas_path,

            SaigeImporters = SaigeImporters,

            HailDocker = HailDocker,
            n_cpu = n_cpu,
            mem = mem
    }

    output {
        Boolean task_complete = run_export.task_complete
    }

}

task run_export {

    input {
        
        String encoding
        String suffix
        String phenotype_id
        String output_path

        Boolean gene_analysis
        Boolean use_drc_pop
        String use_custom_pcs
        Boolean legacy_exponentiate_p
        Boolean remove_low_quality_sites

        String gs_bucket
        String gs_gwas_path

        # helper functions
        File SaigeImporters

        String HailDocker
        Int n_cpu
        Int mem

    }

    String drc = if use_drc_pop then 'drc' else 'custom'
    String gene = if gene_analysis then 'gene' else 'variant'
    String exp = if legacy_exponentiate_p then 'exp' else 'neglog10'
    String rm_low = if remove_low_quality_sites then 'low' else 'not_low'

    command <<<

        set -e

        python3.8 <<CODE
    import hail as hl
    import importlib
    import os, sys
    from datetime import date

    curdate = date.today().strftime("%y%m%d")

    this_temp_path = '/cromwell_root/tmp/'
    hl.init(master=f'local[~{n_cpu}]',
            log='load_results.log', tmp_dir=this_temp_path,
            spark_conf={'spark.driver.memory': '~{mem}g'})

    # import relevant objects
    flpath = os.path.dirname('~{SaigeImporters}')
    scriptname = os.path.basename('~{SaigeImporters}')
    sys.path.append(flpath)
    load_module = importlib.import_module(os.path.splitext(scriptname)[0])
    globals().update(vars(load_module))

    gs_prefix = parse_bucket('~{gs_bucket}')
    gs_gwas_path = os.path.join(gs_prefix, '~{gs_gwas_path}'.lstrip('/'))

    drc_tf = '~{drc}' == 'drc'
    suffix_updated = update_suffix('~{suffix}', drc_tf, '~{use_custom_pcs}')
    gene_analysis = '~{gene}' == 'gene'
    legacy_exponentiate_p = '~{exp}' == 'exp'
    remove_low_quality_sites = '~{rm_low}' == 'low'

    print('Exporting ~{phenotype_id}.')

    mt = hl.read_matrix_table(get_saige_cross_biobank_meta_mt_path(gs_gwas_path, suffix_updated, '~{encoding}', gene_analysis))
    mt = mt.select_entries(meta_analysis=mt.meta_analysis[0]).select_rows()
    if not legacy_exponentiate_p:
        mt = mt.annotate_entries(meta_analysis = mt.meta_analysis.annotate(Pvalue=-1*mt.meta_analysis.Pvalue/hl.log(10),
                                                                           Pvalue_het=-1*mt.meta_analysis.Pvalue_het/hl.log(10)))
    
    dct = pheno_str_to_dict('~{phenotype_id}')
    trait_class = SAIGE_PHENO_TYPES[dct['trait_type']]
    for k, v in dct.items():
        mt = mt.filter_cols(mt[k] == v)
    if mt.count_cols() != 1:
        raise ValueError('ERROR: mt did not have the input phenotype.')

    if trait_class == 'quantitative':
        meta_field_rename_dict = quant_meta_field_rename_dict
        field_rename_dict = quant_field_rename_dict
        print('Using quantitative trait mode.')
    elif trait_class == 'binary':
        meta_field_rename_dict = binary_meta_field_rename_dict
        field_rename_dict = binary_field_rename_dict

    meta_field_rename_dict = rename_dict_for_log10(meta_field_rename_dict, legacy_exponentiate_p)

    meta_fields += ['BETA','SE','Pvalue','Pvalue_het', 'N']

    pop_list = sorted(mt.pheno_data.pop.collect()[0])
    annotate_dict = {}
    annotate_dict.update({'chr': mt.locus.contig,
                          'pos': mt.locus.position,
                          'ref': mt.alleles[0],
                          'alt': mt.alleles[1]})

    print('Exporting meta analysis...')
    for field in meta_fields:
        field_expr = mt.meta_analysis[field]
        annotate_dict.update({f'{meta_field_rename_dict[field]}': hl.if_else(hl.is_nan(field_expr),
                                                                  hl.str(field_expr),
                                                                  hl.str(field_expr) if field=='N' else hl.format('%.3e', field_expr))})
    
    mt_out = mt.annotate_entries(**annotate_dict)
    mt_out = mt_out.key_rows_by().drop('locus','alleles','meta_analysis')
    for x in ['gene', 'annotation', 'overall_AN']:
        if x in mt_out.row:
            mt_out = mt_out.drop(x)
    
    ht = mt_out.key_cols_by().select_cols().entries()
    
    print('Table schema prior to export:')
    ht.describe()
    ht.export('~{output_path}')

    CODE

    >>>

    runtime {
        docker: HailDocker
        memory: mem + ' GB'
        cpu: n_cpu
        preemptible: 5
        disks: 'local-disk ' + '10' + ' HDD'
    }

    output {
        Boolean task_complete = true
    }

}