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

    mt = hl.read_matrix_table(get_saige_sumstats_mt_path(gs_gwas_path, suffix_updated, '~{encoding}', gene_analysis, pop='full'))
    mt = mt.annotate_entries(all_low_conf = hl.all(mt.summary_stats.low_confidence))
    if not legacy_exponentiate_p:
        mt = mt.annotate_entries(summary_stats = mt.summary_stats.map(lambda x: x.annotate(Pvalue=-1*x.Pvalue/hl.log(10))))

    mt_meta = hl.read_matrix_table(get_saige_meta_mt_path(gs_gwas_path, suffix_updated, '~{encoding}', gene_analysis))
    mt_meta = mt_meta.select_entries(meta_analysis=mt_meta.meta_analysis[0])
    if not legacy_exponentiate_p:
        mt_meta = mt_meta.annotate_entries(meta_analysis = mt_meta.meta_analysis.annotate(Pvalue=-1*mt_meta.meta_analysis.Pvalue/hl.log(10),
                                                                                          Pvalue_het=-1*mt_meta.meta_analysis.Pvalue_het/hl.log(10)))
    
    dct = pheno_str_to_dict('~{phenotype_id}')
    trait_class = SAIGE_PHENO_TYPES[dct['trait_type']]
    for k, v in dct.items():
        mt = mt.filter_cols(mt[k] == v)
    if mt.count_cols() != 1:
        raise ValueError('ERROR: mt did not have the input phenotype.')

    if remove_low_quality_sites:
        mt = mt.filter_rows(~hl.agg.all(mt.all_low_conf))
    mt = mt.drop('all_low_conf')

    for k, v in dct.items():
        mt_meta = mt_meta.filter_cols(mt_meta[k] == v)
    use_meta = mt_meta.count_cols() == 1

    if trait_class == 'quantitative':
        meta_fields = quant_meta_fields
        fields = quant_fields
        meta_field_rename_dict = quant_meta_field_rename_dict
        field_rename_dict = quant_field_rename_dict
    elif trait_class == 'binary':
        meta_fields = binary_meta_fields
        fields = binary_fields
        meta_field_rename_dict = binary_meta_field_rename_dict
        field_rename_dict = binary_field_rename_dict

    meta_field_rename_dict = rename_dict_for_log10(meta_field_rename_dict, legacy_exponentiate_p)
    field_rename_dict = rename_dict_for_log10(field_rename_dict, legacy_exponentiate_p)

    meta_fields += ['BETA','SE','Pvalue','Pvalue_het', 'N']
    fields += ['BETA','SE','Pvalue','low_confidence']

    pop_list = sorted(mt.pheno_data.pop.collect()[0])
    annotate_dict = {}
    annotate_dict.update({'chr': mt.locus.contig,
                          'pos': hl.format('%.3e', mt.locus.position),
                          'ref': mt.alleles[0],
                          'alt': mt.alleles[1]})
    if use_meta:
        for field in meta_fields:
            field_expr = mt_meta[mt.row_key,mt.col_key].meta_analysis[field]
            annotate_dict.update({f'{meta_field_rename_dict[field]}': hl.if_else(hl.is_nan(field_expr),
                                                                    hl.str(field_expr),
                                                                    hl.format('%.3e', field_expr))})

    for field in fields:
        for pop_idx, pop in enumerate(pop_list):
            field_expr = mt.summary_stats[field][pop_idx]
            annotate_dict.update({f'{field_rename_dict[field]}_{pop}': hl.if_else(hl.is_nan(field_expr),
                                                                    hl.str(field_expr),
                                                                    hl.str(field_expr) if field=='low_confidence' else hl.format('%.3e', field_expr))})

    mt_out = mt.annotate_entries(**annotate_dict)
    mt_out = mt_out.key_rows_by().drop('locus','alleles','summary_stats')
    for x in ['gene', 'annotation', 'overall_AN']:
        if x in mt_out.row:
            mt_out = mt_out.drop(x)
    
    ht = mt_out.key_cols_by().select_cols().entries()
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