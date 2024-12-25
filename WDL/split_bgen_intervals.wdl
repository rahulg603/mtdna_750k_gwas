version 1.0

import "https://personal.broadinstitute.org/rahul/saige/saige_sparse_grm.wdl" as saige_tools

workflow split_bgen_intervals {

    input {

        String chr
        Int start
        Int end
        String pop

        String bgen_prefix

        Boolean sample_qc
        Boolean variant_qc
        Boolean use_drc_pop
        Boolean mean_impute_missing = true
        Float call_rate_filter
        Int min_ac

        String analysis_type = 'variant'
        String encoding = 'additive'

        File repo_tarball

        Int n_cpu
        String HailDocker = 'us-docker.pkg.dev/mito-wgs/mito-wgs-docker-repo/rgupta-hail-utils:0.2.119'

    }

    call split_bgen {

        input:
            chr = chr,
            start = start,
            end = end,
            pop = pop,

            bgen_prefix = bgen_prefix,
            sample_qc = sample_qc,
            variant_qc = variant_qc,
            use_drc_pop = use_drc_pop,
            mean_impute_missing = mean_impute_missing,
            call_rate_filter = call_rate_filter,
            min_ac = min_ac,

            analysis_type = analysis_type,
            encoding = encoding,
            repo_tarball = repo_tarball,

            n_cpu = n_cpu,
            HailDocker = HailDocker

    }

    call index_bgen {

        input:
            bgen = split_bgen.bgen
    
    }

    call saige_tools.upload as u {
    
        input:
            paths = [split_bgen.bgi_path],
            files = [index_bgen.bgi],
            HailDocker = HailDocker
    
    }

    output {
        String bgen = split_bgen.bgen
        String bgi = split_bgen.bgi_path
    }

}

task split_bgen {

    input {

        String chr
        Int start
        Int end
        String pop

        String bgen_prefix

        Boolean sample_qc
        Boolean variant_qc
        Boolean use_drc_pop
        Boolean mean_impute_missing
        String use_custom_pcs = 'custom'
        Float call_rate_filter
        Int min_ac

        String analysis_type
        String encoding

        File repo_tarball

        Int n_cpu
        String HailDocker

    }

    String drc = if use_drc_pop then 'drc' else 'custom'
    String sqc = if sample_qc then 'qc' else 'no_qc'
    String vqc = if variant_qc then 'qc' else 'no_qc'
    String mean_impute_missing_tf = if mean_impute_missing then 'mean' else 'nomean'

    command <<<
        set -e
        ls -lh
        tar -xzf ~{repo_tarball}
        echo $PWD
        ls -lh
        cd ./saige_aou_wdl

        python3.8 <<CODE
    import hail as hl
    from AoU.munge_genotypes import *
    from utils.SaigeImporters import *

    this_temp_path = '/cromwell_root/tmp/'
    hl.init(
        master='local[~{n_cpu}]',
        tmp_dir=this_temp_path,
        driver_memory="highmem",
        driver_cores=~{n_cpu},
        worker_memory="highmem",
        default_reference="GRCh38",
        log='log.log'
    )

    pop = "~{pop}"
    drc_tf = '~{drc}' == 'drc'
    sample_qc_tf = '~{sqc}' == 'qc'
    variant_qc_tf = '~{vqc}' == 'qc'
    mean_impute_missing_tf = '~{mean_impute_missing}' == 'mean'

    mt = get_filtered_genotype_mt(analysis_type="~{analysis_type}", 
                                  filter_variants=variant_qc_tf, 
                                  filter_samples=sample_qc_tf, 
                                  use_drc_pop=drc_tf, 
                                  use_custom_pcs="~{use_custom_pcs}",
                                  pop=pop)

    # Filter to interval
    interval = hl.interval(start=hl.locus('~{chr}', ~{start}, reference_genome='GRCh38'), 
                           end=hl.locus('~{chr}', ~{end}, reference_genome='GRCh38'),
                           includes_start=True, includes_end=False)
    mt = hl.filter_intervals(mt, [interval])

    mt = mt.select_entries("GT")
    mt = mt.filter_rows(
        hl.agg.count_where(mt.GT.is_non_ref()) > 0
    )  # Filter to non-reference sites
    mt = mt.annotate_rows(
        rsid=mt.locus.contig + ":" + hl.str(mt.locus.position) + "_" + mt.alleles[0] + "/" + mt.alleles[1]
    )  # Annotate rsid

    call_stats_ht = get_call_rate_filtered_variants(pop=pop,
                                                    analysis_type="~{analysis_type}",
                                                    sample_qc=sample_qc_tf,
                                                    use_array_for_variant=False,
                                                    use_drc_pop=drc_tf,
                                                    min_call_rate=~{call_rate_filter},
                                                    ac_filter_override=~{min_ac})
    mt = mt.filter_rows(hl.is_defined(call_stats_ht[mt.row_key]))

    mt = mt.annotate_entries(
        GT=hl.if_else(mt.GT.is_haploid(), hl.call(mt.GT[0], mt.GT[0]), mt.GT)
    )
    if encoding == 'recessive':
        mt = mt.annotate_entries(GP = hl.if_else(mt.GT.is_het(), hl.call(0, 0, phased=False), mt.GT))
    mt = gt_to_gp(mt)
    mt = impute_missing_gp(mt, mean_impute=mean_impute_missing_tf)
    hl.export_bgen(mt, "~{bgen_prefix}", gp=mt.GP, varid=mt.rsid)

    print(hl.utils.range_table(10)._force_count())

    CODE
    >>>

    runtime {
        docker: HailDocker
        cpu: n_cpu
        preemptible: 5
    }

    output {
        String bgen = bgen_prefix + '.bgen'
        String bgi_path = bgen_prefix + '.bgen.bgi'
    }

}

task index_bgen {

    input {
        File bgen
    }

    command <<<
        set -e
        bgenix -index -g ~{bgen}
    >>>

    runtime {
        docker: 'befh/bgen:latest'
        cpu: 1
        preemptible: 5
    }

    output {
        File bgi = bgen + '.bgi'
    }

}