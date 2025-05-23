version 1.0

workflow saige_sparse_grm {

    input {

        File pop_list
        Int num_markers
        Float relatednessCutoff

        # path to gs folders
        String gs_bucket
        String gs_genotype_path

        # helper functions
        File SaigeImporters

        # docker images
        String SaigeDocker = 'rahulg603/saige:1.4.3.3'
        String HailDocker = 'rahulg603/rgupta-hail-utils:0.2.119'

        # options
        Boolean use_plink
        Boolean use_drc_pop = true
        Float min_af = 0.01
        Boolean sample_qc = true
        Int n_cpu = 64

        Boolean use_array = true
        Int? n_common
        Int? n_maf
        Int? n_mac

    }

    Array[String] pops = read_lines(pop_list)

    call get_sparse_grm_paths as get_paths {
        input:
            pops = pops,

            gs_bucket = gs_bucket,
            gs_genotype_path = gs_genotype_path,

            n_markers = num_markers,
            relatedness_cutoff = relatednessCutoff,

            min_af = min_af,
            use_drc_pop = use_drc_pop,
            sample_qc = sample_qc,
            use_plink = use_plink,

            use_array = use_array,
            n_common = n_common,
            n_maf = n_maf,
            n_mac = n_mac,

            SaigeImporters = SaigeImporters,
            HailDocker = HailDocker
    }

    Array[File] bedfile = read_json(get_paths.bed)
    Array[File] bimfile = read_json(get_paths.bim)
    Array[File] famfile = read_json(get_paths.fam)

    scatter (this_struct in zip(zip(zip(bimfile, famfile), bedfile), pops)) {

        call create_sparse_grm {
            input:
                pop = this_struct.right,

                bedfile = this_struct.left.right,
                bimfile = this_struct.left.left.left,
                famfile = this_struct.left.left.right,

                n_markers = num_markers,
                relatedness_cutoff = relatednessCutoff,
                min_af = min_af,

                SaigeImporters = SaigeImporters,
                SaigeDocker = SaigeDocker,
                n_cpu = n_cpu
        }

    }

    Array[String] paths_for_upload = flatten([get_paths.mtx, get_paths.ix])
    Array[File] files_for_upload = flatten([create_sparse_grm.mtx, create_sparse_grm.ix])

    call upload {
        input:
            files = files_for_upload,
            paths = paths_for_upload,
            HailDocker = HailDocker
    }

    Boolean wait_for_upload = upload.task_complete

    output {
    }

}


task get_sparse_grm_paths {

    input {

        Array[String] pops

        String gs_bucket
        String gs_genotype_path

        Int n_markers
        Float relatedness_cutoff

        Float min_af
        Boolean use_drc_pop
        Boolean sample_qc
        Boolean use_plink

        Boolean use_array = true
        Int? n_common
        Int? n_maf
        Int? n_mac

        File SaigeImporters
        String HailDocker

    }

    String drc = if use_drc_pop then 'drc' else 'custom'
    String qc = if sample_qc then 'qc' else 'no_qc'
    String plink = if use_plink then 'plink' else 'no_plink'

    String arr = if use_array then 'arr' else 'no_arr'
    Int n_common_this = select_first([n_common, 0])
    Int n_maf_this = select_first([n_maf, 0])
    Int n_mac_this = select_first([n_mac, 0])

    command <<<
        set -e

        python3.8 <<CODE
    import hail as hl
    import importlib
    import os, sys
    import json

    this_temp_path = '/cromwell_root/tmp/'
    hl.init(log='log.log', tmp_dir=this_temp_path)

    # import relevant objects
    flpath = os.path.dirname('~{SaigeImporters}')
    scriptname = os.path.basename('~{SaigeImporters}')
    sys.path.append(flpath)
    load_module = importlib.import_module(os.path.splitext(scriptname)[0])
    globals().update(vars(load_module))

    gs_prefix = parse_bucket('~{gs_bucket}')
    gs_genotype_path = os.path.join(gs_prefix, '~{gs_genotype_path}'.lstrip('/'))

    drc_tf = '~{drc}' == 'drc'
    sample_qc_tf = '~{qc}' == 'qc'
    plink_tf = '~{plink}' == 'plink'
    arr_tf = '~{arr}' == 'arr'
    pops_to_iter = '~{sep="," pops}'.split(',')

    # generate filenames
    bed = []
    bim = []
    fam = []
    mtx = []
    ix = []

    for pop in pops_to_iter:

        if arr_tf:
            pref = get_ld_pruned_array_data_path(gs_genotype_path, pop=pop, 
                                                 sample_qc=sample_qc_tf,
                                                 use_drc_pop=drc_tf,
                                                 use_plink=plink_tf,
                                                 af_cutoff=~{min_af}, extension='')
            bed.append(f'{pref}bed')
            bim.append(f'{pref}bim')
            fam.append(f'{pref}fam')

            mtx_p, ix_p = get_sparse_grm_path(gs_genotype_path, pop=pop,
                                            n_markers=~{n_markers},
                                            relatedness=~{relatedness_cutoff},
                                            sample_qc=sample_qc_tf,
                                            af_cutoff=~{min_af},
                                            use_drc_pop=drc_tf,
                                            use_plink=plink_tf,
                                            use_array_data='')
        else:
            out_set = get_plink_for_null_path(gs_genotype_path, pop=pop,
                                              sample_qc=sample_qc_tf,
                                              analysis_type='both',
                                              ld_pruned=True,
                                              n_common=~{n_common_this},
                                              n_maf=~{n_maf_this},
                                              n_mac=~{n_mac_this},
                                              use_drc_pop=drc_tf,
                                              use_array_for_variant=False)
            bed.append(out_set[0])
            bim.append(out_set[1])
            fam.append(out_set[2])

            mtx_p, ix_p = get_sparse_grm_path(gs_genotype_path, pop=pop,
                                              n_markers=~{n_markers},
                                              relatedness=~{relatedness_cutoff},
                                              sample_qc=sample_qc_tf,
                                              af_cutoff=~{min_af},
                                              use_drc_pop=drc_tf,
                                              use_plink=plink_tf,
                                              use_array_data='~{n_common_this}_~{n_maf_this}_~{n_mac_this}')
        
        mtx.append(mtx_p)
        ix.append(ix_p)

    # output files
    with open('bed.json', 'w') as f:
        json.dump(bed, f)
    
    with open('bim.json', 'w') as f:
        json.dump(bim, f)
    
    with open('fam.json', 'w') as f:
        json.dump(fam, f)
    
    with open('mtx.json', 'w') as f:
        json.dump(mtx, f)
    
    with open('ix.json', 'w') as f:
        json.dump(ix, f)

    CODE

    >>>

    runtime {
        docker: HailDocker
        memory: '3 GB'
        cpu: '1'
        preemptible: 5
    }

    output {
        File bed = 'bed.json'
        File bim = 'bim.json'
        File fam = 'fam.json'

        Array[String] mtx = read_json('mtx.json')
        Array[String] ix = read_json('ix.json')
    }
}


task create_sparse_grm {
    
    input {
        
        String pop

        File bedfile
        File bimfile
        File famfile
        
        Int n_markers
        Float relatedness_cutoff
        Float min_af

        File SaigeImporters
        String SaigeDocker

        Int? n_cpu

    }

    Int this_cpu = select_first([n_cpu, 8])

    command <<<
        set -e

        Rscript /usr/local/bin/createSparseGRM.R \
        "--bedFile=~{bedfile}" \
        "--famFile=~{famfile}" \
        "--bimFile=~{bimfile}" \
        "--nThreads=~{this_cpu}" \
        "--numRandomMarkerforSparseKin=~{n_markers}" \
        "--relatednessCutoff=~{relatedness_cutoff}" \
        "--outputPrefix=~{pop}" \
        "--minMAFforGRM=~{min_af}"

        ls -lh

    >>>
    
    runtime {
        docker: SaigeDocker
        cpu: this_cpu
        preemptible: 1
    }

    output {
        File mtx = pop + "_relatednessCutoff_" + relatedness_cutoff + "_" + n_markers + "_randomMarkersUsed.sparseGRM.mtx"
        File ix = pop + "_relatednessCutoff_" + relatedness_cutoff + "_" + n_markers + "_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
    }
}


task upload {

    input {

        Array[String] paths
        Array[File] files
        String HailDocker

    }

    Int disk = ceil(size(files, 'G') * 6)

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

    local_files = '~{sep="," files}'.split(',')
    target_paths = '~{sep="," paths}'.split(',')

    for fl, target in zip(local_files, target_paths):
        hl.hadoop_copy(fl, target)
    
    CODE

    >>>

    runtime {
        docker: HailDocker
        memory: '2 GB'
        cpu: '1'
        preemptible: 5
        disks: 'local-disk ' + disk + ' HDD'
    }

    output {
        Boolean task_complete = true
    }
}