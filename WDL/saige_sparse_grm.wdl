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
        String SaigeDocker = 'us-docker.pkg.dev/mito-wgs/mito-wgs-docker-repo/saige:1.3.6'
        String HailDocker = 'us-docker.pkg.dev/mito-wgs/mito-wgs-docker-repo/rgupta-hail-utils:0.2.119'

        # options
        Boolean use_plink
        Boolean use_drc_ancestry_data = true
        Float min_af = 0.05
        Boolean sample_qc = true
        Int n_cpu = 64

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
            use_drc_ancestry_data = use_drc_ancestry_data,
            sample_qc = sample_qc,
            use_plink = use_plink,

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
        Boolean use_drc_ancestry_data
        Boolean sample_qc
        Boolean use_plink

        File SaigeImporters
        String HailDocker

    }

    String drc = if use_drc_ancestry_data then 'drc' else 'custom'
    String qc = if sample_qc then 'qc' else 'no_qc'
    String plink = if use_plink then 'plink' else 'no_plink'

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
    pops_to_iter = '~{sep="," pops}'.split(',')

    # generate filenames
    bed = []
    bim = []
    fam = []
    mtx = []
    ix = []

    for pop in pops:
        mtx_p, ix_p = get_sparse_grm_path(gs_genotype_path, pop=pop,
                                          n_markers=~{n_markers},
                                          relatedness=~{relatedness_cutoff},
                                          sample_qc=sample_qc_tf,
                                          af_cutoff=~{min_af},
                                          use_drc_ancestry_data=drc_tf,
                                          use_plink=plink_tf)
        pref = get_ld_pruned_array_data_path(gs_genotype_path, pop=pop, sample_qc=sample_qc_tf,
                                             use_drc_ancestry_data=drc_tf,
                                             use_plink=plink_tf,
                                             af_cutoff=~{min_af}, extension='')
        bed.append(f'{pref}bed')
        bim.append(f'{pref}bim')
        fam.append(f'{pref}fam')
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

    >>>

    runtime {
        docker: HailDocker
        memory: '4 GB'
        cpu: '2'
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

        File SaigeImporters
        String SaigeDocker

        Int? n_cpu

    }

    Int this_cpu = select_first([n_cpu, 8])

    command <<<
        set -e

        Rscript /usr/local/bin/createSparseGRM.R \
        " --bedfile=~{bedfile}" \
        " --famfile=~{famfile}" \
        " --bimfile=~{bimfile}" \
        " --nThreads=~{this_cpu}" \
        " --numRandomMarkerforSparseKin=~{n_markers}" \
        " --relatednessCutoff=~{relatedness_cutoff}" \
        " --outputPrefix=~{pop}"

        echo ~{pop}.sparseGRM.mtx > mtx.txt
        echo ~{pop}.sampleIDs.txt > ix.txt

    >>>
    
    runtime {
        docker: SaigeDocker
        cpu: this_cpu
    }

    output {
        File mtx = read_string('mtx.txt')
        File ix = read_string('ix.txt')
    }
}


task upload {

    input {

        Array[String] paths
        Array[File] files
        String HailDocker

    }

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
        memory: '4 GB'
        cpu: '2'
    }

    output {
        Boolean task_complete = true
    }
}