version 1.0

import "https://personal.broadinstitute.org/rahul/saige/saige_sparse_grm.wdl" as saige_tools

workflow ld_prune {

    input {

        File bedfile
        File bimfile
        File famfile
        String gs_output_path_root

        Int step_size
        Int window_size
        Float r2

        # helper functions
        File SaigeImporters

        # docker images
        String PlinkDocker = 'rahulg603/plink:1.9'
        String HailDocker = 'rahulg603/rgupta-hail-utils:0.2.119'

        Int n_cpu = 16

    }

    call plink_prune {
        input:
            bedfile = bedfile,
            bimfile = bimfile,
            famfile = famfile,

            window_size = window_size,
            r2 = r2,
            step_size = step_size,

            gs_output_path_root = gs_output_path_root,

            PlinkDocker = PlinkDocker,
            n_cpu = n_cpu
    }

    Array[String] paths_for_upload = [gs_output_path_root + '.txt', gs_output_path_root + '.log']
    Array[File] files_for_upload = [plink_prune.out, plink_prune.log]

    call saige_tools.upload {
        input:
            files = files_for_upload,
            paths = paths_for_upload,
            HailDocker = HailDocker
    }

    Boolean wait_for_upload = upload.task_complete

    output {
    }

}


task plink_prune {

    input {

        File bedfile
        File bimfile
        File famfile

        Int window_size
        Float r2
        Int step_size

        Int n_cpu

        String gs_output_path_root

        String PlinkDocker

    }

    command <<<
        set -e

        plink1.9 \
        --bed ~{bedfile} \
        --bim ~{bimfile} \
        --fam ~{famfile} \
        --indep-pairwise ~{window_size}kb ~{step_size} ~{r2}

    >>>

    runtime {
        docker: PlinkDocker
        memory: '12 GB'
        cpu: n_cpu
    }

    output {
        File out = 'plink.prune.in'
        File log = 'plink.log'
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
        preemptible: 5
    }

    output {
        Boolean task_complete = true
    }
}