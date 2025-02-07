version 1.0

import "https://personal.broadinstitute.org/rahul/saige/saige_sparse_grm.wdl" as saige_tools

workflow saige_tests {

    input {
        String pheno
        String suffix
        String pop
        String encoding

        File null_rda
        File null_var_ratio
        File sample_list
        File? sparse_grm
        File? sparse_grm_ids
        Array[Array[String]] tests

        String? groups
        String? max_maf_for_group

        Float min_mac
        Float min_maf
        Int markers_per_chunk

        # path to gs folders
        String gs_bucket
        String gs_genotype_path
        String gs_output_path
        String? google_project_req_pays

        # options
        Boolean rvas_mode
        Boolean always_use_sparse_grm
        Boolean disable_loco

        # runtime
        Int n_cpu_test

        # helper functions
        File SaigeImporters

        # docker images
        String HailDocker
        String SaigeDocker
    }

    String analysis_type = if rvas_mode then 'gene' else 'variant'

    scatter (this_chr in tests) {

        if (this_chr[1] == "") {

            File this_bgen = this_chr[5] + '.bgen'
            File this_bgi = this_chr[5] + '.bgen.bgi'
            File this_bgen_sample = this_chr[5] + '.sample'

            call run_test {
                # this function will read in a single phenotype flat file, munge them into a correct format, and output the phenotypes to process
                input:
                    phenotype_id = pheno,
                    pop = pop,
                    chr = this_chr[4],
                    segment = this_chr[0],
                    suffix = suffix,
                    encoding = encoding,
                    analysis_type = analysis_type,

                    bgen = this_bgen,
                    bgi = this_bgi,
                    bgen_sample = this_bgen_sample,
                    sparse_grm = sparse_grm,
                    sparse_grm_ids = sparse_grm_ids,

                    group_file = this_chr[6],
                    groups = groups,
                    max_maf_for_group = max_maf_for_group,

                    null_rda = null_rda,
                    null_var_ratio = null_var_ratio,
                    samples = sample_list,
                    
                    gs_bucket = gs_bucket,
                    gs_output_path = gs_output_path,

                    always_use_sparse_grm = always_use_sparse_grm,
                    disable_loco = disable_loco,

                    min_maf = min_maf,
                    min_mac = min_mac,
                    markers_per_chunk = markers_per_chunk,

                    SaigeImporters = SaigeImporters,
                    SaigeDocker = SaigeDocker,

                    n_cpu_test = n_cpu_test
            }

            if (analysis_type == 'variant') {
                call saige_tools.upload as u_variant {
                    input:
                        paths = select_all([run_test.single_test_path, run_test.log_path]),
                        files = select_all([run_test.single_test, run_test.log]),
                        HailDocker = HailDocker
                }
            }

            if (analysis_type == 'gene') {
                call saige_tools.upload as u_gene {
                    input:
                        paths = select_all([run_test.single_test_path, run_test.gene_test_path, run_test.log_path]),
                        files = select_all([run_test.single_test, run_test.gene_test, run_test.log]),
                        HailDocker = HailDocker
                }
            }
        }

        File single_test_output = select_first([run_test.single_test, this_chr[1]])
        File test_log = select_first([run_test.log, this_chr[3]])
        
        if (analysis_type == 'gene') {
            File gene_test_output = select_first([run_test.gene_test, this_chr[2]])
        }

    }

    output {
        Array[File] single_variant = single_test_output
        Array[File] test_logs = test_log
        Array[File?] gene_test = gene_test_output
    }

}

task run_test {

    input {
        String phenotype_id
        String chr
        String? segment
        String suffix
        String pop
        String encoding

        File bgen
        File bgi
        File bgen_sample
        File? sparse_grm
        File? sparse_grm_ids

        File? group_file
        String? groups
        String? max_maf_for_group

        File null_rda
        File null_var_ratio
        File samples

        String gs_bucket
        String gs_output_path

        Float min_maf
        Float min_mac
        Int markers_per_chunk

        Boolean disable_loco
        Boolean always_use_sparse_grm

        String analysis_type

        Int n_cpu_test

        File SaigeImporters
        String SaigeDocker
    }

    String tf_defined_spGRM = if defined(sparse_grm) then "defined" else "not"
    String always_spGRM = if always_use_sparse_grm then 'sp' else 'not'
    String loco = if disable_loco then "noloco" else "loco"
    String seg = select_first([segment, chr])
    String output_prefix = phenotype_id + "." + analysis_type + "." + seg + "." + pop + "." + suffix
    Int disk = ceil(size(bgen, 'G') * 3)

    command <<<
        set -e
        set -o pipefail

        tail -n +3 ~{bgen_sample} | sed 's/ / /' | awk '{print $1}' > bgen_trunc.sample

        /app/bgenix -g ~{bgen} -list | grep success | sed --expression='s/.*success,\stotal\s//g' | sed --expression='s/\svariants.//g' > /app/var_ct.txt
        
        echo "Number of variants in BGEN:"
        cat /app/var_ct.txt
        
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

    gs_prefix = parse_bucket('~{gs_bucket}')
    gs_output_path = os.path.join(gs_prefix, '~{gs_output_path}'.lstrip('/'))

    pheno_dct = pheno_str_to_dict('~{phenotype_id}')
    trait_type = SAIGE_PHENO_TYPES[pheno_dct['trait_type']]
    result_dir = get_result_path(gs_output_path, '~{suffix}', '~{pop}', '~{encoding}')
    pheno_results_dir = get_pheno_output_path(result_dir, pheno_dct, '')
    results_prefix = get_results_prefix(pheno_results_dir, pheno_dct, '~{seg}')
    results_files = get_results_files(results_prefix, '~{analysis_type}')

    saige_step_2 = ['Rscript', '/usr/local/bin/step2_SPAtests.R',
                    '--bgenFile=~{bgen}',
                    '--bgenFileIndex=~{bgi}',
                    '--chrom=~{chr}',
                    '--minMAF=~{min_maf}',
                    '--minMAC=~{min_mac}',
                    '--sampleFile=bgen_trunc.sample',
                    '--GMMATmodelFile=~{null_rda}',
                    '--varianceRatioFile=~{null_var_ratio}',
                    '--AlleleOrder=ref-first',
                    '--SAIGEOutputFile=~{output_prefix + ".result.txt"}',
                    '--markers_per_chunk=~{markers_per_chunk}']

    if "~{analysis_type}" == "variant":
        if ("~{tf_defined_spGRM}" == 'defined') and ("~{always_spGRM}" == "sp"):
            saige_step_2 = saige_step_2 + ['--sparseGRMFile=~{sparse_grm}',
                                           '--sparseGRMSampleIDFile=~{sparse_grm_ids}',
                                           '--is_fastTest=TRUE',
                                           '--LOCO=FALSE']
        else:
            if '~{loco}' == 'loco':
                saige_step_2 = saige_step_2 + ['--LOCO=TRUE']
            else:
                saige_step_2 = saige_step_2 + ['--LOCO=FALSE']
    
    if "~{analysis_type}" == 'gene':
        saige_step_2 = saige_step_2 + ['--cateVarRatioMinMACVecExclude=0.5,1.5,2.5,3.5,4.5,5.5,10.5,15.5,20.5',
                                       '--cateVarRatioMaxMACVecInclude=1.5,2.5,3.5,4.5,5.5,10.5,15.5,20.5',
                                       '--groupFile=~{group_file}',
                                       '--annotation_in_groupTest=~{groups}',
                                       '--is_output_markerList_in_groupTest=TRUE',
                                       '--maxMAF_in_groupTest=~{max_maf_for_group}',
                                       '--is_single_in_groupTest=TRUE']
        if trait_type == 'binary':
            saige_step_2 = saige_step_2 + ['--is_output_moreDetails=TRUE']

        if ("~{tf_defined_spGRM}" == 'defined'):
            saige_step_2 = saige_step_2 + ['--sparseGRMFile=~{sparse_grm}',
                                           '--sparseGRMSampleIDFile=~{sparse_grm_ids}',
                                           '--is_fastTest=TRUE',
                                           '--LOCO=FALSE']

    with open('/app/var_ct.txt') as f:
        var_count = int(f.readline().strip())
        print(f'Variant count in BGEN: {str(var_count)}')
    
    if var_count > 0:
        with open('~{output_prefix}' + ".log", "wb") as f:
            process = subprocess.Popen(saige_step_2, stdout=subprocess.PIPE)
            for c in iter(lambda: process.stdout.read(1), b""):
                sys.stdout.buffer.write(c)
                f.write(c)
    else:
        with open('~{output_prefix}' + '.log', 'w') as f:
            f.write('No variants identified in BGEN. Skipping SAIGE and generating empty files.\n')

        headers = ['CHR', 'POS', 'MarkerID', 'Allele1', 'Allele2', 'AC_Allele2',
                   'AF_Allele2', 'MissingRate', 'BETA', 'SE', 'Tstat', 'var', 'p.value',
                   'N']
        df = pd.DataFrame({x: [] for x in headers})
        df.to_csv('~{output_prefix}' + '.result.txt', sep='\t', index=None)

        if results_files[1] is not None:
            headers = ['Region', 'Group', 'max_MAF', 'Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT',
                       'BETA_Burden', 'SE_Burden', 'MAC', 'Number_rare', 'Number_ultra_rare', 'contig']

            df = pd.DataFrame({x: [] for x in headers})
            df.to_csv('~{output_prefix}' + '.geneAssoc.txt', sep='\t', index=None)
        

    with open('single_test.txt', 'w') as f:
        f.write(results_files[0])

    with open('log.txt', 'w') as f:
        f.write(results_files[2])

    with open('gene_test_path.txt', 'w') as f:
        if results_files[1] is None:
            f.write('')
        else:
            f.write(results_files[1])


    CODE

    ls -lh

    if [[ "~{analysis_type}" == 'gene' ]]; then
        input_length=$(wc -l ~{group_file} | awk '{{print $1}}')
        output_length=$(wc -l ~{output_prefix + ".result.txt"} | awk '{{print $1}}')
        echo 'Got input:' $input_length 'output:' $output_length | tee -a ~{output_prefix + '.log'}
        if [[ $input_length > 0 ]]; then 
            echo 'got input' | tee -a ~{output_prefix + '.log'}
            if [[ $output_length == 1 ]]; then 
                echo 'but not enough output' | tee -a ~{output_prefix + '.log'}
                rm -f ~{output_prefix + '.gene'} 
                exit 1
            fi
        fi
    else
        touch ~{output_prefix + ".geneAssoc.txt"}
    fi

    >>>

    runtime {
        docker: SaigeDocker
        memory: '3 GB'
        cpu: n_cpu_test
        disks: 'local-disk ' + disk + ' HDD'
        preemptible: 5
    }

    output {
        File single_test = output_prefix + '.result.txt'
        File gene_test = output_prefix + '.geneAssoc.txt'
        File log = output_prefix + ".log"
        String single_test_path = read_string('single_test.txt')
        String gene_test_path = read_string('gene_test_path.txt')
        String log_path = read_string('log.txt')
    }
}