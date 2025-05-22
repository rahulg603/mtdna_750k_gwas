version 1.0

workflow clump_sumstats {

  input {
      File sumstats
      String pheno
      String suffix
      String p_col
      Float p_thresh
      Int window_kb

      String chr = 'chr'
      String pos = 'pos'
      String ref = 'ref'
      String alt = 'alt'
      String? conf_col
      Int n_partitions = 20
      
      Boolean exp_p
      File script = 'gs://mito-wgs-public-free/isolate_loci.py'
      File gene_file = 'gs://mito-wgs-public-free/NCBI38_ensembl.gene.loc'
      String ref_genome = 'GRCh38'

      String docker = 'rahulg603/rgupta-hail-utils:0.2.119'

      #Optional runtime arguments
      Int? preemptible_tries
      Int? cpu
  }

  call clump {
    input:
        sumstats = sumstats,
        pheno = pheno,
        suffix = suffix,
        p_col = p_col,
        p_thresh = p_thresh,
        window_kb = window_kb,
        gene_file = gene_file,

        chr = chr, 
        pos = pos, 
        ref = ref, 
        alt = alt, 
        conf_col = conf_col,
        n_partitions = n_partitions, 
        
        exp_p = exp_p,
        script = script,
        ref_genome = ref_genome,

        docker = docker,
        preemptible_tries = preemptible_tries,
        cpu = cpu
  }

  output {
    File clumped_results = clump.out
  }
}

task clump {
  input {
      File sumstats
      String pheno
      String suffix
      String p_col
      Float p_thresh
      Int window_kb
      File gene_file

      String chr
      String pos
      String ref
      String alt
      String? conf_col
      Int n_partitions = 20
      
      Boolean exp_p
      File script
      String ref_genome
      String docker

      #Optional runtime arguments
      Int? preemptible_tries
      Int? cpu
  }
  Float stat_size = size(sumstats, "GB") * 10
  Int disk_size = ceil(stat_size) + 20
  Int machine_cpu = select_first([cpu, 4])

  String conf_field_this = select_first([conf_col, ''])
  String conf_argument = if defined(conf_col) then '--conf-col ' + conf_field_this else ''
  String exponentiate_argument = if exp_p then '--exp-p' else ''

  String file_out = suffix + '_' + pheno  + "_" + p_col + '_' + "~{p_thresh}" + '_window' + "~{window_kb}" + 'kb_loci.tsv.bgz'
  
  meta {
    description: "'Clumps' using proximity"
  }

  command <<<
    set -e

    python3.8 ~{script} \
        --sumstats ~{sumstats} \
        --gene-annot ~{gene_file} \
        --output ~{file_out} \
        --p-col ~{p_col} \
        --p-thresh ~{p_thresh} \
        --window-kb ~{window_kb} \
        --ref-genome ~{ref_genome} \
        --n-threads ~{machine_cpu} \
        --chr ~{chr} \
        --pos ~{pos} \
        --ref ~{ref} \
        --alt ~{alt} \
        --n-partitions ~{n_partitions} \
        ~{exponentiate_argument} \
        ~{conf_argument}

  >>>
  runtime {
    cpu: machine_cpu
    disks: "local-disk " + disk_size + " SSD"
    docker: docker
    preemptible: select_first([preemptible_tries, 5])
  }
  output {
    File out = file_out
  }
}