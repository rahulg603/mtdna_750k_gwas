version 1.0

workflow ManhattanPlotter {

  input {
      File sumstats
      String pop
      String pheno
      String suffix

      String p_col
      String af_col
      String? conf_col 
      # only filters if above is provided

      Int wid
      Int hei
      Float cex
      Float point_size
      Float? af_filter
      
      Boolean exponentiate_p
      Boolean var_as_rsid = true
      Boolean keep_x = false

      File? hq_file
      File plotting_script = 'gs://mito-wgs-public-free/produce_manhattan.R'

      #Optional runtime arguments
      Int? preemptible_tries
      Int? mem
  }

  call RunManhattan {
    input:
        sumstats = sumstats,
        pop = pop,
        pheno = pheno,
        suffix = suffix,

        p_col = p_col,
        af_col = af_col,
        conf_col = conf_col,

        wid = wid,
        hei = hei,
        cex = cex,
        point_size = point_size,
        af_filter = af_filter,

        exponentiate_p = exponentiate_p,
        var_as_rsid = var_as_rsid,
        keep_x = keep_x,

        hq_file = hq_file,
        plotting_script = plotting_script,
        preemptible_tries = preemptible_tries,
        mem = mem
  }

  output {
    File manhattan = RunManhattan.manhattan
    File qq = RunManhattan.qq
    File sugg = RunManhattan.sugg
    File sugg_gene = RunManhattan.sugg_gene
  }
}

task RunManhattan {
  input {
      File sumstats
      String pop
      String pheno
      String suffix

      String p_col
      String af_col
      String? conf_col 
      # only filters if above is provided

      Int wid
      Int hei
      Float cex
      Float point_size
      Float? af_filter
      
      Boolean exponentiate_p
      Boolean var_as_rsid
      Boolean keep_x

      File? hq_file
      File plotting_script

      #Optional runtime arguments
      Int? preemptible_tries
      Int? mem
  }
  Float stat_size = size(sumstats, "GB") * 10
  Int disk_size = ceil(stat_size) + 20
  Int machine_mem = select_first([mem, 20])

  String conf_field_this = select_first([conf_col, 'NA'])
  String conf_mid = if conf_field_this == "NA" then '' else '_high_conf'
  String af_field_this = select_first([af_filter, 'NA'])
  String hq_field_this = select_first([hq_file, 'NA'])
  
  String exponentiate_argument = if exponentiate_p then 'TRUE' else 'FALSE'
  String var_arg = if var_as_rsid then 'TRUE' else 'FALSE'
  String x_arg = if keep_x then 'TRUE' else 'FALSE'

  String prefix = pheno + conf_mid + "_" + pop + '_' + suffix

  String path_manhattan = prefix + "_manhattan.png"
  String path_qq = prefix + "_qq.png"
  String path_sugg = prefix + "_suggestive.tsv"
  String path_sugg_gene = prefix + "_suggestive_gene.txt"

  String til = "~"
  
  meta {
    description: "Produces nicer looking manhattan plots"
  }

  command <<<
    set -e

    cd '/~{til}/'
    tar -xf variants.tsv.tar.gz
    gunzip -c '~{sumstats}' > this_extracted_file.tsv

    Rscript '~{plotting_script}' this_extracted_file.tsv '~{pheno}' '~{p_col}' '~{af_col}' '~{af_field_this}' '~{conf_field_this}' '~{exponentiate_argument}' '~{x_arg}' '~{point_size}' '~{wid}' '~{hei}' '~{cex}' '~{var_arg}' '~{hq_field_this}' '~{pop}_~{suffix}'

    cp ~{path_manhattan} /cromwell_root/~{path_manhattan}
    cp ~{path_qq} /cromwell_root/~{path_qq}
    cp ~{path_sugg} /cromwell_root/~{path_sugg}
    cp ~{path_sugg_gene} /cromwell_root/~{path_sugg_gene}
  >>>
  runtime {
    memory: machine_mem + " GB"
    disks: "local-disk " + disk_size + " SSD"
    docker: "us-docker.pkg.dev/mito-wgs/mito-wgs-docker-repo/r_for_manhattans"
    preemptible: select_first([preemptible_tries, 5])
  }
  output {
    File manhattan = "~{path_manhattan}"
    File qq = "~{path_qq}"
    File sugg = "~{path_sugg}"
    File sugg_gene = "~{path_sugg_gene}"
  }
}