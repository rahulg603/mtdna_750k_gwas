#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(magrittr)
library(stringr)
library(stringi)
library(readr)
library(tibble)
library(ggplot2)
library(dplyr)
library(qqman)
library(data.table)
source('sum_stats.R')

# ARGS
if (length(args) %in% 0:13) {
  stop("At least 13 arguments must be supplied (input file, phenotype, P_col, af_col, af_restrict, low_conf_restrict, exponentiate_p, keep_x, pointsize, var_as_rsid, wid, hei, cex, hq_path).n", call.=FALSE)
} else if (length(args)==14) {
  # default output file
  args[15] = ""
  args[16] = "png"
} else if (length(args)==15) {
  args[16] = "png"
}
sumstat_file <- args[1]
pheno <- args[2]
p_col <- args[3]
af_col <- args[4]
af_restrict <- as.numeric(args[5])
low_conf_restrict <- args[6]
exponentiate_p <- as.logical(args[7])
keep_x <- as.logical(args[8])
pointsize <- as.numeric(args[9])
var_as_rsid <- as.logical(args[10]) # if true, renames variant to rsID.
wid <- as.integer(args[11])
hei <- as.integer(args[12])
cex <- as.numeric(args[13])
hq_file <- as.character(args[14])
suffix <- args[15]

#if(args[16] == 'pdf') { # assumes ppi of 72
wid <- wid / 72
hei <- hei / 72
#}

# ARGS
print(paste0('Producing manhattan plot for ',pheno,'...'))
print(paste0('p column: ', p_col))
print(paste0('AF column: ', af_col))
if(!is.na(af_restrict)) {
  print(paste0('Restricting to MAF > ', af_restrict))
}
if(low_conf_restrict != 'NA') {
  print(paste0('Filtering variants that are TRUE for ', low_conf_restrict))
} else {
  low_conf_restrict <- NA
}
if(hq_file == 'NA') {
  hq_file <- NA
} else {
  print('Restricting variants to high-quality only.')
}
if(exponentiate_p) {
  print('Exponentiating p values.')
}
if(keep_x) {
  print('Retaining X chromosome.')
}
print(paste0('Using plot pointsize: ', pointsize))
print('Using new version of plotting function.')

# LOAD SUMSTATS
if(str_detect(sumstat_file, '.bgz$')) {
  dest_file <- paste0('~/', basename(str_replace(sumstat_file,'.bgz$','.gz')))
  file.copy(sumstat_file, dest_file)
} else { dest_file <- sumstat_file }
dir.create('~/tempfold/')
print(paste0('Reading from ', dest_file, '...'))
if(str_detect(sumstat_file, '.bgz$')) {
  data_file <- fread(dest_file, tmpdir = '~/tempfold/')
} else { 
  data_file <- fread(dest_file)
}

data_file <- data_file %>% filter(!is.nan(get(p_col)))
names(data_file)

if(args[16] == 'pdf') {
  output_plotter <- function(...) pdf(family='Helvetica', ...)
} else {
  output_plotter <- function(...) png(res=300, units='in', ...)
}


# MAKE PLOTS
if (nrow(data_file) == 0){
  print(paste0(dest_file, ' has no non-NA records.'))
} else {
  print(paste0(dest_file, ' imported.'))
  if(!is.na(af_restrict)) {
    data_file <- data_file %>% filter((get(af_col) >= af_restrict) & (get(af_col) <= (1-af_restrict)))
    pheno <- paste0(pheno, '_afrestrict_', af_restrict)
  }
  if(!is.na(low_conf_restrict)) {
    data_file <- data_file %>% filter(!get(low_conf_restrict))
    pheno <- paste0(pheno, '_high_conf')
  }
  if(exponentiate_p) {
    data_file[[p_col]] <- exp(data_file[[p_col]])
  }
  if(('variant' %in% names(data_file)) & all(!c('chr','pos') %in% names(data_file))) {
    mat_split <- stringi::stri_split_fixed(data_file$variant,':',n = 3,simplify = T)
    data_file <- data_file %>% mutate(chr = mat_split[,1],
                                      pos = as.numeric(mat_split[,2]))
    append_chrpos <- FALSE
  } else if (any(!c('chr','pos') %in% names(data_file))) { 
    append_chrpos <- TRUE 
  } else {
    append_chrpos <- FALSE
  }
  table_all_sig <- data_file %>% filter(get(p_col) < 0.05)
  table_all_insig <- data_file %>% filter(get(p_col) >= 0.05) %>% sample_frac(0.1)
  joint_tab <- bind_rows(table_all_sig, table_all_insig)
  if(append_chrpos) {
    varinfo <- read_tsv('variants.tsv') %>% transmute(rsid, chr, pos)
    joint_tab <- joint_tab %>% left_join(y = varinfo%>%filter(rsid%in%joint_tab$SNP), by=c('SNP'='rsid'))
  }
  if(!'SNP' %in% names(joint_tab)) {
    if(all(c('chr','pos','ref','alt') %in% names(joint_tab))) {
      joint_tab[['SNP']] <- paste0(joint_tab[['chr']],':',
                                   joint_tab[['pos']],':',
                                   joint_tab[['ref']],':',
                                   joint_tab[['alt']])
    } else if (var_as_rsid) {
      joint_tab[['SNP']] <- joint_tab[['variant']]
    } else {
      joint_tab[['SNP']] <- NA
    }
  }
  if(!is.na(hq_file)) {
    hq_data <- read_tsv(hq_file)
    print(paste0('Found ', sum(joint_tab[['SNP']] %in% hq_data[['varid']]), ' SNPs of ', nrow(joint_tab), ' in hq table.'))
    hq_data_f <- hq_data %>% filter(high_quality)
    joint_tab <- joint_tab %>% filter(SNP %in% hq_data_f[['varid']])
    print(paste0('Retained ', nrow(joint_tab), ' SNPs after hq filtering.'))
  }
  if(str_detect(joint_tab[['chr']][1], '^chr')) {
    joint_tab[['chr']] <- str_remove(joint_tab[['chr']], '^chr')
  }
  table_all_sugg <- joint_tab %>% filter(get(p_col) < 5e-5)
  write_tsv(table_all_sugg, paste0(pheno, '_',suffix,'_suggestive.tsv'))
  final_tab <- joint_tab %>%
    mutate(chr = ifelse(chr %in% c('X', 'chrX'), 23, as.numeric(chr)), pos = as.numeric(pos),
           !!p_col := ifelse(get(p_col) < 1e-149, 1e-149, get(p_col)))
  if(keep_x) {
    final_tab <- final_tab %>% filter(chr %in% c(1:23))
    labs_chr <- c(1:22, 'X')
  } else {
    final_tab <- final_tab %>% filter(chr %in% c(1:22))
    labs_chr <- as.character(1:22)
  }
  #col <- c('#393d47','#939aa9')
  highest_p <- max(-log10(final_tab[[p_col]]), na.rm=T)
  highest_p_for_plot <- min(-log10(1e-149), highest_p * 1.1)
  print(paste0('Highest -log10 p-value is: ', highest_p))
  print(paste0('Plotting upper bound is: ', highest_p_for_plot))
  output_plotter(file=paste0(pheno,'_',suffix,'_manhattan.',args[16]), 
      width = wid, height = hei, pointsize = pointsize)
  manhattan(final_tab, chr = 'chr', bp = 'pos', snp='SNP', p = p_col, cex.lab=2, cex.axis=1.7, 
            col = c('#3B322C','#8797AF'), chrlabs=labs_chr, yaxt = "n",
            family="Helvetica", ylim=c(0, highest_p_for_plot), cex=cex)
  axis(2, cex.lab=1.25, cex.axis=1.8, las=2, mgp = c(-1, 0.75, 0))
  box(col = "black", lwd=2) 
  dev.off()
  lgc <- median(qchisq(1-data_file[[p_col]],1)) / qchisq(0.5,1)
  arr_data <- data_file %>% 
    arrange(get(p_col)) %>%
    mutate(unif = (1:nrow(.))/(nrow(.)+1))
  qqp <- arr_data %>% 
    filter(get(p_col) < 0.05) %>%
    bind_rows(arr_data%>%filter(get(p_col) >= 0.05)%>%sample_n(5000)) %>%
    ggplot() + geom_point(aes(y = -log10(get(p_col)), x = -log10(unif))) +
    geom_abline(slope=1, intercept=0) +
    xlab('Expected -log10 p') + ylab('Observed -log10 p') + theme_bw() +
    annotate(geom = 'text', x = Inf, y = -Inf, hjust = 1.1, vjust = -1, 
             label = paste0('Lambda GC: ', round(lgc, 4)))
  ggsave(filename = paste0(pheno,'_',suffix,'_qq.',args[16]))
  print(paste0(dest_file, ' finished.'))
}

# ANNOTATE GENES
thisfile <- paste0(pheno, '_',suffix,'_suggestive.tsv')
if(file.exists(thisfile)) {
  summary_statistics <- read_tsv(thisfile, guess_max = 1000000) %>%
    dplyr::rename(pval = p_col, Chr = chr, Pos = pos) %>%
    mutate(Chr = as.character(Chr)) 
} else {
  summary_statistics <- tibble()
}

table_of_genes <- read_tsv('NCBI37_ensembl.gene.loc',
                           col_names = c('ensembl_gene_id', 'chromosome_name',
                                         'start_position', 'end_position', 'strand', 
                                         'external_gene_name'), 
                           col_types = 'ccddcc')

if(nrow(summary_statistics) > 0) {
  if (nrow(summary_statistics) > 4000) {
    summary_statistics <- (summary_statistics %>% arrange(pval))[1:4000,]
  }
  summary_stat_with_gene <- sumstat_annotate_gene_local(sumstats = summary_statistics, method = 'nearest_gene', 
                                                        gene.annot = table_of_genes, apply_unique = F,
                                                        pval_filter = 5e-5)
} else { summary_stat_with_gene <- NA }

if('data.frame' %in% class(summary_stat_with_gene)) {
  write_tsv(summary_stat_with_gene, file = paste0(pheno,'_',suffix,'_suggestive_gene.txt'))
}
