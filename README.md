# Genetic analysis of mtDNA across 750k individuals

This repo contains a full library used for all analyses of mtDNA in AoU. mtDNA variant calls were generated using [mtSwirl](https://github.com/rahulg603/mtSwirl).

Here, we include functions used to:
- perform variant and sample QC on WES/WGS from AoU
- munge the QC'd variants into small BGENs and other input files which allow for fast operation of SAIGE in a parallelized fashion
- run SAIGE using Cromwell across sharded BGENs, phenotypes, and ancestries all in parallel, including running additive/recessive GWAS and rare variant testing
- format all results into summary statistics and perform quality control
- perform cross-ancestry and cross-biobank meta-analysis
- visualize results with Manhattan plots
- format and summarize mtDNA phenotypes from MatrixTables produced in mtSwirl used for downstream visualization

We also include a full Cromwell scheduling and submission library used throughout the above functions for parallization on AoU researcher workbench.