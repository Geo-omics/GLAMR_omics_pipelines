#!/usr/bin/env Rscript

'Generate MaxBin2 contig abundance files

Usage:
  create_maxbin_coverage.R <coverage> [options]

Options:
  coverage    Metabat2 style coverage file
  -h --help   Show this help screen.
' -> doc

library(docopt)

# To process actual command line input
arguments <- docopt(doc)
#print(arguments)

# # for testing interactively
#arguments <- docopt(doc, args = c("data/omics/metagenomes/1f1c6a3f970cbe368b4387b4a9280d14/bins/metabat_style_contig_coverage.tsv"))
#print(arguments)

###########################
### Script starts here ####
###########################

library(tidyverse,quietly = TRUE)


#path <- "data/omics/metagenomes/1f1c6a3f970cbe368b4387b4a9280d14/bins/metabat_style_contig_coverage.tsv"

create_maxbin_coverage <- function(path){
  
  bin_dir <- path %>% str_remove("metabat_style_contig_coverage.tsv")
  
  coverage <- read_tsv(path) %>% 
    select(-ends_with("-var"), -contigLen, -totalAvgDepth)
  
  samples <- colnames(coverage)[-1]
  depths_dir <- paste0(bin_dir,"maxbin","/depths/")
  depth_files <- paste0(depths_dir,samples)
  
  dir.create(depths_dir,recursive = TRUE)
  
  depth_files <- paste0(depths_dir,samples)
  
  for (i in 1:length(depth_files)){
    coverage %>% select(contigName, samples[i]) %>% 
      write_tsv(depth_files[i],col_names = FALSE)
  }

  depth_files_list <- depth_files %>% 
    data.frame(path = .) %>% 
    write_tsv(paste0(bin_dir,"maxbin/","depths.txt"),col_names = FALSE)
  
}

create_maxbin_coverage(arguments$coverage)