#! /usr/bin/env Rscript

'Merge bracken results from multiple samples into one table.

Usage:
  merge_bracken.R [options]

Options:
  -t FILE, --taxonomy=FILE    Full tax. lineage from merge_kraken_tax.R
  -c FILE, --counts-out=FILE  Combined count table.
  -r FILE, --rel-out=FILE     Combined relative abundance table.
  -h --help                   Show this help screen.
' -> doc

library(docopt)

# To process actual command line input
arguments <- docopt(doc)
#print(arguments)

# # for testing interactively
arguments <- docopt(doc, args = c("--taxonomy=data/reference/kraken_tax_info_merged.tsv",
                                   "--counts-out=data/sample_data/JGI_bracken_counts.tsv", 
                                   "--rel-out=data/sample_data/JGI_bracken_rel_abund.tsv"))
# print(arguments)


library(tidyverse)

setwd(here::here("~/GLAMR/")) #set working directory to project root
# List prodigal gene calls for each bin

tax_levels <- c("Domain","Phylum","Class","Order","Family","Genus","Species")

merged_tax_summary <- read_tsv(trimws(arguments$taxonomy))

#list.files(here::here(),pattern = "_bracken.txt",full.names = TRUE,recursive = TRUE) %>% 

bracken_files <- system("cd ~/GLAMR/ & find -type f -name '*_bracken.txt'",intern = TRUE) %>% 
  data.frame(file = .) %>% 
  filter(!str_detect(file,"_for_bracken.txt")) %>% 
  mutate(sample = str_remove(file,".*\\/") %>% str_remove("_bracken.txt"),
         database = str_remove(sample, "_.*"),
         sample = str_remove(sample, ".*_"))

import_file <- bracken_files$file[1]

read_bracken <- function(import_file){
  
  file_info <- bracken_files %>% 
    filter(file == import_file)
  
  bracken <- read_tsv(import_file) %>% 
    mutate(sample = file_info$sample[1],
           database = file_info$database[1],
           taxonomy_id = paste0(database,"_",taxonomy_id))
  
}


abund_mat <- map_df(bracken_files$file,read_bracken) %>% 
  left_join(merged_tax_summary) %>% 
  filter(!(database == "refseq" & str_detect(taxonomy,"d__Archaea")),
         !(database == "refseq" & str_detect(taxonomy,"d__Bacteria"))) %>% 
  group_by(sample) %>% 
  mutate(rel_abund = new_est_reads / sum(new_est_reads))

jgi_samples <- system("ls import/staging/jgi_2022/all_sample_filtered_reads/*_interleaved.fastq.gz", intern = TRUE) %>% 
  str_remove("_interleaved.fastq.gz") %>% str_remove(".*filtered_reads/")

wide_counts <- abund_mat %>% 
  ungroup() %>% 
  filter(sample %in% jgi_samples) %>% 
  select(sample, taxonomy_id, taxonomy, database, all_of(tax_levels),new_est_reads) %>% 
  pivot_wider(names_from = sample, values_from = new_est_reads,values_fill = 0) %>% 
  write_tsv(trimws(arguments$counts_out))

wide_rel_abund <- abund_mat %>% 
  ungroup() %>% 
  filter(sample %in% jgi_samples) %>% 
  select(sample, taxonomy_id, taxonomy, database, all_of(tax_levels),rel_abund) %>% 
  pivot_wider(names_from = sample, values_from = rel_abund,values_fill = 0) %>% 
  write_tsv(trimws(arguments$rel_out))
