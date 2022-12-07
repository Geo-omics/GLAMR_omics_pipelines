#! /usr/bin/env Rscript

'Merge taxonomy information from GTDB and refseq Kraken2 databases.
Requires prior processing with taxonkit to obtain lineages.

Usage:
  merge_kraken2_database_taxonomy.R -g FILE -r FILE -o FILE

Options:
  -g FILE, --gtdb=FILE     GTDB Kraken2 inspection with lineages
  -r FILE, --refseq=FILE   Refseq Kraken2 inspection with lineages
  -o FILE, --out=FILE      Output combined taxonomic information
  -h --help     Show this help screen.
' -> doc

library(docopt)

# To process actual command line input
arguments <- docopt(doc)
#print(arguments)

# # for testing interactively
# arguments <- docopt(doc, args = c("--refseq=/geomicro/data2/kiledal/references/kraken_databases/refseq/inspect_w_lineage.txt",
#                                 "--gtdb=/geomicro/data2/kiledal/references/kraken_databases/gtdb_r202/inspect_w_lineage.txt",
#                                 "--out=~/references/kraken_databases/combined_tax_info.tsv"))
# print(arguments)

########## Script Starts Here ##########

library(tidyverse)

inspect_names <- c("perct_of_minimizers_mapping","sub_minimizers","direct_minimizers","level","taxonomy_id","label","taxonomy")
tax_levels <- c("Domain","Phylum","Class","Order","Family","Genus","Species")

# Import GTDB tax info
gtdb_inspect_tax_taxonkit <- read_tsv(trimws(arguments$gtdb),col_names = inspect_names) %>% 
  separate(taxonomy,tax_levels,";",remove = FALSE) %>% 
  mutate(taxonomy_id = paste0("gtdb_",taxonomy_id))

gtdb_tax_summary <- gtdb_inspect_tax_taxonkit %>% 
  select(taxonomy_id,taxonomy,taxonomy,all_of(tax_levels)) %>% 
  mutate(across(everything(), ~replace_na(.x,"")))

# Import Refseq tax info
refseq_inspect_tax_taxonkit <- read_tsv(trimws(arguments$refseq),col_names = inspect_names) %>% 
  select(-taxonomy) %>% 
  rename(taxonomy ="X8") %>% 
  mutate(taxonomy = str_replace(taxonomy,"k__","d__")) %>% 
  separate(taxonomy,tax_levels,";",remove = FALSE) %>% 
  mutate(across(all_of(tax_levels),~if_else(str_detect(.x,"[a-z]__$"),"",.x)),
         taxonomy_id = paste0("refseq_",taxonomy_id))

refseq_tax_summary <- refseq_inspect_tax_taxonkit %>% 
  select(taxonomy_id,taxonomy,taxonomy,all_of(tax_levels))

# Merge and export the combined taxonomy info. from kraken2 databases
merged_tax_summary <- bind_rows(gtdb_tax_summary,refseq_tax_summary) %>% 
  write_tsv(trimws(arguments$o))
