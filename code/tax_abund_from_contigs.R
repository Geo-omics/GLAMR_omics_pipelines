#!/usr/bin/env Rscript 
"Usage:
  calc_contig_abund_summary.R -l <lca_path> -r <read_mapping_summary_path> -o <tsv_out> -t <taxonkit_path> -d <taxdump_path> [-c <cores>]

Options:
  -l <lca_path>                             Specify the path to the LCA file.
  -r <summary>                              Specify the path to the read mapping summary file.
  -o <tsv_out>                              Specify the path to the tsv output file.
  -c <cores>                                Number of cores to use. Default is 8.
  -t <taxonkit_path>                        Path to the taxonkit executable. Default is ~/miniconda3/envs/taxonkit/bin/taxonkit.
  -d <taxdump_path>                         Path to directory containing TaxDump
  -h --help                                 Show this screen." -> doc

# Load libraries
library(tidyverse)
library(furrr)
library(docopt)

# Increase memory available for multithreading from default (500MB to 8GB)
options(future.globals.maxSize = 8000 * 1024^2)

# For testing only
# arguments <- docopt(doc, args = "-l data/omics/metagenomes/samp_1290/samp_1290_contig_lca.tsv
#             -r data/omics/metagenomes/samp_1290/samp_1290_contig_abund.tsv
#             -o data/omics/metagenomes/samp_1290/samp_1290_contig_abund.tsv
#             -c 8
#             -t ~/miniconda3/envs/taxonkit/bin/taxonkit
#             -d data/reference/ncbi_tax")


# Parse command-line arguments
arguments <- docopt(doc)
lca_path <- arguments[["-l"]]
read_mapping_summary_path <- arguments[["-r"]]
tsv_out <- arguments[["-o"]]
cores <- as.integer(arguments[["-c"]])
if (is.na(cores)) { cores = 8 }
taxonkit_path <- arguments[["-t"]]
if (is.na(taxonkit_path)) { taxonkit_path = "~/miniconda3/envs/taxonkit/bin/taxonkit" }
taxdump_path <- arguments[["-d"]]
if (is.na(taxdump_path)) { taxdump_path = "data/reference/ncbi_tax" }


# lca_path="data/omics/metagenomes/samp_1290/samp_1290_contig_lca.tsv"
# read_mapping_summary_path="data/omics/metagenomes/samp_1290/samp_1290_contig_abund.tsv"
# tsv_out="data/omics/metagenomes/samp_1290/samp_1290_contig_abund.tsv"
#taxonkit_path = "~/miniconda3/envs/taxonkit/bin/taxonkit"
# taxdump_path = "data/reference/ncbi_tax"


# Define function for summarizing abundance
calc_contig_abund_summary <- function(lca_path, read_mapping_summary_path, tsv_out, cores = 64, taxonkit_path = "~/miniconda3/envs/taxonkit/bin/taxonkit", taxdump_path){
  
  contig_tax_cols <- c("Contig", "tax_id", "tax_rank", "tax_name", "num_frags_retained", "num_frags_tax_assigned", "num_frags_agree_w_contig_label", "support", "lineage")
  
  lca_path_parsed <- unglue::unglue_data(lca_path, "data/omics/metagenomes/{SampleID}/{SampleID}_contig_lca.tsv")
  
  contig_abund <- read_tsv(read_mapping_summary_path)
  
  contig_lca <- read_tsv(lca_path,col_names = contig_tax_cols) %>% 
    left_join(contig_abund)
  
  unique_tax <- contig_lca %>% 
    select(tax_id) %>% 
    distinct()
  
  direct_tax_abund <- contig_lca %>% 
    group_by(tax_id) %>% 
    summarize(abund_direct = sum(TPM))
  
  
  calc_tax_sum_w_desc_taxonkit <- function(id,taxonkit=taxonkit_path,taxdump=taxdump_path){
    downstream_tax <- system(str_glue('{taxonkit} list --ids {id} --indent "" --data-dir {taxdump}'),intern = TRUE)
    
    abund <- contig_lca %>% 
      filter(tax_id %in% c(id, downstream_tax)) %>% 
      summarize(abund_w_subtax = sum(TPM)) %>% 
      mutate(tax_id = id)
    
    return(abund)
  }
  
  # Start parallel workers
  plan(multisession, workers = cores)
  
  abund_summary <- future_map_dfr(unique_tax$tax_id, calc_tax_sum_w_desc_taxonkit) %>% 
    left_join(direct_tax_abund) %>% 
    relocate(tax_id)
  
  if(all(abund_summary$abund_direct == abund_summary$abund_w_subtax)){
    stop("Abund with subtax calculation likely failed, should not be identical to direct abundance.")
  }
  
  abund_summary %>% write_tsv(tsv_out)
  
  plan(sequential) # Shut down parallel workers
    
  return(abund_summary)
}


# Call function with the variables obtained from command-line
calc_contig_abund_summary(lca_path, read_mapping_summary_path, tsv_out, cores, taxonkit_path, taxdump_path)
