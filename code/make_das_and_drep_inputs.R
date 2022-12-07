#!/usr/bin/env Rscript

'Generate input files for DAS tool and drep. Requires checkM results and bins.

Usage:
  make_das_and_drep_inputs.R --sample=SAMPLE

Options:
  --sample=SAMPLE    The sample for which bins should be standardized
  -h --help          Show this help screen.
' -> doc

library(docopt)

# To process actual command line input
arguments <- docopt(doc)
#print(arguments)

# # for testing interactively
#arguments <- docopt(doc, args = c("--sample=df6f906562653569adf6544503e1892a"))
#print(arguments)

###########################
### Script starts here ####
###########################
library(tidyverse)

####
# Read in checkM results
cm_res <- suppressMessages(read_tsv(paste0("data/omics/metagenomes/",arguments$sample,"/bins/all_raw_bins/checkm.txt"))) %>%
  mutate(sample = str_remove(`Bin Id`, "_.*"),
         binner = str_remove(`Bin Id`, "[a-z,0-9]*_") %>% str_remove("_.*"),
         bin = str_remove_all(`Bin Id`, ".*_"),
         qual = if_else(Completeness > 90 & Contamination <= 5, "High", NA_character_),
         qual = if_else((Completeness > 50 & Contamination < 10) & is.na(qual), "Medium", qual),
         qual = if_else(Completeness > 30 & Contamination < 50 & is.na(qual), "Low", qual))


####
# Make DAS tool input files

all_bins <- read_tsv(paste0("results/bins/",arguments$sample,"_bins.tsv"))

bins_to_consider <- cm_res %>% 
  mutate(qual = factor(qual, levels = c("Low", "Medium", "High"))) %>% 
  filter(!is.na(qual))

bin_contigs <- read_rds(paste0("data/omics/metagenomes/",arguments$sample,"/bins/contig_bins.rds"))

make_das_inputs <- function(sample_id, binner_name, contig_info){
  
  das_table <- contig_info %>% 
    filter(sample == sample_id, binner == binner_name,
           new_bin_name %in% bins_to_consider$`Bin Id`) %>% 
    select(contig, new_bin_name)
  
  das_dir <- glue::glue("data/omics/metagenomes/{sample_id}/bins/das_tool")
  
  if (!dir.exists(das_dir)){dir.create(das_dir)}
  
  write_tsv(das_table,glue::glue("{das_dir}/{binner_name}_contigs.tsv"),col_names = FALSE)
  
}

sample_and_binner_combos <- expand.grid(binner = unique(all_bins$binner), sample = arguments$sample)

start_time = Sys.time()
for (i in 1:nrow(sample_and_binner_combos)) {
  make_das_inputs(contig_info = bin_contigs, 
                  sample_id = sample_and_binner_combos$sample[i],
                  binner_name = sample_and_binner_combos$binner[i])
}


####
# Make drep input files
bins_for_drep <- all_bins %>% 
  filter(new_bin_name %in% bins_to_consider$`Bin Id`) %>% 
  left_join(cm_res %>% rename(new_bin_name = "Bin Id")) %>%
  filter(!is.na(qual)) %>% 
  mutate(drep_input_bin_path = str_replace(per_sample_combined_path, "all_raw_bins", "bins_for_drep"))


# Make dirs to collect bins in sample folders
bin_dirs <- paste0("data/omics/metagenomes/",arguments$sample,"/bins/bins_for_drep")
bin_dirs %>% map(., unlink, recursive = TRUE)
bin_dirs %>% map(., dir.create)


# Link standardized bins to collection dir in sample folders
suppressWarnings(file.link(bins_for_drep$per_sample_combined_path, bins_for_drep$drep_input_bin_path))

# Write out .touch saying bins are linked for snakemake
arguments$sample %>% paste0("data/omics/metagenomes/",.,"/bins/bins_for_drep/.bins_linked") %>% file.remove()
arguments$sample %>% paste0("data/omics/metagenomes/",.,"/bins/bins_for_drep/.bins_linked") %>% file.create()

#Crete genome info files to prevent rerunning checkM unncescesarily

all_genome_info <- bins_for_drep %>% 
  select(genome = "new_bin_name",
         completeness = "Completeness",
         contamination = "Contamination", sample) %>% 
  mutate(genome = paste0(genome,".fa"))

for (samp in unique(all_genome_info$sample)){
  all_genome_info %>% 
    filter(sample == samp) %>% 
    select(-sample) %>% 
    write_csv(glue::glue("data/omics/metagenomes/{samp}/bins/bins_for_drep/genome_info.csv"))
}