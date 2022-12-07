#!/usr/bin/env Rscript

'Standardize bins for each sample

Usage:
  standardize_bins.R --sample_dir=SAMPLE

Options:
  --sample_dir=DIR    The sample for which bins should be standardized
  -h --help          Show this help screen.
' -> doc

library(docopt)

# To process actual command line input
arguments <- docopt(doc)
#print(arguments)

# # for testing interactively
arguments <- docopt(doc, args = c("--sample_dir=data/projects/2022_geomicro_JGI_CSP/metagenomes/fa88397911af2a92053e4ca8e7499894"))
#print(arguments)

###########################
### Script starts here ####
###########################



library(metacoder)
library(tidyverse)
library(furrr)

project_dir <- str_remove(arguments$sample_dir, "/metagenomes/.*") %>% unique()
project <- project_dir %>% str_remove("data/projects/")

metabat_bins <- system(paste0("ls ",arguments$sample_dir,"/bins/METABAT2/*.fa"), intern = TRUE) %>% 
  data.frame(orig_bin_path = .) %>% 
  mutate(sample = orig_bin_path %>% str_extract("(?<=metagenomes/).*(?=/bins/METABAT2/.*)"),
         orig_bin_name_w_ext = orig_bin_path %>% str_remove(".*METABAT2/"),
         orig_bin_name = orig_bin_name_w_ext %>% str_remove("\\.[a-z]*$"),
         binner = "metabat2")

concoct_bins <- system(paste0("ls ",arguments$sample_dir,"/bins/CONCOCT/output/fasta_bins/*.fa"), intern = TRUE) %>% 
  data.frame(orig_bin_path = .) %>% 
  mutate(sample = orig_bin_path %>% str_extract("(?<=metagenomes/).*(?=/bins/CONCOCT/output/fasta_bins/.*)"),
         orig_bin_name_w_ext = orig_bin_path %>% str_remove(".*fasta_bins/"),
         orig_bin_name = orig_bin_name_w_ext %>% str_remove("\\.[a-z]*$"), 
         binner = "concoct")

maxbin_bins <- system(paste0("ls ",arguments$sample_dir,"/bins/maxbin/*.fasta"), intern = TRUE) %>% 
  data.frame(orig_bin_path = .) %>% 
  mutate(sample = orig_bin_path %>% str_extract("(?<=metagenomes/).*(?=/bins/maxbin/.*)"),
         orig_bin_name_w_ext = orig_bin_path %>% str_remove(".*maxbin/"),
         orig_bin_name = orig_bin_name_w_ext %>% str_remove("\\.[a-z]*$"), 
         binner = "maxbin")

semibin_bins <- system(paste0("ls ",arguments$sample_dir,"/bins/semibin/output_bins/*.fa"), intern = TRUE) %>% 
  data.frame(orig_bin_path = .) %>% 
  mutate(sample = orig_bin_path %>% str_extract("(?<=metagenomes/).*(?=/bins/semibin/output_bins/.*)"),
         orig_bin_name_w_ext = orig_bin_path %>% str_remove(".*output_bins/"),
         orig_bin_name = orig_bin_name_w_ext %>% str_remove("\\.[a-z]*$"), 
         binner = "semibin")

VAMB_bins <- system(paste0("ls ",arguments$sample_dir,"/bins/VAMB/bins/*.fna"), intern = TRUE) %>% 
  data.frame(orig_bin_path = .) %>% 
  mutate(sample = orig_bin_path %>% str_extract("(?<=metagenomes/).*(?=/bins/VAMB/.*)"),
         orig_bin_name_w_ext = orig_bin_path %>% str_remove(".*VAMB/"),
         orig_bin_name = orig_bin_name_w_ext %>% str_remove("\\.[a-z]*$"), 
         binner = "VAMB")

metadecoder_bins <- system(paste0("ls ",arguments$sample_dir,"/bins/metadecoder/bins/*.fasta"), intern = TRUE) %>% 
  data.frame(orig_bin_path = .) %>% 
  mutate(sample = orig_bin_path %>% str_extract("(?<=metagenomes/).*(?=/bins/metadecoder/.*)"),
         orig_bin_name_w_ext = orig_bin_path %>% str_remove(".*metadecoder/bins/"),
         orig_bin_name = orig_bin_name_w_ext %>% str_remove("\\.[a-z]*$"), 
         binner = "metadecoder")


all_bins <- bind_rows(metabat_bins, concoct_bins, maxbin_bins, semibin_bins, VAMB_bins, metadecoder_bins) %>% 
  filter(sample != "coassembly") %>% 
  mutate(create_time = file.mtime(orig_bin_path)) %>% 
  arrange(create_time) %>% 
  mutate(bin_num = row_number(),
         new_bin_name = glue::glue("{sample}_{binner}_{bin_num}"),
         per_sample_combined_path = glue::glue("{arguments$sample_dir}/bins/all_raw_bins/{new_bin_name}.fa"),
         all_combined_bin_path = glue::glue("{project_dir}/metagenomes/metagenome_bins/raw_combined_bins/{new_bin_name}.fa")) %>% 
  write_tsv(paste0("results/bins/",project,"__",arguments$sample,"_bins.tsv"))


# Check that all binners produced bins for all samples
binner_completion <- all_bins %>% 
  select(sample,binner) %>% 
  distinct() %>% 
  mutate(complete = TRUE,
         n_binners = length(unique(binner))) %>% 
  group_by(sample) %>% 
  mutate(n_binners_complete = n()) %>% 
  pivot_wider(names_from = binner, values_from = complete) %>% 
  write_tsv(paste0("results/bins/",arguments$sample,"_binner_completion.tsv"))

fully_binned_samples <- binner_completion %>% 
  filter(n_binners_complete == n_binners)

contig_info_paths <- system(paste0("ls ",arguments$sample_dir,"/assembly/megahit/contigs_info.tsv"),intern = TRUE)

contig_lengths <- map_df(contig_info_paths, read_tsv) %>% 
  dplyr::rename(contig = "contig_id",
                assembler_est_cov = "approx_cov")

#path <- all_bins$orig_bin_path[1]

#bin_id <- 1
get_bin_contigs <- function(bin_id){
  
  bin_path <- all_bins %>% filter(bin_num == bin_id) %>% pull(orig_bin_path)
  
  bin <- Biostrings::readDNAStringSet(bin_path)
  
  if (length(names(bin)) > 0) {
    
    contig_bins <- data.frame(contig = names(bin), bin_num = bin_id) %>% 
      left_join(all_bins %>% select(binner, sample, bin_num, new_bin_name), by = "bin_num")
  } else {
    contig_bins <- data.frame(contig = "bad bin", bin_num = bin_id) %>% 
      left_join(all_bins %>% select(binner, sample, bin_num, new_bin_name), by = "bin_num")
  }
}

#future::plan(strategy = "multisession",workers = 4)

#bin_contigs <- future_map_dfr(all_bins$bin_num,get_bin_contigs) %>% 
#  left_join(contig_lengths)

bin_contigs <- map_dfr(all_bins$bin_num,get_bin_contigs) %>% 
  left_join(contig_lengths)


bin_contigs %>% write_rds(paste0("ls ",arguments$sample_dir,"/bins/contig_bins.rds"))

#bin_contigs <- read_rds("data/omics/metagenomes/metagenome_bins/contig_bins.rds")

bin_sizes <- bin_contigs %>%
  group_by(bin_num, new_bin_name) %>% 
  summarise(genome_length = sum(length), 
            mean_est_cov = mean(assembler_est_cov),
            median_est_cov = median(assembler_est_cov))

bad_bins_no_contigs <- bin_contigs %>% 
  select(contig, bin_num) %>% 
  filter(contig == "bad bin") %>% 
  left_join(all_bins, by = "bin_num") %>% 
  select(-contig) %>% 
  distinct()

more_bad_bins <- all_bins %>% 
  filter(str_detect(orig_bin_name, "lowDepth") | str_detect(orig_bin_name, "unbinned"))

bad_bins <- bind_rows(bad_bins_no_contigs, more_bad_bins)

good_bins <- all_bins %>% 
  filter(!bin_num %in% bad_bins$bin_num)

good_bins_complete_samples <- good_bins %>% 
  filter(sample %in% fully_binned_samples$sample)


# Make dirs to collect bins in sample folders
unique(dirname(good_bins$per_sample_combined_path)) %>% map(., dir.create)

# Delete old linked bins
system(paste0("rm ",arguments$sample_dir,"/bins/all_raw_bins/*.fa"))

# Link standardized bins to collection dir in sample folders
suppressWarnings(file.link(good_bins_complete_samples$orig_bin_path, good_bins_complete_samples$per_sample_combined_path))

# Write out .touch saying bins are linked for snakemake
unique(dirname(good_bins %>% 
                 #filter(sample %in% fully_binned_samples$sample) %>% 
                 pull( per_sample_combined_path))) %>% paste0("/.bins_linked") %>% file.create()