#!/usr/bin/env Rscript

'Summarizing the tblout results of a nhmmscan and identifying the gene region

Usage: 
  summarize_and_identify_gene.R --input FILE [--output OUTPUT] [--directory DIRECTORY] [--fasta FILE]

Options:
  -i --input=<N>        Space delimited table of domain hits
  -o --output=<N>       Output Directory to save the summary of the hmmscan. Defaults to input basename + tsv.
  -d --directory=<N>    Output Directory to save the summary of the hmmscan
  -f --fasta=<N>        Fasta file associated with the input table
  -h --help             Show this screen
' -> doc

library(docopt)
library(tidyverse)

# the below arguments for processing from the command line
arguments <- docopt(doc)

# uncomment the below arguments for testing
#arguments <- docopt(doc, args = "--input data/omics/amplicons/samp_832/detect_region/fwd.txt --directory data/omics/amplicons/samp_832/detect_region/summary_fwd")

file <- arguments$input
#fs::dir_create(path = arguments$directory)

# read and parses the domtblout file
table <- read.table(file, header = FALSE, sep = "", comment.char = "#", fill = TRUE)

# combines the description columns
table$description <- apply(table[, 16:ncol(table)], 1, paste, collapse = " ")

# combines the description column with the rest of the columns and names each column
table <- table[, c(1:15, 22)]
colnames(table) <- c("target_name", "target_acc", "query_name", "acc", "hmm_from", "hmm_to", "align_from", "align_to", 
                     "env_from", "env_to", "modlen", "strand", "evalue", "score", "bias", "description")

summarized_table <- table %>%
  group_by(query_name) %>% # group by original sequence to find which model best matches
  slice_min(evalue, # Only consider model with lowest e-value
            n = 1,
            with_ties = FALSE # If two models have equal e-value (unlikely), pick one randomly to keep accurate sequence counts
  ) %>%
  dplyr::rename(hmm_model = "target_name") %>% # change column name
  group_by(hmm_model) %>%  # summarize on a per model basis
  summarise(n_seqs = n(), # Number of seqs for which this model had lowest e-value
            hmm_start_median = median(hmm_from),
            hmm_end_median = median(hmm_to),
            e_value_median = median(evalue),
            score_median = median(score),
            seq_start_median = median(align_from),
            seq_end_median = median(align_to)) %>%
  # More efficient to pull out tax_group & gene at the end where it only has to parse a few rows rather than one for each sequence
  mutate(tax_group = case_when(str_detect(hmm_model, ".*_bac") ~ "bacteria",
                               str_detect(hmm_model, ".*_arc") ~ "archaea",
                               str_detect(hmm_model, ".*_euk") ~ "eukaryote",
                               str_detect(hmm_model, ".*_mito") ~ "mitochondria",
                               .default = NA_character_), # NA when no match
         gene = case_when(str_detect(hmm_model, "16S") ~ "16S_rRNA",
                          str_detect(hmm_model, "12S") ~ "12S_rRNA",
                          str_detect(hmm_model, "18S") ~ "18S_rRNA",
                          str_detect(hmm_model, "28S") ~ "28S_rRNA",
                          str_detect(hmm_model, "23S") ~ "23S_rRNA",
                          str_detect(hmm_model, "5S") ~ "5S_rRNA",
                          str_detect(hmm_model, "5_8S") ~ "5.8S_rRNA",
                          .default = NA_character_)) %>%  # NA when no match
  arrange(desc(n_seqs)) # sort table by number of seqs

if (is.null(arguments$output)) {
  tsv_path <- stringr::str_glue("{fs::path_ext_remove(arguments$input)}_summary.tsv")
} else {
  tsv_path <- arguments$output
}

write_tsv(summarized_table,tsv_path)


# gene_regions <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9")
# 
# for (i in gene_regions) {
#   base_path <- stringr::str_glue("{arguments$directory}/16S")
#   fs::dir_create(file.path(base_path, i))
# }
# 
# for (i in gene_regions) {
#   base_path <- stringr::str_glue("{arguments$directory}/18S")
#   fs::dir_create(file.path(base_path, i))
# }
# 
# if(str_detect(summarized_table[1, "gene"], "16S")) {
#   
#   if (str_detect(summarized_table[1, "tax_group"], "mito")) {
#     fs::dir_create(stringr::str_glue("{arguments$directory}/16S/mito"))
#     fasta_output_path <- file.path(stringr::str_glue("{arguments$directory}/16S/mito"))
#     full_fasta_output_path <- file.path(fasta_output_path, basename(arguments$fasta))
#     file.copy(from = arguments$fasta, to = full_fasta_output_path, overwrite = TRUE)
#   }
#   else { # the bac and arc tax groups
#     summary_table_w_region <- summarized_table[1,] %>%
#       mutate(region = case_when(hmm_start_median < 100 ~ "V1",
#                                 hmm_start_median < 400 ~ "V2",
#                                 hmm_start_median < 500 ~ "V3",
#                                 hmm_start_median < 800 ~ "V4",
#                                 hmm_start_median < 900 ~ "V5",
#                                 hmm_start_median < 1100 ~ "V6",
#                                 hmm_start_median < 1200 ~ "V7",
#                                 hmm_start_median < 1300 ~ "V8",
#                                 .default = "V9"))
#     
#     fasta_output_path <- file.path(stringr::str_glue("{arguments$directory}/16S"), summary_table_w_region[1, "region"])
#     full_fasta_output_path <- file.path(fasta_output_path, basename(arguments$fasta))
#     file.copy(from = arguments$fasta, to = full_fasta_output_path, overwrite = TRUE)
#   }
# }
# 
# if(str_detect(summarized_table[1, "gene"], "18S")) {
#   
#   summary_table_w_region <- summarized_table[1,] %>%
#     mutate(region = case_when(hmm_start_median < 100 ~ "V1",
#                               hmm_start_median < 300 ~ "V2",
#                               hmm_start_median < 600 ~ "V3",
#                               hmm_start_median < 800 ~ "V4",
#                               hmm_start_median < 1100 ~ "V5",
#                               hmm_start_median < 1350 ~ "V6",
#                               hmm_start_median < 1500 ~ "V7",
#                               hmm_start_median < 1650 ~ "V8",
#                               .default = "V9"))
#   
#   fasta_output_path <- file.path(stringr::str_glue("{arguments$directory}/18S"), summary_table_w_region[1, "region"])
#   full_fasta_output_path <- file.path(fasta_output_path, basename(arguments$fasta))
#   file.copy(from = arguments$fasta, to = full_fasta_output_path, overwrite = TRUE)
# }
# 
# if (str_detect(summarized_table[1, "gene"], "12S")) {
#   fs::dir_create(stringr::str_glue("{arguments$directory}/12S"))
#   fasta_output_path <- file.path(stringr::str_glue("{arguments$directory}/12S"))
#   full_fasta_output_path <- file.path(fasta_output_path, basename(arguments$fasta))
#   file.copy(from = arguments$fasta, to = full_fasta_output_path, overwrite = TRUE)
# }
# 
# if (str_detect(summarized_table[1, "gene"], "5S")) {
#   fs::dir_create(stringr::str_glue("{arguments$directory}/5S"))
#   fasta_output_path <- file.path(stringr::str_glue("{arguments$directory}/5S"))
#   full_fasta_output_path <- file.path(fasta_output_path, basename(arguments$fasta))
#   file.copy(from = arguments$fasta, to = full_fasta_output_path, overwrite = TRUE)
# }
# 
# if (str_detect(summarized_table[1, "gene"], "5.8S")) {
#   fs::dir_create(stringr::str_glue("{arguments$directory}/5_8S"))
#   fasta_output_path <- file.path(stringr::str_glue("{arguments$directory}/5_8S"))
#   full_fasta_output_path <- file.path(fasta_output_path, basename(arguments$fasta))
#   file.copy(from = arguments$fasta, to = full_fasta_output_path, overwrite = TRUE)
# }
# 
# if (str_detect(summarized_table[1, "gene"], "23S")) {
#   fs::dir_create(stringr::str_glue("{arguments$directory}/23S"))
#   fasta_output_path <- file.path(stringr::str_glue("{arguments$directory}/23S"))
#   full_fasta_output_path <- file.path(fasta_output_path, basename(arguments$fasta))
#   file.copy(from = arguments$fasta, to = full_fasta_output_path, overwrite = TRUE)
# }
# 
# if (str_detect(summarized_table[1, "gene"], "28S")) {
#   fs::dir_create(stringr::str_glue("{arguments$directory}/28S"))
#   fasta_output_path <- file.path(stringr::str_glue("{arguments$directory}/28S"))
#   full_fasta_output_path <- file.path(fasta_output_path, basename(arguments$fasta))
#   file.copy(from = arguments$fasta, to = full_fasta_output_path, overwrite = TRUE)
# }