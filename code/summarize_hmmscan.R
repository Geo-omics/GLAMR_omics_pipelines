#!/usr/bin/env Rscript

'Summarizing the tblout results of a nhmmscan

Usage: 
  summarize_hmmscan.R --input FILE [--output OUTPUT]

Options:
  -i --input=<N>        Space delimited table of domain hits
  -o --output=<N>       Output Directory to save the summary of the hmmscan. Defaults to input basename + tsv.
  -h --help             Show this screen
' -> doc

library(docopt)
library(tidyverse)

# the below arguments for processing from the command line
arguments <- docopt(doc)

file <- arguments$input

# Read and parses the domtblout file: this is a hmmer human-readable table with
# columns aligned by whitespace. The last column is text, sometimes multiple
# words also separated by space. Have to specify exactly how many columns there
# are for read.table not to make a mess (its parsing magic only looks at the
# first five rows.)
num_cols = max(count.fields(file))
table <- read.table(
    file,
    header=FALSE,
    sep="",
    comment.char="#",
    fill=TRUE,
    col.names=paste0('V', seq_len(num_cols)),
)

# combines the description columns
table$description <- apply(table[, 16:ncol(table)], 1, paste, collapse = " ")

# combines the description column with the rest of the columns and names each column
table <- table[, c(1:15, ncol(table))]
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


write_tsv(summarized_table, arguments$output %||% stdout())
