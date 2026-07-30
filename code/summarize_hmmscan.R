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

# Function for pasing hmmer output 

parse_domtblout_tidy <- function(file_path) {
  lines <- readLines(file_path)
  data_lines <- lines[!grepl("^#", lines)]
  
  # Convert to tibble for easier manipulation
  tibble(line = data_lines) %>%
    mutate(
      # Split each line into exactly 16 parts (15 fixed + description)
      parts = str_split(str_trim(line), "\\s+", n = 16)
    ) %>%
    # Only keep rows that have at least 15 parts
    filter(map_int(parts, length) >= 15) %>%
    mutate(
      target_name = map_chr(parts, ~.x[1]),
      target_acc = map_chr(parts, ~.x[2]),
      query_name = map_chr(parts, ~.x[3]),
      acc = map_chr(parts, ~.x[4]),
      hmm_from = map_dbl(parts, ~as.numeric(.x[5])),
      hmm_to = map_dbl(parts, ~as.numeric(.x[6])),
      align_from = map_dbl(parts, ~as.numeric(.x[7])),
      align_to = map_dbl(parts, ~as.numeric(.x[8])),
      env_from = map_dbl(parts, ~as.numeric(.x[9])),
      env_to = map_dbl(parts, ~as.numeric(.x[10])),
      modlen = map_dbl(parts, ~as.numeric(.x[11])),
      strand = map_chr(parts, ~.x[12]),
      evalue = map_dbl(parts, ~as.numeric(.x[13])),
      score = map_dbl(parts, ~as.numeric(.x[14])),
      bias = map_dbl(parts, ~as.numeric(.x[15])),
      description = map_chr(parts, ~ifelse(length(.x) > 15 && .x[16] != "-", .x[16], NA_character_))
    ) %>%
    select(-line, -parts)
}

# the below arguments for processing from the command line
arguments <- docopt(doc)

COLUMNS <- c("target_name", "target_acc", "query_name", "acc", "hmm_from", "hmm_to", "align_from", "align_to",
             "env_from", "env_to", "modlen", "strand", "evalue", "score", "bias", "description")

file <- arguments$input

# Read and parses the domtblout file: this is a hmmer human-readable table with
# columns aligned by whitespace. The last column is text, sometimes multiple
# words also separated by space. Have to specify exactly how many columns there
# are for read.table not to make a mess (its parsing magic only looks at the
# first five rows.)
num_cols_all_rows = count.fields(file)
if (is.null(num_cols_all_rows)) {
    # no hits  / file parsing failed and would keep failing with read.table()
    # so create empty table
    table = data.frame(matrix(ncol=length(COLUMNS), nrow=0))
} else {
    table <- read.table(
        file,
        header=FALSE,
        sep="",
        comment.char="#",
        fill=TRUE,
        col.names=paste0('V', seq_len(max(num_cols_all_rows))),
    )

    if (ncol(table) > 16) {
        # some descriptions are multiple words
        # combines the description columns
        table$description <- apply(table[, 16:ncol(table)], 1, paste, collapse = " ")

        # remove the per-word columns
        table <- table[, c(1:15, ncol(table))]
    }
}

colnames(table) <- COLUMNS

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
