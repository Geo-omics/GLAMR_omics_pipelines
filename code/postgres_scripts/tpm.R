library(tidyverse)
library(DBI)

pg <- DBI::dbConnect(RPostgres::Postgres(),dbname = "glamr_data", host = "localhost", port = "5432", user = "glamr_admin", password = "glamr2023")

pg_community <- tbl(pg, "read_mapping_LCA_summary")
pg_sample_info <- tbl(pg, "glamr_samples")
pg_read_mapping <- tbl(pg, "read_mapping_to_uniref")
pg_u100 <- tbl(pg, "uniref100_info")
pg_u90 <- tbl(pg, "uniref90_info")
pg_umrad <- tbl(pg, "umrad_uniref")
mmseqs_uniref_index <- tbl(pg, "uniref100_index_from_mmseqs_db") 
mmseqs_uniref_key <- tbl(pg, "uniref100_ids_from_mmseqs_db") 
read_counts <- tbl(pg, "read_count")


# Make TPM table
pg_read_mapping %>%
  #filter(sample == "samp_447") %>% # for testing on a single sample
  distinct() %>% 
  dplyr::select(sample, target, uniref100_id, num_seqs_aligned) %>%
  left_join(left_join(mmseqs_uniref_index, mmseqs_uniref_key) %>% dplyr::select(target = "uniref100", length)) %>%
  #left_join(read_counts %>% filter(read_state == "decon_reads" & direction == "fwd") %>% dplyr::select(sample, total_reads = "count") %>% distinct()) %>% 
  group_by(sample) %>%
  mutate(per_sample_total_seqs_aligned = sum(num_seqs_aligned),
         rpkm = (num_seqs_aligned * 1e9) / (per_sample_total_seqs_aligned * length),
         tpm_numerator = as.numeric(num_seqs_aligned) / as.numeric(length),
         tpm_denom = sum(tpm_numerator),
         tpm = 1e6 * (tpm_numerator / tpm_denom)
         ) %>%
  dplyr::select(sample, target, uniref100_id, num_seqs_aligned, tpm, rpkm) %>% 
  compute("tpm2", # Name of new table
          temporary = FALSE, # Keep table
          indexes = list("sample","uniref100_id", "target") # Create indices
          )
