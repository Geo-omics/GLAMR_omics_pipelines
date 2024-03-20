#!/usr/bin/env Rscript 

# Load packages
library(tidyverse)
library(googlesheets4)
library(googledrive)
library(DBI)

# Open database connection
pg <- DBI::dbConnect(RPostgres::Postgres(),dbname = "glamr_data", host = "localhost", port = "5432", user = "glamr_admin", password = "glamr2023")

mmseqs_tpm <- tbl(pg, "tpm2") # create pointer to postgres table with TPM values

get_pg_table_samples <- function(db, table, sampleid_col){
  # Makes use of the SQL index, *much* faster to run
  query <- sql(str_glue("WITH RECURSIVE t AS (
     SELECT min(\"{sampleid_col}\") AS sample FROM \"{table}\"
     UNION ALL
     SELECT (SELECT min(\"{sampleid_col}\") FROM \"{table}\" WHERE \"{sampleid_col}\" > t.sample)
     FROM t WHERE t.sample IS NOT NULL
     )
  SELECT sample FROM t WHERE sample IS NOT NULL
  UNION ALL
  SELECT null WHERE EXISTS(SELECT 1 FROM \"{table}\" WHERE \"{sampleid_col}\" IS NULL);"))
  
  mapped_samples <- DBI::dbGetQuery(db, query)
  return(mapped_samples)
}

tpm_samples <- get_pg_table_samples(pg, "tpm2", "sample")

existing_tpm_exports <- Sys.glob("data/omics/metagenomes/*/*_tophit_TPM.tsv") %>% 
  data.frame(path = .) %>% 
  unglue::unglue_unnest(path, "data/omics/metagenomes/{sample}/{sample2}_tophit_TPM.tsv",remove = FALSE)

samples_to_export <- tpm_samples %>% 
  filter(!sample %in% existing_tpm_exports$sample)

write_tpm_table <- function(SampleID){
  
  message(str_glue("writing TPM values for {SampleID} to: data/omics/metagenomes/{SampleID}/{SampleID}_tophit_TPM.tsv"))
  
  mmseqs_tpm %>%
    filter(sample == SampleID) %>%
    dplyr::select(-target) %>%
    arrange(desc(tpm)) %>%
    as.data.frame() %>%
    write_tsv(str_glue("data/omics/metagenomes/{SampleID}/{SampleID}_tophit_TPM.tsv"))
}

walk(samples_to_export$sample,write_tpm_table)
