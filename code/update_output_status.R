#!/usr/bin/env Rscript

setwd("~/GLAMR")
library(tidyverse)
library(DBI)
library(googlesheets4)

 # googlesheets4::gs4_deauth()
 googlesheets4::gs4_auth(path = ".secrets/glamr-425619-f6508150aa53.json")

output_key <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1rkLra_xXPHjYinXzY1UpC6Mg5nzybcvG3HcGQaW8Tvc/edit#gid=0",) %>% 
  separate_longer_delim(sample_types,delim = ",")

pg <- DBI::dbConnect(RPostgres::Postgres(),dbname = "glamr_data", host = "localhost", port = "5432", user = "glamr_admin", password = "glamr2023")
samples <- tbl(pg, "glamr_samples")

samp_to_check <- samples %>% 
  dplyr::select(SampleID, StudyID,sample_type) %>% 
  collect()

outputs_to_check <- samp_to_check %>% 
  cross_join(output_key ) %>% 
  filter(sample_type == sample_types)

output_status <- outputs_to_check %>% 
  mutate(project = StudyID) %>% 
  rowwise() %>% 
  mutate(full_path = str_glue(example_path)) %>% 
  ungroup() %>% 
  mutate(output_exists = fs::file_exists(full_path),
         checked = now())

dbWriteTable(pg,"output_status",output_status,overwrite = TRUE)

outputs_wide <- output_status %>%
  distinct() %>% 
  dplyr::select(-full_path, -example_path, -module, -standard_product, -checked) %>% 
  pivot_wider(id_cols = SampleID:sample_type,names_from = product, values_from = output_exists)

dbWriteTable(pg,"output_status_wide",outputs_wide,overwrite = TRUE)

