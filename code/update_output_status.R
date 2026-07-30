#!/usr/bin/env Rscript

setwd("~/GLAMR")
library(tidyverse)
library(DBI)
library(googlesheets4)
library(glue)

 # googlesheets4::gs4_deauth()
 googlesheets4::gs4_auth(path = ".secrets/glamr-425619-f6508150aa53.json")

output_key <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1rkLra_xXPHjYinXzY1UpC6Mg5nzybcvG3HcGQaW8Tvc/edit#gid=0",) %>%
  separate_longer_delim(sample_types,delim = ",")

pg <- DBI::dbConnect(RPostgres::Postgres(),dbname = "glamr_data", host = "cayman.earth.lsa.umich.edu", port = "5432", user = "glamr_admin", password = "glamr2023")
samples <- tbl(pg, "glamr_samples")

samp_to_check <- samples %>%
  dplyr::select(SampleID, StudyID,sample_type) %>%
  collect()

outputs_to_check <- samp_to_check %>%
  inner_join(output_key, by = c("sample_type" = "sample_types"), relationship = "many-to-many")

# Vectorized per-template glue instead of rowwise()
output_status <- outputs_to_check %>%
  mutate(project = StudyID) %>%
  group_by(example_path) %>%
  mutate(full_path = as.character(glue_data(pick(everything()), first(example_path)))) %>%
  ungroup()

# Directory-listing based existence check instead of one fs::file_exists() stat
# per file: list each unique directory once, then check membership in memory.
# Faster than individual NFS calls.
output_status <- output_status %>%
  mutate(dir = fs::path_dir(full_path))

unique_dirs <- unique(output_status$dir)
dir_listings <- setNames(
  # all = TRUE is required: several products use dotfile markers (e.g.
  # .drep_done, .done_GTDB) which fs::dir_ls() silently omits by default.
  # tryCatch is required: unlike base dir.exists(), fs::dir_exists() hard-errors
  # (rather than returning FALSE) on broken symlinks - e.g. project dirs whose
  # target was moved/deleted on the underlying NFS mount.
  lapply(unique_dirs, function(d) {
    tryCatch(
      if (fs::dir_exists(d)) as.character(fs::dir_ls(d, all = TRUE)) else character(0),
      error = function(e) character(0)
    )
  }),
  unique_dirs
)

output_status <- output_status %>%
  rowwise() %>%
  mutate(output_exists = full_path %in% dir_listings[[dir]]) %>%
  ungroup() %>%
  dplyr::select(-dir) %>%
  mutate(checked = now())

dbWriteTable(pg,"output_status",output_status,overwrite = TRUE)

outputs_wide <- output_status %>%
  distinct() %>%
  dplyr::select(-full_path, -example_path, -module, -standard_product, -checked) %>%
  pivot_wider(id_cols = SampleID:sample_type,names_from = product, values_from = output_exists)

dbWriteTable(pg,"output_status_wide",outputs_wide,overwrite = TRUE)
