#!/usr/bin/env Rscript

# This script creates a log of samples imported into GLAMR. Where possible, it sources information 
# from logs created at the time samples are actually imported, but also attempts to discern this information
# from sample directories for samples which were imported without expected logging.
# 
# This script is run daily via CRON on Alpena; a new log should only be produced when changes are identified.
# Crontab entry: 30 0 * * * /bin/bash -c "singularity run docker://eandersk/r_microbiome /geomicro/data2/kiledal/GLAMR/code/check_import_status.R >> /geomicro/data2/kiledal/GLAMR/logs/check_import_status/$(date +\%Y\%m\%d).log 2>&1"

message(str_glue("Started at {now()}"))

# Requires tidyverse tools and unglue package. Typically run w/ Docker image with dependencies.
suppressMessages(library(tidyverse))

setwd("/geomicro/data2/kiledal/GLAMR")

# Read in sample and study table downloaded from Google Drive
samples <- readxl::read_excel("import/Great_Lakes_Omics_Datasets.xlsx",sheet = "samples",guess_max = 3000)

message("Sample info loaded succesfully")

# Read in log built when (some) samples/files were actually imported
suppressMessages(
  import_log <- read_tsv("data/import_log.tsv", show_col_types = FALSE) %>%
    group_by(SampleID) %>% 
    filter(previously_imported == FALSE, 
           import_sucess == TRUE) %>% 
    mutate(import_time = lubridate::parse_date_time(import_time, orders = "a b d H M S y")) %>% 
    slice_max(import_time,n = 1,with_ties = FALSE)
)

message("Import log loaded succesfully")

# Build import log table
suppressMessages(imported_samples <- Sys.glob("data/omics/*/*") %>% # list sample directories
  data.frame(sample_dir = .) %>% 
  bind_cols(., .$sample_dir %>% unglue::unglue_data("data/omics/{sample_type}s/{SampleID}")) %>% # parse sample dir paths
  left_join(import_log %>% select(SampleID, import_time)) %>% 
  mutate(imported_from = if_else(file.exists(file.path(sample_dir, "reads/accession")), "SRA", "local_reads"), # imported from SRA if accession present, otherwise assume manually added
         raw_reads_dir = str_glue("data/omics/{sample_type}s/{SampleID}/reads"),
         imported_reads_fp = file.path(sample_dir, "reads/raw_fwd_reads.fastq.gz"), # path to raw reads
         decon_reads_fp = file.path(sample_dir, "reads/decon_fwd_reads_fastp.fastq.gz"), # path to decontaminated reads, important because re-downloadable raw reads are deleted after decon to save storage,
         raw_reads_fp = file.path(sample_dir, "reads/raw_fwd_reads.fastq.gz"),
         accession_fp = str_glue("data/omics/{sample_type}s/{SampleID}/reads/accession"),
         import_success = fs::file_exists(imported_reads_fp) | fs::file_exists(decon_reads_fp), # Import successful if reads present, either raw or decontam (check for both because raw deleted after QC for samples downloaded from external sources)
         import_time = if_else(is.na(import_time), # Use import time record if possible, otherwise oldest mod. time of either reads or accession file
                               min(file.info(raw_reads_fp)$mtime,
                                   file.info(accession_fp)$mtime,na.rm = TRUE),
                               import_time)) %>% 
  filter(SampleID %in% samples$SampleID) %>% 
  left_join(samples %>% select(SampleID, StudyID, ProjectID)) %>% 
  relocate(SampleID, StudyID, ProjectID, sample_type, sample_dir, 
           imported_from, import_success, import_time, imported_reads_fp) %>% 
    arrange(across(everything())) %>% 
    write_tsv(str_glue('/tmp/{format(now(), "%Y%m%d")}_sample_status.tsv')) # write log to /tmp
)

message(str_glue("Found {length(unique(imported_samples$SampleID))} sample directories, of which {length(unique(imported_samples %>% filter(import_success == TRUE) %>% pull(SampleID)))} were succesfully imported"))

# Find previous logs
suppressMessages(
  prev_logs <- Sys.glob("data/import_logs/*_sample_status.tsv") %>% 
    data.frame(path = .) %>% 
    unglue::unglue_unnest(path, "data/import_logs/{date}_sample_status.tsv",remove = FALSE) %>% 
    mutate(date = lubridate::ymd(date))
)
  
# Read in the most recent log
suppressMessages(
  most_recent_log <- prev_logs %>% 
    slice_max(order_by = date, n=1, with_ties = FALSE) %>% 
    pull(path) %>% 
    read_tsv() %>% 
    arrange(across(everything()))
)

# Read in current log to match processing w/ previous logs
suppressMessages(
  current_log <- read_tsv(str_glue('/tmp/{format(now(), "%Y%m%d")}_sample_status.tsv')) %>% 
    arrange(across(everything()))
)

# Check if current log is the same as the last exported
log_changed = all.equal(most_recent_log, current_log)

# If new log is different, move it to import log directory, otherwise delete

if (!log_changed) {
  
  fs::file_move(str_glue('/tmp/{format(now(), "%Y%m%d")}_sample_status.tsv'), 
                str_glue('data/import_logs/{format(now(), "%Y%m%d")}_sample_status.tsv'))
  
  message(str_glue("New changes detected since {max(prev_logs$date)}, writing new import log."))

  } else {
  
    fs::file_delete(str_glue('/tmp/{format(now(), "%Y%m%d")}_sample_status.tsv'))
    message(str_glue("No new changes detected since {max(prev_logs$date)}, no new log written."))
                  
}

