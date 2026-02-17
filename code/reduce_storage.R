library(tidyverse)

# Clean-up large intermediate files from GLAMR
setwd(here::here("~/GLAMR"))

dir_size_generic = function(dir){
  data.frame(path = dir %>% fs::dir_ls(recurse = TRUE)) %>% 
    mutate(size = fs::file_size(path)) %>% 
    pull("size") %>% 
    sum()
}

# Vectorize the function
dir_size_generic_vec <- Vectorize(dir_size_generic, "dir")

# MegaHit noNorm (GLAMR's primary assembly) intermediate contigs
intermediate_contig_dirs <- Sys.glob("data/omics/metagenomes/*/assembly/megahit_noNORM/intermediate_contigs") %>% 
  data.frame(dir_path = ., unglue::unglue_data(.,"data/omics/metagenomes/{SampleID}/assembly/megahit_noNORM/intermediate_contigs")) %>% 
  mutate(assembly_name = paste0(str_remove(dir_path,"/intermediate_contigs"),"/final.contigs.renamed.fa"),
         assembly_finished = file.exists(assembly_name),
         dir_size = dir_size_generic_vec(dir_path) %>%
           fs::as_fs_bytes()
  )

str_glue("Can clean up {sum(intermediate_contig_dirs %>% filter(assembly_finished == TRUE) %>% pull(dir_size))} of extraneous megahit assembly output")

## Delete intermediate contig directories where the completed assembly exists (as indicated by the presence of the renamed contigs file)
intermediate_contig_dirs %>% 
  filter(assembly_finished == TRUE) %>% 
  pull(dir_path) %>% 
  fs::dir_delete()


# MegaHit normalized intermediate contigs
intermediate_contig_dirs_normalized <- Sys.glob("data/omics/metagenomes/*/assembly/megahit/intermediate_contigs") %>% 
  data.frame(dir_path = ., unglue::unglue_data(.,"data/omics/metagenomes/{SampleID}/assembly/megahit_noNORM/intermediate_contigs")) %>% 
  mutate(assembly_done_file = paste0(str_remove(dir_path,"/intermediate_contigs"),"/done"),
         assembly_finished = file.exists(assembly_done_file),
         dir_size = dir_size_generic_vec(dir_path) %>%
           fs::as_fs_bytes()
  )

str_glue("Can clean up {sum(intermediate_contig_dirs_normalized %>% filter(assembly_finished == TRUE) %>% pull(dir_size))} of extraneous megahit assembly output")

## Delete intermediate contig directories where the completed assembly exists (as indicated by the presence of the renamed contigs file)
intermediate_contig_dirs_normalized %>% 
  filter(assembly_finished == TRUE) %>% 
  pull(dir_path) %>% 
  fs::dir_delete()



### MetaSpades extra outputs ###

# Function for calculating directory size
dir_size = function(dir){
  data.frame(path = c(Sys.glob(str_glue("{dir}/K*/")), 
                      Sys.glob(str_glue("{dir}/corrected/"))) %>% 
               fs::dir_ls(recurse = TRUE)) %>% 
    mutate(size = fs::file_size(path)) %>% 
    pull("size") %>% 
    sum()
}

# Vectorize the function
dir_size_vec <- Vectorize(dir_size, "dir")

spades_temp_dir_size <- Sys.glob("data/omics/*/*/assembly/*/spades.log") %>% 
  data.frame(log_path = .,
             dir_path = dirname(.)) %>% 
  unglue::unglue_unnest(log_path, "data/omics/{sample_type}/{SampleID}/assembly/{assembly_name}/spades.log",remove = FALSE) %>% 
  #slice_head(n = 3) %>% 
  mutate(assembly_name = paste0(str_remove(dir_path,"/spades.log"),"/scaffolds.fasta"),
         assembly_finished = file.exists(assembly_name),
         dir_size = dir_size_vec(dir_path) %>%
           fs::as_fs_bytes()
         )

str_glue("Can clean up {sum(spades_temp_dir_size$dir_size)} of extraneous SPAdes assembly output")

# Get a list of the actual directories to delete
spades_assemblies_to_clean_up <- spades_temp_dir_size %>% 
  filter(assembly_finished == TRUE 
         # & dir_size > 0
         ) %>% 
  mutate(delete_dirs = list(data.frame(delete_dir = c(Sys.glob(str_glue("{dir_path}/K*/")), 
                         Sys.glob(str_glue("{dir_path}/corrected/")))))) %>% 
  unnest(delete_dirs)

fast_delete_dir <- function(path){
  system(str_glue("rm -rf {path}"),wait = FALSE)
}

spades_assemblies_to_clean_up %>% dplyr::select(delete_dir) %>% distinct() %>% write_tsv("tmp/dirs_to_delete.txt",col_names = FALSE)

system("~/scripts/delete_dirs_in_file.sh tmp/dirs_to_delete.txt > ")

#walk(spades_assemblies_to_clean_up$delete_dir, fast_delete_dir,.progress = TRUE)

#fs::dir_delete(spades_assemblies_to_clean_up$delete_dir) # Actually delete the directories


# Remove raw reads in cases where they were downloaded from SRA and decontaminated reads have already been produced

## Make sure there's an accession number recorded for any sample where raw reads might be deleted
samples <- readxl::read_excel("import/Great_Lakes_Omics_Datasets.xlsx",sheet = "sequencing",guess_max = 3000) %>% 
  filter(!is.na(`SRA accession`))

## Make sure reads were downloaded (.reads_downloaded is present) and that decon reads have already been made
read_status_summary <- Sys.glob("data/omics/metagenomes/*/reads/raw_fwd_reads.fastq.gz") %>% 
  data.frame(fwd_raw_read_path = ., unglue::unglue_data(.,"data/omics/metagenomes/{SampleID}/reads/raw_fwd_reads.fastq.gz")) %>% 
  mutate(rev_raw_read_path = str_glue("{dirname(fwd_raw_read_path)}/raw_rev_reads.fastq.gz"),
         reads_downloaded = str_glue("{dirname(fwd_raw_read_path)}/.reads_downloaded") %>% fs::file_exists(),
         decontaminated_reads_present = str_glue("{dirname(fwd_raw_read_path)}/decon_fwd_reads_fastp.fastq.gz") %>% fs::file_exists(),
         read_counts_present = str_glue("{dirname(fwd_raw_read_path)}/{SampleID}_read_count_fastp.tsv") %>% fs::file_exists()) %>% 
  filter(reads_downloaded == TRUE & decontaminated_reads_present & read_counts_present & SampleID %in% samples$SeqSampleID) %>% 
  mutate(raw_fwd_size = fs::file_size(fwd_raw_read_path),
         raw_rev_size = fs::file_size(rev_raw_read_path),
         raw_read_size = raw_fwd_size + raw_rev_size)

str_glue("{sum(read_status_summary$raw_read_size)} of storage space can be reclaimed by deleting raw reads that have already been QCd and are availabe in SRA")

fs::file_delete(c(read_status_summary$fwd_raw_read_path, read_status_summary$rev_raw_read_path))

# Remove fastp reads once decontaminated reads have been produced & reads have been counted at key QC steps

fastp_status_summary <- Sys.glob("data/omics/metagenomes/*/reads/fastp_fwd_reads.fastq.gz") %>% 
  data.frame(fwd_fastp_read_path = ., unglue::unglue_data(.,"data/omics/metagenomes/{SampleID}/reads/fastp_fwd_reads.fastq.gz")) %>% 
  mutate(rev_fastp_read_path = str_glue("{dirname(fwd_fastp_read_path)}/fastp_rev_reads.fastq.gz"),
         decontaminated_reads_present = str_glue("{dirname(fwd_fastp_read_path)}/decon_fwd_reads_fastp.fastq.gz") %>% fs::file_exists(),
         reads_counted = str_glue("{dirname(fwd_fastp_read_path)}/{SampleID}_read_count_fastp.tsv") %>% fs::file_exists()) %>% 
  filter(decontaminated_reads_present & reads_counted & str_detect(SampleID,"^samp_")) %>% 
  mutate(fastp_fwd_size = fs::file_size(fwd_fastp_read_path),
         fastp_rev_size = fs::file_size(rev_fastp_read_path),
         fastp_read_size = fastp_fwd_size + fastp_rev_size,
         delete_date = lubridate::now()) %>% 
  write_tsv("data/storage_reclamation/fastp_removal.tsv",append = TRUE)

str_glue("{sum(fastp_status_summary$fastp_read_size)} of storage space can be reclaimed by deleting fastp intermediate reads")
  
fs::file_delete(c(fastp_status_summary$fwd_fastp_read_path, fastp_status_summary$rev_fastp_read_path))

# Remove phix-removed reads once decontaminated reads have been produced & reads have been counted at key QC steps

phix_fastp_status_summary <- Sys.glob("data/omics/*/*/reads/phix_fwd_reads_fastp.fastq.gz") %>% 
  data.frame(fwd_fastp_read_path = ., unglue::unglue_data(.,"data/omics/{sample_type}/{SampleID}/reads/phix_fwd_reads_fastp.fastq.gz")) %>% 
  mutate(rev_fastp_read_path = str_glue("{dirname(fwd_fastp_read_path)}/phix_rev_reads_fastp.fastq.gz"),
         decontaminated_reads_present = str_glue("{dirname(fwd_fastp_read_path)}/decon_fwd_reads_fastp.fastq.gz") %>% fs::file_exists(),
         reads_counted = str_glue("{dirname(fwd_fastp_read_path)}/{SampleID}_read_count_fastp.tsv") %>% fs::file_exists()) %>% 
  filter(decontaminated_reads_present & reads_counted & str_detect(SampleID,"^samp_")) %>% 
  mutate(fastp_fwd_size = fs::file_size(fwd_fastp_read_path),
         fastp_rev_size = fs::file_size(rev_fastp_read_path),
         fastp_read_size = fastp_fwd_size + fastp_rev_size,
         delete_date = lubridate::now()) %>% 
  write_tsv("data/storage_reclamation/phix_removal.tsv",append = TRUE)

str_glue("{sum(phix_fastp_status_summary$fastp_read_size)} of storage space can be reclaimed by deleting fastp intermediate reads")

fs::file_delete(c(phix_fastp_status_summary$fwd_fastp_read_path, phix_fastp_status_summary$rev_fastp_read_path))

# Remove kraken unclassified intermediate reads

kraken_unclassified_reads <- Sys.glob("data/omics/*/*/kraken_fastp/*.fasta") %>% 
  data.frame(path = ., unglue::unglue_data(.,"data/omics/{sample_type}/{SampleID}/kraken_fastp/{read_name}.fasta")) %>% 
  mutate(report_exists = str_glue("{dirname(path)}/{SampleID}_braken_metacodeR.pdf") %>% fs::file_exists()) %>% 
  mutate(path_size = fs::file_size(path)) 

str_glue("{sum(kraken_unclassified_reads$path_size)} of storage space can be reclaimed by deleting unclassified kraken reads")

fs::file_delete(c(kraken_unclassified_reads$path))


# Delete bbmap assembly indices
assembly_index_bbmap <- Sys.glob("data/omics/metagenomes/*/assembly/megahit_noNORM/ref") %>% 
  data.frame(index_dir_path = ., unglue::unglue_data(.,"data/omics/metagenomes/{SampleID}/assembly/megahit_noNORM/ref")) %>% 
  mutate(read_mapping_to_contig_res_present = str_glue("data/omics/metagenomes/{SampleID}/assembly/{SampleID}_READSvsCONTIGS.sam.gz") %>% fs::file_exists()) %>% 
  filter(read_mapping_to_contig_res_present) %>% 
  mutate(dir_size = dir_size_generic_vec(index_dir_path) %>% fs::as_fs_bytes())

str_glue("{sum(assembly_index_bbmap$dir_size)} of storage space can be reclaimed by deleting bbmap indices")

unlink(assembly_index_bbmap$index_dir_path,recursive = TRUE)

# Delete bbmap gene indices
gene_index_bbmap <- Sys.glob("data/omics/metagenomes/*/genes/ref") %>% 
  data.frame(index_dir_path = ., unglue::unglue_data(.,"data/omics/metagenomes/{SampleID}/genes/ref")) %>% 
  mutate(read_mapping_to_genes_present = str_glue("data/omics/metagenomes/{SampleID}/genes/{SampleID}_READSvsGENES.rpkm") %>% fs::file_exists()) %>% 
  filter(read_mapping_to_genes_present) %>% 
  mutate(dir_size = dir_size_generic_vec(index_dir_path) %>% fs::as_fs_bytes())

str_glue("{sum(gene_index_bbmap$dir_size)} of storage space can be reclaimed by deleting bbmap indices")

unlink(gene_index_bbmap$index_dir_path,recursive = TRUE)


# Delete bam files once bins are made
bam_dirs <- Sys.glob("data/omics/metagenomes/*/bins/bam") %>% 
  data.frame(bam_dir_path = ., unglue::unglue_data(.,"data/omics/metagenomes/{SampleID}/bins/bam")) %>% 
  group_by(bam_dir_path) %>% 
  mutate(mapping_completed = str_glue("data/omics/metagenomes/{SampleID}/bins/bam/.snakemake_timestamp") %>% fs::file_exists(),
         binning_completed = str_glue("data/omics/metagenomes/{SampleID}/bins/bins_for_drep/.bins_linked") %>% fs::file_exists(),
         gtdb_completed = str_glue("data/omics/metagenomes/{SampleID}/bins/.done_GTDB") %>% fs::file_exists(),
         drep_completed = str_glue("data/omics/metagenomes/{SampleID}/bins/.drep_done") %>% fs::file_exists(),
         dir_size = dir_size_generic_vec(bam_dir_path) %>% fs::as_fs_bytes()) 

completed_binning_bams <- bam_dirs %>% 
  filter(binning_completed & gtdb_completed & drep_completed) 

str_glue("{sum(completed_binning_bams$dir_size)} of storage space can be reclaimed by deleting bams from binning")

completed_binning_bams %>% 
  pull(bam_dir_path) %>% 
  unlink(recursive = TRUE)


# Delete maxbin depths
maxbinDepth_dirs <- Sys.glob("data/omics/metagenomes/*/bins/maxbin/depths") %>% 
  data.frame(maxbinDepths_dir_path = ., unglue::unglue_data(.,"data/omics/metagenomes/{SampleID}/bins/maxbin/depths")) %>% 
  group_by(maxbinDepths_dir_path) %>% 
  mutate(mapping_completed = str_glue("data/omics/metagenomes/{SampleID}/bins/bam/.snakemake_timestamp") %>% fs::file_exists(),
         binning_completed = str_glue("data/omics/metagenomes/{SampleID}/bins/bins_for_drep/.bins_linked") %>% fs::file_exists(),
         gtdb_completed = str_glue("data/omics/metagenomes/{SampleID}/bins/.done_GTDB") %>% fs::file_exists(),
         drep_completed = str_glue("data/omics/metagenomes/{SampleID}/bins/.drep_done") %>% fs::file_exists(),
         dir_size = dir_size_generic_vec(maxbinDepths_dir_path) %>% fs::as_fs_bytes()) 

completed_maxbinDepths <- maxbinDepth_dirs %>% 
  filter(binning_completed & gtdb_completed & drep_completed) 

str_glue("{sum(completed_maxbinDepths$dir_size)} of storage space can be reclaimed by deleting MaxBin depths")

completed_maxbinDepths %>% 
  pull(maxbinDepths_dir_path) %>% 
  unlink(recursive = TRUE)


# Delete checkM working data
checkm_temp_data <- Sys.glob("data/omics/metagenomes/*/bins/all_raw_bins/checkm/") %>% 
  data.frame(checkm_workdir = ., unglue::unglue_data(.,"data/omics/metagenomes/{SampleID}/bins/all_raw_bins/checkm/")) %>% 
  mutate(checkm_done = str_glue("data/omics/metagenomes/{SampleID}/bins/all_raw_bins/checkm.txt") %>% fs::file_exists()) %>% 
  filter(checkm_done) %>% 
  mutate(dir_size = dir_size_generic_vec(checkm_workdir) %>% fs::as_fs_bytes())

str_glue("{sum(checkm_temp_data$dir_size)} of storage space can be reclaimed by deleting checkM working data")

unlink(checkm_temp_data$checkm_workdir,recursive = TRUE)


# Delete VAMB working data
vamb_temp_data <- Sys.glob("data/omics/metagenomes/*/bins/VAMB/*.npz") %>% 
  data.frame(vamb_file = ., unglue::unglue_data(.,"data/omics/metagenomes/{SampleID}/bins/VAMB/{file}.npz")) %>% 
  mutate(checkm_done = str_glue("data/omics/metagenomes/{SampleID}/bins/all_raw_bins/checkm.txt") %>% fs::file_exists()) %>% 
  filter(checkm_done) %>% 
  mutate(file_size = fs::file_size(vamb_file))

str_glue("{sum(vamb_temp_data$file_size)} of storage space can be reclaimed by deleting VAMB working data")

fs::file_delete(vamb_temp_data$vamb_file)

# Delete concoct working data
concoct_temp_data <- Sys.glob("data/omics/metagenomes/*/bins/CONCOCT/cut_contigs_10K.fa") %>% 
  data.frame(concoct_fa = ., unglue::unglue_data(.,"data/omics/metagenomes/{SampleID}/bins/CONCOCT/cut_contigs_10K.fa")) %>% 
  mutate(checkm_done = str_glue("data/omics/metagenomes/{SampleID}/bins/all_raw_bins/checkm.txt") %>% fs::file_exists()) %>% 
  filter(checkm_done) %>% 
  mutate(file_size = fs::file_size(concoct_fa))

str_glue("{sum(concoct_temp_data$file_size)} of storage space can be reclaimed by deleting concoct working data")

fs::file_delete(concoct_temp_data$concoct_fa)


# Delete failed clumpify data

clumpify_partials <- Sys.glob("data/omics/*/*/reads/clump*.fastq.gz") %>% 
  data.frame(file = ., unglue::unglue_data(.,"data/omics/{sample_type}/{SampleID}/reads/{file_name}")) %>% 
  mutate(qc_done = str_glue("data/omics/{sample_type}/{SampleID}/reads/decon_fwd_reads_fastp.fastq.gz") %>% fs::file_exists()) %>% 
  #filter(qc_done) %>% 
  mutate(file_size = fs::file_size(file))

str_glue("{sum(clumpify_partials$file_size)} of storage space can be reclaimed by deleting failed clumpify runs. These may indicate data corruption, investigate samples.")

fs::file_delete(clumpify_partials$file)


# Delete import intermediate reads

import_intermediates <- Sys.glob("import/*/*.fastq.gz") %>% 
  data.frame(file = ., unglue::unglue_data(.,"import/{sample_type}/{prefix}_{samp_num}_{direction}.fastq.gz")) %>% 
  mutate(decon_done = str_glue("data/omics/{sample_type}/{prefix}_{samp_num}/reads/decon_fwd_reads_fastp.fastq.gz") %>% fs::file_exists()) %>% 
  filter(decon_done) %>% 
  mutate(file_size = fs::file_size(file))

str_glue("{sum(import_intermediates$file_size)} of storage space can be reclaimed by deleting import intermediates files")

fs::file_delete(import_intermediates$file)


# Compress kraken GTDB read classifications

kraken_reads <- Sys.glob("data/omics/*/*/kraken_fastp/*_out.txt") %>% 
  data.frame(file = .) %>% 
  mutate(file_size = fs::file_size(file)) %>% 
  bind_cols(fs::file_info(.$file)) %>% 
  arrange(desc(file_size)) %>% 
  mutate(time_since_edit = time_length(now() - modification_time, unit = "hours")) %>% 
  filter(time_since_edit > 1)

walk(kraken_reads$file, ~ system(str_glue("zstd -T24 -12 --force --rm {.x}")))





