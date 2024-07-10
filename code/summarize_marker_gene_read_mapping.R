#!/usr/bin/env Rscript

'Summarize read mapping to Microcystis marker genes. 
Takes a .bam file of reads mapped to marker genes as input.

Usage:
  summarize_marker_gene_read_mapping.R -i <bam> --clade-summary <file> --marker-summary <file> --prefix <prefix> --info <file> --read-counts <file> [options]

Options:
  -i bam, --input bam                   Bam file of reads mapped to markers
  -m <file>, --marker-summary <file>    Per marker read mapping summary, a .tsv file
  -c <file>, --clade-summary <file>     Read mapping summarized per clade, a .tsv file
  -p <prefix>, --prefix=prefix          Prefix to use for contig names (e.g. sample name)
  --read-counts <file>                  Read counts for the sample, used for RPKM calculation. A .tsv file.
  --info <file>                         Marker clade membership file, also includes colors, etc. A .tsv file.
  -h, --help                            Show this help screen.
' -> doc

library(docopt)
library(tidyverse)

# To process actual command line input
arguments <- docopt(doc)
#print(arguments)

# for testing interactively
#arguments <- docopt(doc, args = c(" -i data/omics/metagenomes/samp_447/microcystis_markers/samp_447--lgt__516_mapped.bam --clade-summary data/omics/metagenomes/samp_447/microcystis_markers/samp_447--lgt__516_clade_summary.tsv --marker-summary data/omics/metagenomes/samp_447/microcystis_markers/samp_447--lgt__516_summary.tsv --prefix samp_447 --info data/reference/microcystis_markers/info/20240709_groupings.tsv --read-counts data/omics/metagenomes/samp_447/reads/samp_447_read_count_fastp.tsv"))
#print(arguments)

### Script starts here ###

# Read in genome / marker clade membership for summarizing
microcystis_groups <- read_tsv(arguments$info) %>% 
  arrange(genome) #%>% 
  #write_tsv("data/reference/microcystis_markers/info/groupings.tsv")

# Get total number of reads for sample
decon_reads <- read_tsv(arguments$read_counts) %>% 
  filter(read_state == "decon_reads") %>% 
  pull(fwd_read_count)

# Get number of reads mapped to markers
bam_counts <- Rsamtools::idxstatsBam(arguments$input)

# Per sequence RPKM
rpkms <- bam_counts %>% 
  filter(seqnames != "\\*") %>% 
  mutate(sample_read_count = decon_reads,
         per_million = sample_read_count/1000000,
         rpm = mapped/per_million,
         rpkm = rpm/seqlength,
         seqnames = seqnames %>% str_remove(";.*") %>% str_replace_all("-", "_")) %>%
  left_join(microcystis_groups %>% dplyr::rename(seqnames = genome, clade = group) %>% dplyr::mutate(seqnames = str_replace_all(seqnames, "-", "_"))) %>% 
  filter(!is.na(clade)) %>% 
  arrange(seqnames) %>% 
  write_tsv(arguments$marker_summary)
  

# RPKM summarized by sequence type/group
rpkms_summarized <- rpkms %>% 
  group_by(clade) %>%
  summarise(rpkm = sum(rpkm),
            mapped_reads = sum(mapped)) %>% 
  write_tsv(arguments$clade_summary)
