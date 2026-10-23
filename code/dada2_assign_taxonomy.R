#!/usr/bin/env Rscript

'Assign taxonomy to DADA2 ASVs using a DADA2-format reference database

Usage:
  dada2_assign_taxonomy.R --ref <PATH> --seqs <PATH> --out <PATH> [--cpus <N>]

Options:
  -h --help         Show this screen
  --ref <PATH>      Path to DADA2-format taxonomy training set FASTA (may be gzipped)
  --seqs <PATH>     Path to representative sequences FASTA produced by dada2
  --out <PATH>      Output TSV file path
  --cpus <N>        Number of threads [default: 1]
' -> doc

library(dada2)
library(docopt)
library(readr)
library(tibble)
library(stringr)
message(str_glue("DADA2 version: {packageVersion('dada2')}"))

args <- docopt(doc)
cpus <- as.integer(args$cpus)

seqs <- Biostrings::readDNAStringSet(args$seqs)
asv_ids <- names(seqs)
seq_strings <- as.character(seqs)
names(seq_strings) <- asv_ids

cat("Assigning taxonomy to", length(seqs), "ASVs using:", args$ref, "\n")
taxa <- assignTaxonomy(
    seq_strings,
    args$ref,
    multithread = cpus,
    verbose = TRUE,
)

taxa |>
    as.data.frame() |>
    rownames_to_column("asv_id") |>
    write_tsv(args$out)

cat("\nTaxonomy written to:", args$out, "\n")
