#!/usr/bin/env Rscript

'Assign taxonomy to DADA2 ASVs using a DADA2-format reference database

Usage:
  dada2_assign_taxonomy.R --ref <PATH> --seqs <PATH> --out <PATH> [--cpus <N>] [--tax-levels <LIST>] [--no-rc]

Options:
  -h --help             Show this screen
  --ref <PATH>          Path to DADA2-format taxonomy training set FASTA (may be gzipped)
  --seqs <PATH>         Path to representative sequences FASTA produced by dada2
  --out <PATH>          Output TSV file path
  --cpus <N>            Number of threads [default: 1]
  --tax-levels <LIST>   Comma-separated rank names, one per rank in the reference
                        (e.g. PR2 has 9) [default: Kingdom,Phylum,Class,Order,Family,Genus,Species]
  --no-rc               Only classify ASVs in the given orientation (by default
                        the reverse complement is also tried)
' -> doc

library(dada2)
library(docopt)
library(readr)
library(tibble)
library(stringr)
message(str_glue("DADA2 version: {packageVersion('dada2')}"))

args <- docopt(doc)
cpus <- as.integer(args$cpus)
tax_levels <- str_split_1(args$tax_levels, ",")
try_rc <- !args$no_rc

seqs <- Biostrings::readDNAStringSet(args$seqs)
asv_ids <- names(seqs)
seq_strings <- as.character(seqs)
names(seq_strings) <- asv_ids

cat("Assigning taxonomy to", length(seqs), "ASVs using:", args$ref, if (try_rc) "(also trying reverse complement)", "\n")
taxa <- assignTaxonomy(
    seq_strings,
    args$ref,
    taxLevels = tax_levels,
    tryRC = try_rc,
    multithread = cpus,
    verbose = TRUE,
)

taxa |>
    as.data.frame() |>
    rownames_to_column("asv_id") |>
    write_tsv(args$out)

cat("\nTaxonomy written to:", args$out, "\n")
