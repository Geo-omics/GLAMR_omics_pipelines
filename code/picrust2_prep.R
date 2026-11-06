#!/usr/bin/env Rscript

'Prepare DADA2 ASV outputs as PICRUSt2 input

Strips descriptions from FASTA headers so IDs match the abundance table and
transposes the sample x ASV table into the ASV x sample layout PICRUSt2 expects.
ASVs with zero total abundance are dropped from both outputs.

Usage:
  picrust2_prep.R --seqs <PATH> --abund <PATH> --out-seqs <PATH> --out-abund <PATH>

Options:
  -h --help             Show this screen
  --seqs <PATH>         ASV FASTA (asvs.fasta from amplicon_asv_check)
  --abund <PATH>        DADA2 sample x ASV abundance table (asv_table.tsv)
  --out-seqs <PATH>     Output FASTA for PICRUSt2
  --out-abund <PATH>    Output ASV x sample TSV for PICRUSt2
' -> doc

library(docopt)
library(tidyverse)

args <- docopt(doc)

seqs <- Biostrings::readDNAStringSet(args$seqs)
names(seqs) <- word(names(seqs), 1)

abund <- read_tsv(args$abund, col_types = cols(sample = "c", .default = "d")) |>
    pivot_longer(-sample, names_to = "asv_id", values_to = "count") |>
    group_by(asv_id) |>
    filter(sum(count) > 0) |>
    ungroup() |>
    filter(asv_id %in% names(seqs))

if (nrow(abund) == 0) stop("No ASV IDs shared between FASTA and abundance table")

abund_wide <- abund |>
    pivot_wider(names_from = sample, values_from = count, values_fill = 0)

kept <- seqs[abund_wide$asv_id]
cat("Keeping", length(kept), "of", length(seqs), "ASVs across",
    ncol(abund_wide) - 1, "samples\n")

Biostrings::writeXStringSet(kept, args$out_seqs, width = 20000)
write_tsv(abund_wide, args$out_abund)
