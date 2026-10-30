#!/usr/bin/env Rscript

'Flag COI ASVs that are likely pseudogenes (NUMTs) or errors

COI is protein coding, so a real ASV should translate without stop codons, and
its length should differ from the typical ASV length by a multiple of 3.  All
ASVs of a target start at the same primer, so they share a reading frame: the
frame in which most ASVs translate without stops is used for all of them.  ASVs
are translated with the invertebrate (5) and vertebrate (2) mitochondrial codes
and the lower stop codon count is reported.  ASVs are only flagged, the ASV
table is not changed.

Usage:
  coi_asv_check.R --seqs <PATH> --out <PATH>

Options:
  -h --help       Show this screen
  --seqs <PATH>   ASV FASTA (asvs.fasta from the amplicon_asv_check rule)
  --out <PATH>    Output TSV, one row per ASV
' -> doc

library(docopt)
library(readr)
library(dplyr)
library(purrr)
library(stringr)

args <- docopt(doc)

GENETIC_CODES <- c(invertebrate_mito = "5", vertebrate_mito = "2")

seqs <- Biostrings::readDNAStringSet(args$seqs)
names(seqs) <- str_remove(names(seqs), " .*")

# stop codons per forward frame, fewest over both genetic codes
count_stops <- function(seq) {
    expand.grid(frame = 0:2, code = names(GENETIC_CODES), stringsAsFactors = FALSE) |>
        mutate(stops = map2_int(frame, code, function(frame, code) {
            s = Biostrings::subseq(seq, start = frame + 1)
            s = Biostrings::subseq(s, end = length(s) - length(s) %% 3)
            aa = Biostrings::translate(s, genetic.code = Biostrings::getGeneticCode(GENETIC_CODES[[code]]),
                                       if.fuzzy.codon = "X")
            str_count(as.character(aa), fixed("*"))
        })) |>
        slice_min(stops, n = 1, with_ties = FALSE, by = frame)
}

lengths <- setNames(Biostrings::width(seqs), names(seqs))
modal_length <- as.integer(names(which.max(table(lengths))))

stops <- map(as.list(seqs), count_stops) |>
    list_rbind(names_to = "asv_id")

target_frame <- stops |>
    summarise(n_stop_free = sum(stops == 0), .by = frame) |>
    slice_max(n_stop_free, n = 1, with_ties = FALSE) |>
    pull(frame)

checks <- stops |>
    filter(frame == target_frame) |>
    mutate(length = lengths[asv_id],
           length_offset = length - modal_length,
           in_frame = length_offset %% 3 == 0,
           likely_pseudogene = stops > 0 | !in_frame) |>
    select(asv_id, length, length_offset, in_frame, min_stop_codons = stops,
           frame, best_code = code, likely_pseudogene)

write_tsv(checks, args$out)

cat("Modal ASV length:", modal_length, " reading frame:", target_frame, "\n")
cat(sum(checks$likely_pseudogene), "of", nrow(checks), "ASVs flagged as likely pseudogenes/errors\n")
cat("Written to:", args$out, "\n")
