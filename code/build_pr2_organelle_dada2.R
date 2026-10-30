#!/usr/bin/env Rscript

'Build a DADA2-format taxonomy training set of organelle 16S sequences from PR2

PR2 includes plastid and mitochondrial 16S rRNA sequences, labeled with the
taxonomy of their eukaryotic host, but only distributes its 18S as a DADA2
training set. This pulls the organelle 16S sequences out of the PR2 "merged"
metadata table and writes them with PR2 9-rank headers (Domain..Species), with
the organelle suffixes (":plas", ":mito", ...) removed from the rank names.

Usage:
  build_pr2_organelle_dada2.R --merged <PATH> --organelles <LIST> --out <PATH>

Options:
  -h --help             Show this screen
  --merged <PATH>       PR2 merged metadata table (pr2_version_*_merged.xlsx)
  --organelles <LIST>   Comma-separated values of the PR2 organelle column to keep,
                        e.g. plastid,apicoplast,chromatophore or mitochondrion
  --out <PATH>          Output DADA2 training set FASTA (.fasta.gz)
' -> doc

library(docopt)
library(readr)
library(dplyr)
library(stringr)

args <- docopt(doc)
organelles <- str_split_1(args$organelles, ",")
ranks <- c("domain", "supergroup", "division", "subdivision", "class", "order", "family", "genus", "species")

pr2 <- readxl::read_excel(args$merged, guess_max = 1e5)

refs <- pr2 |>
    filter(gene == "16S_rRNA", organelle %in% organelles, !is.na(sequence)) |>
    mutate(across(all_of(ranks), \(x) str_remove(x, ":[a-z]+$")),
           sequence = str_to_upper(str_remove_all(sequence, "[^A-Za-z]"))) |>
    filter(if_all(all_of(ranks), \(x) !is.na(x)))

if (nrow(refs) == 0) stop("No PR2 16S sequences for organelle(s): ", args$organelles)

cat("Keeping", nrow(refs), "PR2 16S sequences for organelle(s):", args$organelles, "\n")
print(count(refs, organelle))

headers <- paste0(">", do.call(paste, c(refs[ranks], sep = ";")), ";")
write_lines(c(rbind(headers, refs$sequence)), args$out)

cat("Training set written to:", args$out, "\n")
