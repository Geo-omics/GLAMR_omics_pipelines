#!/usr/bin/env Rscript

'Assign taxonomy to ASVs from vsearch global alignment top hits

Each ASV gets the consensus taxonomy of its top hits (all hits tied at the
best identity, per part for concatenated ASVs): a rank is assigned if at least --min-consensus of the top hits
agree, tolerating some mislabeled references.  A rank is also only assigned if
the best identity reaches that rank\'s threshold, e.g. for COI genus needs 95%
and species 98%.  Output has the same format as dada2_assign_taxonomy.R
(sequences in the first column).

Usage:
  vsearch_lca_taxonomy.R --seqs <PATH> --hits <PATH> --out <PATH> --tax-levels <LIST> --min-identity <LIST> [--min-consensus <FRAC>]

Options:
  -h --help              Show this screen
  --seqs <PATH>          ASV FASTA, the one given to vsearch
  --hits <PATH>          vsearch --userout table with fields query+id+target,
                         target labels are DADA2-style "rank1;rank2;...;",
                         queries of concatenated ASV parts are <ASV>|<part>
  --out <PATH>           Output TSV file path
  --tax-levels <LIST>    Comma-separated rank names, one per rank in the reference
  --min-identity <LIST>  Comma-separated minimum identity (0-1) per rank, same
                         order as --tax-levels
  --min-consensus <FRAC> Minimum fraction of top hits agreeing on a taxon [default: 0.8]
' -> doc

suppressPackageStartupMessages({
    library(docopt)
    library(tidyverse)
})

args <- docopt(doc)
tax_levels <- str_split_1(args$tax_levels, ",")
min_consensus <- as.numeric(args$min_consensus)
min_identity <- set_names(as.numeric(str_split_1(args$min_identity, ",")), tax_levels)
stopifnot(length(min_identity) == length(tax_levels), !anyNA(min_identity))

seqs <- Biostrings::readDNAStringSet(args$seqs)
asvs <- tibble(
    query = str_remove(names(seqs), " .*"),  # vsearch truncates labels at the first space
    asv_id = as.character(seqs),
)
# Concatenated (non-overlapping) ASVs are searched as separate parts, split at
# the N spacer, see rule amplicon_vsearch_hits, with queries named <ASV>|<part>
parts <- asvs %>%
    mutate(part = str_split(asv_id, "N{10}")) %>%
    unnest_longer(part, indices_to = "part_num") %>%
    transmute(query, part_num, part_len = nchar(part))

hits <- read_tsv(args$hits, col_names = c("part_query", "identity", "target"), col_types = "cdc") %>%
    separate_wider_delim(part_query, "|", names = c("query", "part_num"), too_few = "align_start") %>%
    mutate(identity = identity / 100, part_num = as.integer(replace_na(part_num, "1"))) %>%
    separate_wider_delim(target, ";", names = tax_levels, too_many = "drop", too_few = "align_start") %>%
    group_by(query, part_num) %>%
    filter(identity == max(identity)) %>%
    ungroup()
cat("Top hits for", n_distinct(hits$query), "of", nrow(asvs), "ASVs\n")

# identity over all parts, weighted by part length, parts without hits count as 0
identity <- parts %>%
    left_join(distinct(hits, query, part_num, identity), by = c("query", "part_num")) %>%
    group_by(query) %>%
    summarize(identity = sum(replace_na(identity, 0) * part_len) / sum(part_len))

# fraction of top hits per taxon, averaged over parts (all parts need to agree)
consensus <- hits %>%
    pivot_longer(all_of(tax_levels), names_to = "rank", values_to = "taxon") %>%
    mutate(taxon = na_if(taxon, "")) %>%
    count(query, part_num, rank, taxon) %>%
    group_by(query, part_num, rank) %>%
    mutate(frac = n / sum(n)) %>%
    left_join(count(parts, query, name = "num_parts"), by = "query") %>%
    group_by(query, rank, taxon) %>%
    summarize(frac = sum(frac) / first(num_parts), .groups = "drop_last") %>%
    slice_max(frac, n = 1, with_ties = FALSE) %>%
    ungroup()

lca <- consensus %>%
    left_join(identity, by = "query") %>%
    mutate(taxon = if_else(frac >= min_consensus & identity >= min_identity[rank], taxon, NA_character_)) %>%
    select(query, rank, taxon) %>%
    pivot_wider(names_from = rank, values_from = taxon)

# a lower rank can have a consensus below an unset higher rank, cut those
taxa <- asvs %>%
    left_join(lca, by = "query") %>%
    select(asv_id, all_of(tax_levels))
for (i in seq_along(tax_levels)[-1]) {
    taxa[[tax_levels[i]]][is.na(taxa[[tax_levels[i - 1]]])] <- NA_character_
}

for (rank in tax_levels) {
    cat(sprintf("  %-10s assigned: %d ASVs\n", rank, sum(!is.na(taxa[[rank]]))))
}
write_tsv(taxa, args$out)
cat("\nTaxonomy written to:", args$out, "\n")
