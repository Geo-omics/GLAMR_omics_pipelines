#! /usr/bin/env Rscript

'Plot community composition using contig abundance and LCA annotations from UniRef100 mapping.

Usage:
  plot_contig_abund_uniref_LCA_single_sample.R [options]

Options:
  -a FILE, --abund=FILE         Contig LCA abundance file.
  -s SAMPLE, --sample=SAMPLE    The sample to plot.
  -o FILE, --output=FILE        The plot file.
  -h --help                     Show this help screen.
' -> doc

library(docopt)
library(metacoder)
library(tidyverse)

# To process actual command line input
arguments <- docopt(doc)
#print(arguments)

# # for testing interactively
# arguments <- docopt(doc, args = c("--abund=data/omics/metagenomes/samp_447/samp_447_lca_abund_summarized.tsv",
#                                   "--sample=samp_447",
#                                   "--output=test_metacoder.pdf"))
# print(arguments)

pg <- DBI::dbConnect(RPostgres::Postgres(),dbname = "glamr_data", host = "cayman.earth.lsa.umich.edu", port = "5432", user = "glamr_admin", password = "glamr2023")

tax_info <- tbl(pg, "tax_info") %>% collect()

abund_from_contigs <- read_tsv(arguments$abund) %>% 
  left_join(tax_info) 

mmseqs_for_metacoder <- abund_from_contigs %>% 
  mutate(rel_abund = abund_direct / sum(abund_direct)) %>% 
  filter(!tax_name %in% c("unclassified","  cellular organisms")) %>% 
  mutate(lineage = str_glue("r__Root;{std_lineage}"),
         lineage = str_remove_all(lineage, "(;[a-z]__)*$"),
         lineage = factor(lineage, levels = unique(lineage), ordered = TRUE)) %>% 
  filter(!str_detect(lineage, ";NA")) %>% 
  select(rel_abund, lineage) %>% 
  group_by(lineage) %>% 
  summarise(rel_abund = sum(rel_abund))

obj <- parse_tax_data(mmseqs_for_metacoder, class_cols = "lineage", class_sep = ";",
                      class_key = c(tax_rank = "taxon_rank", tax_name = "taxon_name"),
                      class_regex = "^(.+)__(.*)$")

#summing per-taxon counts
obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data")

obj$data$tax_abund$S3 <- obj$data$tax_abund$rel_abund


#All domains
(tree <- obj %>%
    filter_taxa(obj$data$tax_abund$rel_abund > 0.001,supertaxa = TRUE, reassign_obs = TRUE) %>%
    heat_tree(node_label = taxon_names,
              node_size = S3,
              node_size_range = c(0.00175,0.045),
              node_label_size_range = c(.015,.025),
              node_size_axis_label = "OTU count",
              initial_layout = "reingold-tilford", layout = "davidson-harel",
              overlap_avoidance = 10,
              node_label_max = 75,
              node_color = S3,
              node_color_range = c("gray","gray","gray"),
              node_color_axis_label = "Relative abundance") +
    labs(subtitle = arguments$sample)
)
ggsave(plot = tree,filename = trimws(arguments$output), device = cairo_pdf, width = 12, height = 12, dpi = 600)

