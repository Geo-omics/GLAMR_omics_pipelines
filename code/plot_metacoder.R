#! /usr/bin/env Rscript

'Plot community composition using metacoder.

Usage:
  plot_metacoder.R [options]

Options:
  -a FILE, --abund=FILE       Combined relative abundance table
  -s SAMPLE, --sample=SAMPLE  The sample to plot.
  -o FILE, --output=FILE      The plot file.
  -h --help                   Show this help screen.
' -> doc

library(docopt)

# To process actual command line input
arguments <- docopt(doc)
#print(arguments)

# # for testing interactively
# arguments <- docopt(doc, args = c("--abund=data/sample_data/bracken_rel_abund.tsv",
#                                   "--sample=0245e67227092c6391e231e929be7edd",
#                                   "--output=test_metacoder.pdf"))
# print(arguments)


library(metacoder)
library(tidyverse)

rel.abund <- read_tsv(trimws(arguments$abund)) 

sample_abund <- rel.abund %>%
  dplyr::select(arguments$sample,taxonomy) %>%
  filter(trimws(arguments$sample) > 0)


obj <- parse_tax_data(sample_abund, class_cols = "taxonomy", class_sep = ";",
                      class_key = c(tax_rank = "taxon_rank", tax_name = "taxon_name"),
                      class_regex = "^(.+)__(.*)$")

#converting to relative abundance
obj$data$tax_data <- calc_obs_props(obj, "tax_data")

#summing per-taxon counts
obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data")

obj$data$tax_abund <- obj$data$tax_abund %>% rename(rel_abund = trimws(arguments$sample))

#obj2 <- obj %>% filter_taxa(taxon_ranks == "g", supertaxa = TRUE, reassign_obs = TRUE, drop_obs = TRUE)

(tree <- obj %>%
    metacoder::filter_taxa(taxon_ranks == "g", supertaxa = TRUE,reassign_obs = TRUE, subtaxa = FALSE) %>% # subset to the genus rank
    metacoder::filter_obs("tax_abund", rel_abund > 0.005,drop_taxa = TRUE,supertaxa = TRUE,reassign_obs = TRUE) %>% #filter to only plot taxa that represent at least 0.5% of the community
    heat_tree(node_label = taxon_names,
              node_size = rel_abund,
              node_size_range = c(0.00175,0.045),
              node_label_size_range = c(.015,.025),
              node_size_axis_label = "OTU count",
              initial_layout = "reingold-tilford", layout = "davidson-harel",
              overlap_avoidance = 10,
              node_label_max = 75,
              node_color = rel_abund,
              node_color_range = c("gray","gray","gray"),
              node_color_axis_label = "Relative abundance")
)
ggsave(plot = tree,filename = trimws(arguments$output), device = cairo_pdf(), width = 12, height = 12, dpi = 600)