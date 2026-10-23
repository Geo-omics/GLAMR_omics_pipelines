#! /usr/bin/env Rscript

'Plot community composition using metacoder.

Usage:
  plot_metacoder.R [options]

Options:
  -r FILE, --abund-refseq=FILE  Bracken Refseq file.
  -g FILE, --abund-gtdb=FILE    Bracken GTDB file.
  -t FILE, --tax-ref=FILE       The reference taxonomy file.
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
# arguments <- docopt(doc, args = c("--abund-refseq=data/omics/metagenomes/samp_447/kraken_fastp/refseq_E20150029_bracken.txt",
#                                   "--abund-gtdb=data/omics/metagenomes/samp_447/kraken_fastp/gtdb_E20150029_bracken.txt",
#                                   "--tax-ref=data/reference/kraken_tax_info_merged.tsv",
#                                   "--sample=samp_447",
#                                   "--output=test_metacoder.pdf"))
# print(arguments)

tax_levels <- c("Domain","Phylum","Class","Order","Family","Genus","Species")

merged_tax_summary <- read_tsv(trimws(arguments$tax_ref))

#list.files(here::here(),pattern = "_bracken.txt",full.names = TRUE,recursive = TRUE) %>% 

abund_mat <- bind_rows(read_tsv(arguments$abund_gtdb) %>% 
                         mutate(database = "gtdb",
                                taxonomy_id = str_glue("{database}_{taxonomy_id}")),
                       read_tsv(arguments$abund_refseq) %>% 
                         mutate(database = "refseq",
                                taxonomy_id = str_glue("{database}_{taxonomy_id}"))) %>% 
  left_join(merged_tax_summary) %>% 
  filter(!(database == "refseq" & str_detect(taxonomy,"d__Archaea")),
         !(database == "refseq" & str_detect(taxonomy,"d__Bacteria"))) %>% 
  mutate(rel_abund = new_est_reads / sum(new_est_reads),
         sample = arguments$sample)


wide_counts <- abund_mat %>% 
  ungroup() %>% 
  select(sample, taxonomy_id, taxonomy, database, all_of(tax_levels),new_est_reads) %>% 
  pivot_wider(names_from = sample, values_from = new_est_reads,values_fill = 0)

wide_rel_abund <- abund_mat %>%
  ungroup() %>% 
  select(sample, taxonomy_id, taxonomy, database, all_of(tax_levels),rel_abund) %>% 
  pivot_wider(names_from = sample, values_from = rel_abund,values_fill = 0) 

sample_abund <- wide_rel_abund %>%
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