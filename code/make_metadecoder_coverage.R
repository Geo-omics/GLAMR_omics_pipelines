#!/usr/bin/env Rscript

'Convert CoverM contig coverage profiles (tsv) to metadecoder coverage file format

Usage:
  make_metadecover_coverage.R --coverage=COVERAGE --contigs=CONTIGS --out=OUT [options]

Options:
  --coverage=COVERAGE    CoverM mean contig coverage file
  --contigs=CONTIGS      Contigs fasta file
  --out=OUT              Output coverage formated for metadecoder
  -h --help   Show this help screen.
' -> doc

library(docopt)

# To process actual command line input
arguments <- docopt(doc)
#print(arguments)

# # for testing interactively
#arguments <- docopt(doc, args = c("--coverage=data/omics/metagenomes/fa88397911af2a92053e4ca8e7499894/bins/contig_coverage.tsv --contigs=data/omics/metagenomes/fa88397911af2a92053e4ca8e7499894/assembly/megahit/final.contigs.renamed.fa --out=data/omics/metagenomes/fa88397911af2a92053e4ca8e7499894/bins/metadecoder/coverage.tsv"))
#print(arguments)

###########################
### Script starts here ####
###########################

library(tidyverse,quietly = TRUE)

#coverm_path <- "data/omics/metagenomes/fa88397911af2a92053e4ca8e7499894/bins/contig_coverage.tsv"
#contigs_path <- "data/omics/metagenomes/fa88397911af2a92053e4ca8e7499894/assembly/megahit/final.contigs.renamed.fa"
#out_path <- "data/omics/metagenomes/fa88397911af2a92053e4ca8e7499894/bins/metadecoder/coverage.tsv"

### Test files
#met_cov <- read_tsv("data/omics/metagenomes/fa88397911af2a92053e4ca8e7499894/bins/bam/COVERAGE")


make_metadecoder_cov <- function(coverm_path, contigs_path, out_path){
  contig_lengths <- Biostrings::readDNAStringSet(contigs_path) %>% 
    data.frame(Contig = names(.), seq = .) %>% 
    mutate(bin_size = nchar(seq)) %>% 
    select(Contig, bin_size)
  
  coverM <- read_tsv(coverm_path)
  
  n_covs <- length(colnames(coverM)) -1
  newcols <- paste0("coverage",1:n_covs)
  
  meta_jgi_cov <- coverM %>% 
    rename_with(~newcols,-Contig) %>% 
    mutate("bin index" = 1) %>% 
    left_join(contig_lengths) %>% 
    relocate(Contig, `bin index`,bin_size) %>% 
    rename("sequence id" = "Contig",
           "bin size" = "bin_size") %>% 
    write_tsv(out_path)
}

make_metadecoder_cov(coverm_path = arguments$coverage,
                     contigs_path = arguments$contigs, 
                     out_path = arguments$out)
