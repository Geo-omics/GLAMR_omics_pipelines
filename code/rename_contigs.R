#!/usr/bin/env Rscript

'Rename MegaHit contigs & export information from header as a .tsv file

Usage:
  rename_megahit_contigs.R -i contigs  [options]

Options:
  -i CONTIGS  Uncompressed contigs .fa produced by MegaHit
  -o RENAMED  --out=RENAMED     Renamed assembly
  -s INFO     --supp=INFO       Supplemental information reported by the assembler (.tsv file)
  -p PREFIX   --prefix=PREFIX   Prefix to use for contig names (e.g. sample name)
  -h --help                     Show this help screen.
' -> doc

library(docopt)

# To process actual command line input
arguments <- docopt(doc)
#print(arguments)

# # for testing interactively
#arguments <- docopt(doc, args = c(" -i ~/GLAMR/data/omics/metagenomes/0308a1e0034cfa73ce078604e79e9ee7/assembly/megahit/final.contigs.fa -o ~/GLAMR/data/omics/metagenomes/0308a1e0034cfa73ce078604e79e9ee7/assembly/megahit/final.contigs.renamed.fa -s supplemental_info.tsv -p samp1"))
#print(arguments)

###########################
### Script starts here ####
###########################

#contigs <- "data/projects/2022_geomicro_JGI_CSP/metagenomes/samp_447/assembly/metaspades_noNORM/contigs.fasta"

library(tidyverse,quietly = TRUE)

if (!is.null(arguments$`--contigs`)) {
  contig_path <- arguments$`--contigs`
} else if(!is.null(arguments$contigs)){
  contig_path <- arguments$contigs
} else if(!is.null(arguments$i)){
  contig_path <- arguments$i
}

rename_contigs <- function(contigs){
  
  sample <- str_remove(contigs, ".*metagenomes/") %>% str_remove("/assembly.*")
  contig_info_fp <- paste0(dirname(contigs),"/contigs_info.tsv")
  
  contigs_df <- Biostrings::readDNAStringSet(contigs) %>% 
    data.frame(header = names(.),
               seq = .)
  getwd()
  
  if(contigs_df$header[1] == paste0(sample,"_","1")){
    print("Contigs have already been processed")
    quit(save = "no", status = 1)
  }
  
  # Check if megahit assembly
  if(str_detect(contigs_df$header[1], "^k")){
    contigs_df <- contigs_df %>% 
      separate(header, into = c("orig_id","flag","approx_cov", "length"),sep = " ") %>% 
      mutate(across(c(flag,approx_cov, length), ~str_remove(.,".*="))) %>% 
      type_convert() %>%
      mutate(sample = sample,
             contig_id = paste0(sample,"_",row_number())) %>% 
      relocate(contig_id)
    
    contig_info <- contigs_df %>% 
      select(contig_id, orig_id, flag, approx_cov, length) %>% 
      write_tsv(contig_info_fp)
    
    out_path <- paste0(str_remove(contigs,"\\.fa"),".renamed.fa")
    
  } else if(str_detect(contigs_df$header[1], "^NODE")){   # Check if metaSpades assembly
    contigs_df <- contigs_df %>% 
      separate(header, into = c("remove1","contig_num","remove2", "length", "remove3", "approx_cov"),sep = "_") %>% 
      select(-starts_with("remove")) %>%  
      type_convert() %>%
      mutate(sample = sample,
             contig_id = paste0(sample,"_",row_number())) %>% 
      relocate(contig_id)
    
    contig_info <- contigs_df %>% 
      select(contig_id, contig_num, approx_cov, length) %>% 
      write_tsv(contig_info_fp)
  
    out_path <- paste0(str_remove(contigs,"\\.fasta"),".renamed.fasta")
    }
  
  for_export <- Biostrings::DNAStringSet(contigs_df$seq)
  names(for_export) <- contigs_df$contig_id
  
  Biostrings::writeXStringSet(for_export,out_path)
}

rename_contigs(contig_path)