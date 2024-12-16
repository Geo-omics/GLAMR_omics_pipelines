#!/usr/bin/env Rscript

'Using dada2 with amplicon sequencing data

Usage:
  amplicon_quality_with_trunc.R [--input FILE] [--directory DIRECTORY] --quality <N>

Options:
  -h --help             Show this screen
  -d --directory=<N>    Output Directory to save the filter and trim files, quality plots, seqtab.nochim, representative sequences, tracking read numbers
  -q --quality=<N>      Quality threshold
  -i --input=<N>        Fastq file
' -> doc

library(dada2)
library(docopt)
library(fs)
library(dplyr)
library(readr)

# the below arguments for processing from the command line
#arguments <- docopt(doc)

# uncomment the below arguments for testing
arguments <- docopt(doc, args = "--input ~/amplicon_sequencing/amplicons --directory ~/amplicon_sequencing/dataTest2 --quality 20")

inPath <- arguments$input

# samp_1/detect_region/fwd.txt

file_list <- fs::dir_ls(path = inPath, recurse = TRUE) %>%
  data_frame(path = .) %>%
  filter(str_detect(path, ".fastq.gz")) %>%
  unglue::unglue_unnest(path, "{reads_dir}/samp_{sample_num}/reads/raw_{dir}_reads.fastq.gz", remove = FALSE) %>%
  filter(!str_detect(sample_num, "NA"),
         !is.na(sample_num))

wide_file_list <- file_list %>%
  pivot_wider(id_cols = c("reads_dir", "sample_num"), names_from = dir, values_from = path) %>%
  mutate(filt_reads_fwd = str_glue("{reads_dir}/filt_reads/{sample_num}_F_filt.fastq.gz"),
         filt_reads_rev = str_glue("{reads_dir}/filt_reads/{sample_num}_R_filt.fastq.gz"))

forwardReads <- wide_file_list$fwd

reverseReads <- wide_file_list$rev

# extract sample names
namesOfSamples <- wide_file_list$sample_num

fs::dir_create(path = arguments$directory)

# visualize and save quality of forward and reverse reads
(forwardPlot <- plotQualityProfile(forwardReads[1:2]))
plotDataFwd <- forwardPlot$data
ggsave(filename = str_glue("{arguments$directory}/forward_quality_plot.pdf"),
       plot = forwardPlot, width = 5, height = 3, scale = 2)

(reversePlot <- plotQualityProfile(reverseReads[1:2]))
plotDataRev <- reversePlot$data
ggsave(filename = str_glue("{arguments$directory}/reverse_quality_plot.pdf"),
       plot = reversePlot, width = 5, height = 3, scale = 2)

base_dir <- dirname(dirname(forwardReads[3]))
fwdScanSummary <- file.path(base_dir, "detect_region", "fwd_summary.tsv")
revScanSummary <- file.path(base_dir, "detect_region", "rev_summary.tsv")

fwd_summary <- read.delim(fwdScanSummary, header = TRUE, sep = "\t")
rev_summary <- read.delim(revScanSummary, header = TRUE, sep = "\t")

resultMeanFwd <- plotDataFwd %>%
  group_by(Cycle) %>%
  mutate(mean = sum(Score * Count) / sum(Count)) %>%
  slice_max(Count, with_ties = FALSE)

resultMeanRev <- plotDataRev %>%
  group_by(Cycle) %>%
  mutate(mean = sum(Score * Count) / sum(Count)) %>%
  slice_max(Count, with_ties = FALSE)

threshold <- arguments$quality

fwd_last <- fwd_summary$hmm_end_median # hmm_end_median
fwd_first <- fwd_summary$hmm_start_median # hmm_start_median
rev_last <- rev_summary$hmm_start_median # hmm_start_median
rev_first <- rev_summary$hmm_end_median # hmm_end_median

if (rev_first < fwd_first || fwd_last < rev_last) {
  # print out: No overlap.
  # stop the program
  stop("No overlap")
}

if (rev_last < fwd_first) {
  # means sequence is as follows. Top line is reverse. Bottom line is forward
  #      -------------
  #          -------------
  if (abs(rev_first - fwd_first) < 20) {
    stop(paste("Not enough overlap. Expected overlap is 20+. Overlap is", abs(fwd_last - rev_last)))
  }

  begin_overlap_hmm <- fwd_first
  end_overlap_hmm <- rev_first

  begin_overlap_seq <- begin_overlap_hmm - fwd_first + 1
  end_overlap_seq <- end_overlap_hmm - fwd_first + 1

  forStart <- begin_overlap_seq
  revStart <- end_overlap_seq

} else {
  # below assuming normal two sequences. Top line is reverse. Bottom line is forward
  #           ---------------
  #     ----------------
  if (abs(fwd_last - rev_last) < 20) {
    stop(paste("Not enough overlap. Expected overlap is 20+. Overlap is", abs(fwd_last - rev_last)))
  }

  begin_overlap_hmm <- rev_last
  end_overlap_hmm <- fwd_last

  begin_overlap_seq <- begin_overlap_hmm - fwd_first + 1
  end_overlap_seq <- end_overlap_hmm - fwd_first + 1

  forStart <- begin_overlap_seq
  revStart <- nrow(resultMeanRev)
}

sum <- 0

# obtain first rolling avg
for (i in 0:19) {
  sum <- sum + resultMeanFwd$mean[i + forStart] + resultMeanRev$mean[revStart - i]
}

revTrunc <- revStart
forTrunc <- forStart + 20
optimal_sum <- sum

for (i in 1:(end_overlap_seq - begin_overlap_seq - 20 + 1)) {
  sum <- sum - resultMeanFwd$mean[i + forStart - 1] - resultMeanRev$mean[revStart - i + 1]
  sum <- sum + resultMeanFwd$mean[i + forStart + 20 - 1] + resultMeanRev$mean[revStart - i - 20 + 1]

  if (optimal_sum < sum) {
    optimal_sum <- sum
    revTrunc <- revStart - i
    forTrunc <- forStart + 20 + i - 1
  }
}

continue <- TRUE
# extend to the "right"
while (forTrunc + 1 <= end_overlap_seq & continue) {
  continue <- FALSE
  sum <- resultMeanFwd$mean[forTrunc] + resultMeanFwd$mean[forTrunc - 1] + resultMeanFwd$mean[forTrunc + 1]
  if (sum / 3 > threshold) {
    continue <- TRUE
    forTrunc <- forTrunc + 1
  }
}

continue <- TRUE
# extend to the "left
while (revTrunc + 1 <= revStart & continue) {
  continue <- FALSE
  sum <- resultMeanRev$mean[revTrunc] + resultMeanRev$mean[revTrunc - 1] + resultMeanRev$mean[revTrunc + 1]
  if (sum / 3 > threshold) {
    continue <- TRUE
    revTrunc <- revTrunc + 1
  }
}

paste("Using forward truncation", forTrunc, "and reverse truncation", revTrunc)

# filter and trim
filtAndTrimForward <- wide_file_list$filt_reads_fwd
filtAndTrimReverse <- wide_file_list$filt_reads_rev

names(filtAndTrimForward) <- namesOfSamples
names(filtAndTrimReverse) <- namesOfSamples

out <- filterAndTrim(forwardReads, filtAndTrimForward, reverseReads, filtAndTrimReverse,
                     truncLen=c(forTrunc, revTrunc),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # on macos multithread = true

# saving the filter and trim to .tsv file
out |>
  as.data.frame() |>
  rownames_to_column("file") |>
  write_tsv(str_glue("{arguments$directory}/filt_and_trim.tsv"))

# learning errors
errorForward <- learnErrors(filtAndTrimForward, multithread=TRUE)
errorReverse <- learnErrors(filtAndTrimReverse, multithread=TRUE)

# saving errors to an rds file
errors <- list(errorForward, errorReverse)
write_rds(errors, str_glue("{arguments$directory}/errors.rds"))

# saving the forward error plot to a file
forwardErrorPlot <- plotErrors(errorForward, nominalQ=TRUE)
ggsave(filename = str_glue("{arguments$directory}/forward_error_plot.pdf"),
       plot = forwardErrorPlot, width = 5, height = 3, scale = 2)

# sample inference
dadaForward <- dada(filtAndTrimForward, err=errorForward, multithread=TRUE)
dadaReverse <- dada(filtAndTrimReverse, err=errorReverse, multithread=TRUE)

# merging forward and reverse reads
mergeReads <- mergePairs(dadaForward, filtAndTrimForward, dadaReverse, filtAndTrimReverse, verbose=TRUE)

# construct sequence table
seqtable <- makeSequenceTable(mergeReads)

# remove chimeras, saved to .tsv file
seqtable_nochimeras <- removeBimeraDenovo(seqtable, method="consensus", multithread=TRUE, verbose=TRUE)
seqtable_nochimeras |>
  as.data.frame() |>
  `colnames<-`(openssl::md5(colnames(seqtab.nochim))) |>
  write_tsv(str_glue("{arguments$directory}/asv_table.tsv"))

# output representative sequences to fasta
headers <- openssl::md5(colnames(seqtable_nochimeras))
Biostrings::DNAStringSet(colnames(seqtable_nochimeras)) |>
  `names<-`(headers) |>
  Biostrings::writeXStringSet(str_glue("{arguments$directory}/rep_seqs.fasta"))

# track reads through pipelines
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaForward, getN), sapply(dadaReverse, getN), sapply(mergeReads, getN), rowSums(seqtable_nochimeras))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- namesOfSamples

# track reads saved to .tsv file
track |>
  as.data.frame() |>
  rownames_to_column("sample") |>
  write_tsv(str_glue("{arguments$directory}/track_read_counts.tsv"))
