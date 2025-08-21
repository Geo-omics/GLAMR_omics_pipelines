#!/usr/bin/env Rscript

'Using dada2 with amplicon sequencing data

Usage:
  amplicon_quality_with_trunc.R [--outdir <PATH>] --quality <N> <sample_listing>

The positional argument <sample_listing> is the path to a file made by the
amplicon-dispatch script.

Options:
  -h --help             Show this screen
  --outdir=<N>          Output Directory to save the filter and trim files,
                        quality plots, seqtab.nochim, representative sequences,
                        tracking read numbers
  -q --quality=<N>      Quality threshold
' -> doc

library(dada2)
library(docopt)
library(fs)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(tibble)

arguments <- docopt(doc)

wide_file_list <- read.delim(
    arguments$sample_listing,
    sep='\t',
    header=FALSE,
    col.names=c('sample', 'sample_dir', 'swapped', 'fwd', 'rev'),
) %>% tibble() %>%
      mutate(filt_reads_fwd=str_glue('{arguments$outdir}/filtered/{sample}_fwd.fastq.gz')) %>%
      mutate(filt_reads_rev=str_glue('{arguments$outdir}/filtered/{sample}_rev.fastq.gz'))

count <- 0

for (i in seq_len(nrow(wide_file_list))) {
  sample <- wide_file_list$sample[i]
  path <- wide_file_list$sample_dir[i]

  fwdScanPath <- file.path(path, "detect_region", "fwd_summary.tsv")
  revScanPath <- file.path(path, "detect_region", "rev_summary.tsv")

  fwdScan <- read.delim(fwdScanPath, header = TRUE, sep = "\t")
  revScan <- read.delim(revScanPath, header = TRUE, sep = "\t")

  if (fwdScan$seq_start_median[1] > fwdScan$seq_end_median[1] && revScan$seq_start_median[1] < revScan$seq_end_median[1]) {
    print(paste0("Sample ", sample, " seems to have swapped forward and reverse reads"))
    count <- count + 1

    fwd_temp <- wide_file_list$fwd[i]
    rev_temp <- wide_file_list$rev[i]
    filt_fwd_temp <- wide_file_list$filt_reads_fwd[i]
    filt_rev_temp <- wide_file_list$filt_reads_rev[i]

    wide_file_list$fwd[i] <- rev_temp
    wide_file_list$rev[i] <- fwd_temp
    wide_file_list$filt_reads_fwd[i] <- filt_rev_temp
    wide_file_list$filt_reads_rev[i] <- filt_fwd_temp
  }
}

print(paste0(count, "/", nrow(wide_file_list), " of the samples are swapped"))

forwardReads <- wide_file_list$fwd
reverseReads <- wide_file_list$rev

# extract sample names
namesOfSamples <- wide_file_list$sample

fs::dir_create(path = arguments$outdir)

# visualize and save quality of forward and reverse reads
(forwardPlot <- plotQualityProfile(forwardReads[1:2]))
plotDataFwd <- forwardPlot$data
ggsave(filename = str_glue("{arguments$outdir}/forward_quality_plot.pdf"),
       plot = forwardPlot, width = 5, height = 3, scale = 2)

(reversePlot <- plotQualityProfile(reverseReads[1:2]))
plotDataRev <- reversePlot$data
ggsave(filename = str_glue("{arguments$outdir}/reverse_quality_plot.pdf"),
       plot = reversePlot, width = 5, height = 3, scale = 2)

base_dir <- dirname(dirname(forwardReads[3]))
fwdScanSummary <- file.path(base_dir, "detect_region", "fwd_summary.tsv")
revScanSummary <- file.path(base_dir, "detect_region", "rev_summary.tsv")

fwd_summary_t <- read.delim(fwdScanSummary, header = TRUE, sep = "\t")
rev_summary_t <- read.delim(revScanSummary, header = TRUE, sep = "\t")

if (fwd_summary_t$seq_start_median[1] > fwd_summary_t$seq_end_median[1]
    && rev_summary_t$seq_start_median[1] < rev_summary_t$seq_end_median[1]) {
      fwd_summary <- rev_summary_t
      rev_summary <- fwd_summary_t
} else {
  fwd_summary <- fwd_summary_t
  rev_summary <- rev_summary_t
}

# using weighted mean
resultMeanFwd <- plotDataFwd %>%
  group_by(Cycle) %>%
  mutate(mean = sum(Score * Count) / sum(Count)) %>%
  slice_max(Count, with_ties = FALSE)

resultMeanRev <- plotDataRev %>%
  group_by(Cycle) %>%
  mutate(mean = sum(Score * Count) / sum(Count)) %>%
  slice_max(Count, with_ties = FALSE)

threshold <- arguments$quality

# fwd_last is meant to be the location on the hmm model that the forward read ends reading
# fwd_first is meant to be the location on the hmm model that the forward read starts reading
# rev_last is meant to be the location on the hmm model that the reverse read ends reading
# rev_first is meant to be the location on the hmm model that the reverse read starts reading
# so rev_last < rev_first and fwd_last > fwd_first
#       rev_last --V          V-- rev_first
#                  ------------
#             -------------
# fwd_first --^           ^-- fwd_last
fwd_last <- fwd_summary$hmm_end_median # hmm_end_median
fwd_first <- fwd_summary$hmm_start_median # hmm_start_median
rev_last <- rev_summary$hmm_start_median # hmm_start_median
rev_first <- rev_summary$hmm_end_median # hmm_end_median


# Note: length of start and end on hmm does not necessarily correlate with length of start and end on seq

# means the sequence is at follows
# rev_first < fwd_first scenario:
#       ---------                   <-- reverse read
#                  ----------       <-- forward read

# fwd_last < rev_last scenario:
#               ----------          <-- reverse read
#    --------                       <-- forward read
if (rev_first < fwd_first || fwd_last < rev_last) {
  # print out: No overlap.
  # stop the program
  stop("No overlap")
}

if (rev_last < fwd_first && rev_first < fwd_last) { # Case 1
  # means sequence is as follows
  #      -------------              <-- reverse read
  #          -------------          <-- forward read
  if (abs(rev_first - fwd_first) < 20) {
    stop(paste("Not enough overlap. Expected overlap is 20+. Overlap is", abs(rev_first - fwd_first)))
  }

  begin_overlap_hmm <- fwd_first
  end_overlap_hmm <- rev_first

  # left to right
  begin_overlap_seq_fwd <- fwd_summary$seq_start_median
  end_overlap_seq_fwd <- rev_first - fwd_first + 1

  # left to right
  begin_overlap_seq_rev <- rev_summary$seq_start_median - (fwd_first - rev_last)
  end_overlap_seq_rev <- rev_summary$seq_end_median

} else if (rev_last < fwd_first && fwd_last < rev_first) { # Case 2
  # below assuming two sequences
  #     ---------------               <-- reverse read
  #        -------                    <-- forward read

  if (fwd_last - fwd_first < 20) {
    stop(paste("Not enough overlap. Expected overlap is 20+. Overlap is", abs(fwd_last - fwd_first)))
  }

  begin_overlap_hmm <- fwd_first
  end_overlap_hmm <- fwd_last

  # left to right
  begin_overlap_seq_fwd <- fwd_summary$seq_start_median
  end_overlap_seq_fwd <- fwd_summary$seq_end_median

  # left to right
  begin_overlap_seq_rev <- rev_summary$seq_start_median - (fwd_first - rev_last)
  end_overlap_seq_rev <- begin_overlap_seq_rev - (end_overlap_seq_fwd - begin_overlap_seq_fwd)

} else if (fwd_first < rev_last && rev_first < fwd_last) { # Case 3
  # below assuming two sequences
  #         --------                  <-- reverse read
  #     --------------                <-- forward read
  if (rev_first - rev_last < 20) {
    stop(paste("Note enough overlap. Expected overlap is 20+. Overlap is", rev_first - rev_last))
  }

  begin_overlap_hmm <- rev_last
  end_overlap_hmm <- rev_first

  # left to right
  begin_overlap_seq_fwd <- rev_last - fwd_first + 1
  end_overlap_seq_fwd <- rev_first - fwd_first + 1

  # left to right
  begin_overlap_seq_rev <- rev_summary$seq_start_median
  end_overlap_seq_rev <- rev_summary$seq_end_median

} else { # Case 4
  # below assuming normal two sequences
  #           ---------------         <-- reverse read
  #     ----------------              <-- forward read
  if (abs(fwd_last - rev_last) < 20) {
    stop(paste("Not enough overlap. Expected overlap is 20+. Overlap is", abs(fwd_last - rev_last)))
  }

  begin_overlap_hmm <- rev_last
  end_overlap_hmm <- fwd_last

  # left to right
  begin_overlap_seq_fwd <- rev_last - fwd_first + 1
  end_overlap_seq_fwd <- fwd_summary$seq_end_median

  # left to right
  begin_overlap_seq_rev <- rev_summary$seq_start_median
  end_overlap_seq_rev <- rev_summary$seq_start_median - (end_overlap_seq_fwd - begin_overlap_seq_fwd)
}

sum <- 0

# obtain first rolling avg
for (i in 0:19) {
  sum <- sum + resultMeanFwd$mean[i + begin_overlap_seq_fwd] + resultMeanRev$mean[begin_overlap_seq_rev - i]
}

revTrunc <- begin_overlap_seq_rev
forTrunc <- begin_overlap_seq_fwd + 20
optimal_sum <- sum

for (i in 1:(end_overlap_seq_fwd - begin_overlap_seq_fwd - 20 + 1)) {
  sum <- sum - resultMeanFwd$mean[i + begin_overlap_seq_fwd - 1] - resultMeanRev$mean[begin_overlap_seq_rev - i + 1]
  sum <- sum + resultMeanFwd$mean[i + begin_overlap_seq_fwd + 20 - 1] + resultMeanRev$mean[begin_overlap_seq_rev - i - 20 + 1]

  if (optimal_sum < sum) {
    optimal_sum <- sum
    revTrunc <- begin_overlap_seq_rev - i
    forTrunc <- begin_overlap_seq_fwd + 20 + i - 1
  }
}

continue <- TRUE
# extend to the "right"
while (forTrunc + 1 <= end_overlap_seq_fwd && continue) {
  continue <- FALSE
  sum <- resultMeanFwd$mean[forTrunc] + resultMeanFwd$mean[forTrunc - 1] + resultMeanFwd$mean[forTrunc + 1]
  if (sum / 3 > threshold) {
    continue <- TRUE
    forTrunc <- forTrunc + 1
  }
}

continue <- TRUE
# extend to the "left
while (revTrunc + 1 <= begin_overlap_seq_rev & continue) {
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
  write_tsv(str_glue("{arguments$outdir}/filt_and_trim.tsv"))

# learning errors
errorForward <- learnErrors(filtAndTrimForward, multithread=TRUE)
errorReverse <- learnErrors(filtAndTrimReverse, multithread=TRUE)

# saving errors to an rds file
errors <- list(errorForward, errorReverse)
write_rds(errors, str_glue("{arguments$outdir}/errors.rds"))

# saving the forward error plot to a file
forwardErrorPlot <- plotErrors(errorForward, nominalQ=TRUE)
ggsave(filename = str_glue("{arguments$outdir}/forward_error_plot.pdf"),
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
  `colnames<-`(openssl::md5(colnames(seqtable_nochimeras))) |>
  write_tsv(str_glue("{arguments$outdir}/asv_table.tsv"))

# output representative sequences to fasta
headers <- openssl::md5(colnames(seqtable_nochimeras))
Biostrings::DNAStringSet(colnames(seqtable_nochimeras)) |>
  `names<-`(headers) |>
  Biostrings::writeXStringSet(str_glue("{arguments$outdir}/rep_seqs.fasta"))

# track reads through pipelines
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaForward, getN), sapply(dadaReverse, getN), sapply(mergeReads, getN), rowSums(seqtable_nochimeras))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- namesOfSamples

# track reads saved to .tsv file
track |>
  as.data.frame() |>
  rownames_to_column("sample") |>
  write_tsv(str_glue("{arguments$outdir}/track_read_counts.tsv"))
