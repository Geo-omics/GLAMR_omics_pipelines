#!/usr/bin/env Rscript

'Using dada2 with amplicon sequencing data

Usage:
  amplicon_quality_with_trunc.R [--glamr-root <PATH>] --quality <N> --outdir <PATH> [--cpus <N>] --assignments <PATH> --samples <PATH> --targets <PATH> <target>

The positional argument <sample_listing> is the path to a file made by the
amplicon-dispatch script.

Options:
  -h --help             Show this screen
  --outdir <PATH>       Output Directory to save the filter and trim files,
                        quality plots, seqtab.nochim, representative sequences,
                        tracking read numbers
  --assignments <PATH>  Path to amplicon target assignment listing as produced
                        by the amplicon-dispatch script.
  --samples <PATH>      Path to the samples and files listing as produced by
                        the amplicon-dispatch script.
  --targets <PATH>      Path to the target_info table.
                        the amplicon-dispatch script.
  --glamr-root <PATH>   GLAMR root directory, certain paths read from the other
                        input files will be interpreted as relative to this root directory
                        [default: ./]
  -q --quality=<N>      Quality threshold
  --cpus=<N>            Number of threads/cpus [default: 1]
' -> doc

library(Rcpp)
library(dada2)
library(docopt)
library(fs)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(tibble)

MIN_OVERLAP = 20


args <- docopt(doc)
cpus = as.integer(args$cpus)


assignments <- read.delim(
    args$assignments,
    header=TRUE,
    fill=TRUE,
) %>% tibble() %>%
    mutate(assigned=if_else(override == "", target, override, missing=target)) %>%
    filter(assigned == args$target)

if (nrow(assignments) == 0) {
    stop('no samples where assigned to given target')
}

target_tab <- read.delim(args$targets, header=TRUE)

samples <- read.delim(
    args$samples,
    header=TRUE,
) %>% tibble() %>%
    semi_join(assignments, by=join_by('sample')) %>%
    left_join(target_tab, by=join_by('sample')) %>%
    mutate(filt_reads_fwd=str_glue('{args$outdir}/filtered/{sample}_fwd.fastq.gz')) %>%
    mutate(filt_reads_rev=str_glue('{args$outdir}/filtered/{sample}_rev.fastq.gz')) %>%
    # below: paste glamr_root and paths, but respect absolute paths in data,
    # make some effort to avoid double slashes
    mutate(sample_dir=if_else(str_starts(sample_dir, '/'),
                              sample_dir,
                              file.path(str_remove(args$glamr_root, '/$'), sample_dir))
    ) %>%
    mutate(fwd_fastq=if_else(str_starts(sample_dir, '/'),
                              fwd_fastq,
                              file.path(str_remove(args$glamr_root, '/$'), fwd_fastq))
    ) %>%
    mutate(rev_fastq=if_else(str_starts(sample_dir, '/'),
                              rev_fastq,
                              file.path(str_remove(args$glamr_root, '/$'), rev_fastq))
    )

if (nrow(assignments) != nrow(samples)) {
    stop('the sample file listing is missing rows for some samples??')
}

# Create symlinks to our fastq files, soley so that we can have the sample id
# as part of the file name This will make the quality plotting happy, as every
# file will get its own plot.
tmpdir = tempdir()
samples = samples |>
    add_column(fwd_fastq_tmp=samples$fwd_fastq) |>
    add_column(rev_fastq_tmp=samples$rev_fastq) |>
    mutate(fwd_fastq_tmp=str_glue('{tmpdir}/{sample}_fwd.fastq.gz')) |>
    mutate(rev_fastq_tmp=str_glue('{tmpdir}/{sample}_rev.fastq.gz'))
link_create(fs::path_real(samples$fwd_fastq), samples$fwd_fastq_tmp, symbolic=TRUE)
link_create(fs::path_real(samples$rev_fastq), samples$rev_fastq_tmp, symbolic=TRUE)

forwardReads <- samples$fwd_fastq_tmp
reverseReads <- samples$rev_fastq_tmp

# extract sample names
namesOfSamples <- samples$sample


# visualize and save quality of forward and reverse reads
cat('Plotting fwd reads quality for', nrow(samples), 'samples...\n')
forwardPlot <- plotQualityProfile(forwardReads)
plotDataFwd <- forwardPlot$data
fs::dir_create(path = args$outdir)
fwd_plot_path = str_glue("{args$outdir}/forward_quality_plot.pdf")
ggsave(filename=fwd_plot_path, plot=forwardPlot, width=5, height=3, scale=2)
cat('Saved as:', fwd_plot_path, '\n')

cat('Plotting rev reads quality for', nrow(samples), 'samples...\n')
reversePlot <- plotQualityProfile(reverseReads)
plotDataRev <- reversePlot$data
rev_plot_path = str_glue("{args$outdir}/reverse_quality_plot.pdf")
ggsave(filename=rev_plot_path, plot=reversePlot, width=5, height=3, scale=2)
cat('Saved as:', rev_plot_path, '\n')

cat('Computing truncation parameters...\n')
# using weighted mean
fwd_quality <- plotDataFwd %>%
  group_by(Cycle) %>%
  mutate(mean = sum(Score * Count) / sum(Count)) %>%
  slice_max(Count, with_ties = FALSE)

rev_quality <- plotDataRev %>%
  group_by(Cycle) %>%
  mutate(mean = sum(Score * Count) / sum(Count)) %>%
  slice_max(Count, with_ties = FALSE)

threshold <- args$quality

# reads in HMM coordinates
fwd_start <- round(median(samples$fwd_hmmfrom))
fwd_end <- round(median(samples$fwd_hmmto))
rev_start <- round(median(samples$rev_hmmfrom))
rev_end <- round(median(samples$rev_hmmto))
begin_overlap_hmm <- round(median(samples$overlap_start))
end_overlap_hmm <- round(median(samples$overlap_end))
overlap_len = end_overlap_hmm - begin_overlap_hmm
cat('  hmm/fwds:', fwd_start, '-', fwd_end, '\n')
cat('  hmm/revs:', rev_start, '-', rev_end, '\n')
cat('  hmm/over:', begin_overlap_hmm, '-', end_overlap_hmm, ' (length:', overlap_len, ')', '\n')
if (overlap_len < MIN_OVERLAP) {
    stop(paste("Not enough overlap. Expected ", MIN_OVERLAP, "+. Overlap is:", overlap_len))
}
# Note: length of start and end on hmm does not necessarily correlate with length of start and end on seq

# translate into read coordinates
begin_overlap_seq_fwd = begin_overlap_hmm - fwd_start
end_overlap_seq_fwd = min(begin_overlap_seq_fwd + overlap_len, max(fwd_quality$Cycle))
begin_overlap_seq_rev = min(rev_end - begin_overlap_hmm, max(rev_quality$Cycle))
end_overlap_seq_rev = rev_end - end_overlap_hmm
cat('  begin_overlap_seq_fwd:', begin_overlap_seq_fwd, '\n')
cat('    end_overlap_seq_fwd:', end_overlap_seq_fwd, '\n')
cat('  begin_overlap_seq_rev:', begin_overlap_seq_rev, '\n')
cat('    end_overlap_seq_rev:', end_overlap_seq_rev, '\n')

# find 20-base window with maximum quality score
fwd_window = begin_overlap_seq_fwd:(begin_overlap_seq_fwd + MIN_OVERLAP - 1)
rev_window = (begin_overlap_seq_rev - MIN_OVERLAP + 1):begin_overlap_seq_rev
high_score = -1
while (last(fwd_window) < end_overlap_seq_fwd && first(rev_window) > end_overlap_seq_rev) {
    fwd_window = fwd_window + 1
    rev_window = rev_window - 1
    score = sum(
        fwd_quality$mean[fwd_window],
        rev_quality$mean[rev_window]
    )
    if (score > high_score) {
        high_score = score
        forTrunc = last(fwd_window)
        revTrunc = last(rev_window)
    }
}
cat('  Initial trunc params:: fwd:', forTrunc, 'rev:', revTrunc, '\n')

# extend to the "right"
while (mean(fwd_quality$mean[(forTrunc - 1):(forTrunc + 1)], na.rm=TRUE) > threshold) {
    if (forTrunc >= end_overlap_seq_fwd) {
      break
    }
    forTrunc <- forTrunc + 1
}

# extend to the "left
while (mean(rev_quality$mean[(revTrunc - 1):(revTrunc + 1)], na.rm=TRUE) > threshold) {
    if (revTrunc >= begin_overlap_seq_rev) {
      break
    }
    revTrunc <- revTrunc + 1
}
cat(" Extended trunc params:: fwd:", forTrunc, "rev:", revTrunc, '\n')
final_overlap = overlap_len - (end_overlap_seq_fwd - forTrunc) - (begin_overlap_seq_rev - revTrunc)
cat(" Estimated overlap:", final_overlap, '\n')

# filter and trim
cat("Filtering reads...\n")
filtAndTrimForward <- samples$filt_reads_fwd
filtAndTrimReverse <- samples$filt_reads_rev

names(filtAndTrimForward) <- namesOfSamples
names(filtAndTrimReverse) <- namesOfSamples

filterstats <- filterAndTrim(forwardReads, filtAndTrimForward, reverseReads, filtAndTrimReverse,
                     truncLen=c(forTrunc, revTrunc),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=cpus)

# saving the filter and trim to .tsv file
rownames(filterstats) = namesOfSamples
filterstats = as.data.frame(filterstats) |> rownames_to_column('sample')
write_tsv(filterstats, str_glue("{args$outdir}/filt_and_trim.tsv"))

# If all reads for a sample got filtered out, the corresponding (empty) file is not saved
filtAndTrimForward = filtAndTrimForward[file.exists(filtAndTrimForward) == TRUE]
filtAndTrimReverse = filtAndTrimReverse[file.exists(filtAndTrimReverse) == TRUE]
num_good = length(filtAndTrimForward)
if (num_good != length(filtAndTrimReverse)) { stop("fwd/rev good count mismatch"); }
if (num_good < nrow(samples)) {
    cat("[WARNING]", nrow(samples) - num_good, "samples out of", nrow(samples), "removed because (presumably) all their reads got filtered out.\n")
}

# learning errors
cat("Learning errors...\n")
errorForward <- learnErrors(filtAndTrimForward, multithread=cpus)
errorReverse <- learnErrors(filtAndTrimReverse, multithread=cpus)

# saving errors to an rds file
errors <- list(errorForward, errorReverse)
write_rds(errors, str_glue("{args$outdir}/errors.rds"))

# saving the forward error plot to a file
ggsave(
    filename=str_glue("{args$outdir}/forward_error_plot.pdf"),
    plot=plotErrors(errorForward, nominalQ=TRUE),
    width=5, height=3, scale=2,
)
ggsave(
    filename=str_glue("{args$outdir}/reverse_error_plot.pdf"),
    plot=plotErrors(errorReverse, nominalQ=TRUE),
    width=5, height=3, scale=2,
)

# sample inference
cat("Running dada...\n")
dadaForward <- dada(filtAndTrimForward, err=errorForward, multithread=cpus)
dadaReverse <- dada(filtAndTrimReverse, err=errorReverse, multithread=cpus)

# merging forward and reverse reads
cat("Merging pairs...\n")
mergeReads <- mergePairs(dadaForward, filtAndTrimForward, dadaReverse, filtAndTrimReverse, verbose=TRUE)

# construct sequence table
seqtable <- makeSequenceTable(mergeReads)

# remove chimeras
cat("Removing chimeras...\n")
seqtable_nochimeras <- removeBimeraDenovo(seqtable, method="consensus", multithread=cpus, verbose=TRUE)

# save to .tsv file
asv_ids <- paste('tasv', 1:ncol(seqtable_nochimeras), sep='')
seqtable_nochimeras |>
  as.data.frame() |>
  `colnames<-`(asv_ids) |>
  rownames_to_column('sample') |>
  write_tsv(str_glue("{args$outdir}/asv_table.tsv"))

# output representative sequences to fasta
Biostrings::DNAStringSet(colnames(seqtable_nochimeras)) |>
  `names<-`(asv_ids) |>
  Biostrings::writeXStringSet(str_glue("{args$outdir}/rep_seqs.fasta"))

# track reads through pipelines and save as .tsv
getN <- function(x) sum(getUniques(x))
fwd_counts = as.data.frame(sapply(dadaForward, getN))    |> rownames_to_column('sample')
rev_counts = as.data.frame(sapply(dadaReverse, getN))    |> rownames_to_column('sample')
merge_cnts = as.data.frame(sapply(mergeReads,  getN))    |> rownames_to_column('sample')
nochim_cts = as.data.frame(rowSums(seqtable_nochimeras)) |> rownames_to_column('sample')
track = filterstats |>
    left_join(fwd_counts, by='sample') |>
    left_join(rev_counts, by='sample') |>
    left_join(merge_cnts, by='sample') |>
    left_join(nochim_cts, by='sample')
colnames(track) <- c("sample", "input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
write_tsv(track, str_glue("{args$outdir}/track_read_counts.tsv"))
