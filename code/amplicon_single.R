#!/usr/bin/env Rscript

'Using dada2 with single-ended amplicon sequencing data

Usage:
  amplicon_single.R [--glamr-root <PATH>] --quality <N> --outdir <PATH> [--cpus <N>] --assignments <PATH> --samples <PATH> --targets <PATH> <target_spec>

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


valid_targets = system2(
    'python3',
    args=shQuote(c('-m', 'pypelib.amplicon.dispatch', 'spec2targets', args$target_spec)),
    stdout=TRUE,
)
status = attr(valid_targets, 'status')
if (!is.null(status)) {
    stop(paste('called process failed with exit status ', status))
}
cat('Targets:', valid_targets, '\n')

assignments <- read.delim(
    args$assignments,
    header=TRUE,
    fill=TRUE,
) |> tibble() |>
    mutate(assigned=if_else(override == "", target, override, missing=target)) |>
    filter(assigned %in% valid_targets)

if (nrow(assignments) == 0) {
    stop('no samples where assigned to given target')
}

target_tab <- read.delim(args$targets, header=TRUE)

samples <- read.delim(
    args$samples,
    header=TRUE,
) |> tibble() |>
    semi_join(assignments, by=join_by('sample')) %>%
    left_join(target_tab, by=join_by('sample')) %>%
    mutate(filt_reads_single=str_glue('{args$outdir}/filtered/{sample}_single.fastq.gz')) %>%
    # below: paste glamr_root and paths, but respect absolute paths in data,
    # make some effort to avoid double slashes
    mutate(sample_dir=if_else(str_starts(sample_dir, '/'),
                              sample_dir,
                              file.path(str_remove(args$glamr_root, '/$'), sample_dir))
    ) |>
    mutate(single_fastq=if_else(str_starts(sample_dir, '/'),
                              single_fastq,
                              file.path(str_remove(args$glamr_root, '/$'), single_fastq))
    )

if (nrow(assignments) != nrow(samples)) {
    stop('the sample file listing is missing rows for some samples??')
}

# Create symlinks to our fastq files, soley so that we can have the sample id
# as part of the file name This will make the quality plotting happy, as every
# file will get its own plot.
tmpdir = tempdir()
samples = samples |>
    add_column(single_fastq_tmp=samples$single_fastq) |>
    mutate(single_fastq_tmp=str_glue('{tmpdir}/{sample}_single.fastq.gz'))
link_create(fs::path_real(samples$single_fastq), samples$single_fastq_tmp, symbolic=TRUE)

singleReads <- samples$single_fastq_tmp

# extract sample names
namesOfSamples <- samples$sample


# visualize and save quality of single reads
plot_cache = Sys.getenv('DADA2_REUSE_QUAL_PLOT')  # for debugging/testing only
if (file.exists(plot_cache)) {
    cat('Loading saved plot object (from', plot_cache, ')...')
    load(plot_cache)
} else {
    # env var not set / old_plot is empty str or non-existing filename
    cat('Plotting single reads quality for', nrow(samples), 'samples...')
    singlePlot <- plotQualityProfile(singleReads)
    if (plot_cache != "") {
        cat('\nSaving quality plot obj to', plot_cache)
        save(singlePlot, file=plot_cache)
    }
}
cat('[OK]\n')
fs::dir_create(path = args$outdir)
single_plot_path = str_glue("{args$outdir}/single_quality_plot.pdf")
ncols = 4
nrows = 5
num_pages = ceiling(length(singleReads) / (ncols * nrows))
cat('Saving plot as:', single_plot_path, 'with', num_pages, 'pages:')
pdf(single_plot_path, paper='default')
for (page in 1:num_pages) {
    print(singlePlot + ggforce::facet_wrap_paginate(~ file, ncol=ncols, nrow=nrows, page=page))
    cat(page, ',', sep='')
}
dev.off() -> .  # be quiet
cat('[OK]\n')


cat('Computing truncation parameters...\n')
# using weighted mean
single_quality <- singlePlot$data %>%
  group_by(Cycle) %>%
  mutate(mean = sum(Score * Count) / sum(Count)) %>%
  slice_max(Count, with_ties = FALSE)

threshold <- args$quality

# TODO: automatic truncation length or equivalent or similar

# filter and trim
cat("Filtering reads...\n")
filtAndTrimSingle <- samples$filt_reads_single

names(filtAndTrimSingle) <- namesOfSamples
filterstats <- filterAndTrim(
    singleReads, filtAndTrimSingle,
    truncLen=380,  # TODO
    maxEE=2,
    maxN=0, truncQ=2, rm.phix=TRUE,
    compress=TRUE,
    multithread=cpus,
)

# saving the filter and trim to .tsv file
rownames(filterstats) = namesOfSamples
filterstats = as.data.frame(filterstats) |> rownames_to_column('sample')
write_tsv(filterstats, str_glue("{args$outdir}/filt_and_trim.tsv"))

# If all reads for a sample got filtered out, the corresponding (empty) file is not saved
filtAndTrimSingle = filtAndTrimSingle[file.exists(filtAndTrimSingle) == TRUE]
num_good = length(filtAndTrimSingle)
if (num_good < nrow(samples)) {
    cat("[WARNING]", nrow(samples) - num_good, "samples out of", nrow(samples), "removed because (presumably) all their reads got filtered out.\n")
}

# learning errors
cat("Learning errors...\n")
errorSingle <- learnErrors(filtAndTrimSingle, multithread=cpus)

# saving errors to an rds file
#errors <- list(errorForward, errorReverse
write_rds(errorSingle, str_glue("{args$outdir}/errors.rds"))

# saving the error plot to a file
ggsave(
    filename=str_glue("{args$outdir}/single_error_plot.pdf"),
    plot=plotErrors(errorSingle, nominalQ=TRUE),
    width=5, height=3, scale=2,
)

# sample inference
cat("Running dada...\n")
dadaSingle <- dada(filtAndTrimSingle, err=errorSingle, multithread=cpus)

# construct sequence table
seqtable <- makeSequenceTable(dadaSingle)

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
single_counts = as.data.frame(sapply(dadaSingle, getN))  |> rownames_to_column('sample')
nochim_cts = as.data.frame(rowSums(seqtable_nochimeras)) |> rownames_to_column('sample')
track = filterstats |>
    left_join(single_counts, by='sample') |>
    left_join(nochim_cts, by='sample')
colnames(track) <- c("sample", "input", "filtered", "denoised", "nonchim")
write_tsv(track, str_glue("{args$outdir}/track_read_counts.tsv"))
