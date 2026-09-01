#!/usr/bin/env Rscript

'QC and filter an MZmine aligned feature table using the run manifest

Each step is conditional on the manifest actually containing the roles it
needs, so a project with no blanks or no QC pools still completes (with the
report noting what was skipped) instead of failing - not every future
metabolomics project will have this one QC design.

Usage:
  metabolomics_feature_qc.R --features <PATH> --manifest <PATH> --group <NAME> [options]

Options:
  -h --help                    Show this screen
  --features <PATH>            MZmine CSV feature table (feature_table.csv)
  --manifest <PATH>            run_manifest.tsv for the dataset
  --group <NAME>                processing group_id (filters manifest rows to this group)
  --blank-fold-change <N>      min sample/blank intensity ratio to keep a feature [default: 5]
  --qc-max-cv <N>               max coefficient of variation across qc_pool runs [default: 0.3]
  --min-qc-for-drift <N>        min qc_pool runs with signal needed to fit a per-feature drift curve [default: 5]
  --out-table <PATH>            output filtered feature table TSV [default: feature_table_filtered.tsv]
  --out-report <PATH>           output QC report HTML [default: qc_report.html]
' -> doc

library(docopt)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)

args <- docopt(doc)
blank_fold_change <- as.numeric(args$blank_fold_change)
qc_max_cv <- as.numeric(args$qc_max_cv)
min_qc_for_drift <- as.numeric(args$min_qc_for_drift)

message("Reading manifest: ", args$manifest)
manifest <- read_tsv(args$manifest, show_col_types = FALSE)

message("Reading feature table: ", args$features)
features <- read_csv(args$features, show_col_types = FALSE)

# MZmine's "full feature table" CSV export (CSVExportModularModule) is wide
# per METRIC per sample - "datafile:{sample}.mzML:area", ":height", ":rt",
# etc. - not a single intensity column per sample. Confirmed against a real
# feature table (2073 features x 6 samples): there is no bare "samp_11305"
# column, only "datafile:samp_11305.mzML:area" and friends. We use :area as
# the abundance measure (integrates signal over the peak; matches the
# default used by mzmine's own GNPS-FBMN export).
sample_area_col <- function(sample) paste0("datafile:", sample, ".mzML:area")

manifest <- manifest |> mutate(area_col = sample_area_col(sample))
present <- manifest$area_col %in% colnames(features)
sample_map <- manifest |> filter(present)
if (nrow(sample_map) == 0) {
    stop(
        "No feature table columns matched any sample ID in the manifest. ",
        "Expected columns named 'datafile:{sample}.mzML:area' (MZmine's ",
        "full-feature-table export naming)."
    )
}
message("Matched ", nrow(sample_map), " sample columns in the feature table")

skipped_steps <- c()
kept <- features
n_start <- nrow(kept)
all_cols <- sample_map$area_col

role_counts <- table(sample_map$ms_role)
message("Roles present in this group: ",
        paste(names(role_counts), role_counts, sep = "=", collapse = ", "))

zero_na <- function(x) replace(x, is.na(x), 0)

qc_cv_of <- function(tbl, cols) {
    if (length(cols) < 2) return(rep(NA_real_, nrow(tbl)))
    vals <- tbl |> select(all_of(cols)) |> as.matrix()
    apply(vals, 1, function(x) {
        x <- x[!is.na(x) & x > 0]
        if (length(x) < 2) return(NA_real_)
        sd(x) / mean(x)
    })
}

# --- Step 1: intensity drift correction (QC-RLSC) ---------------------------
# Corrects systematic run-order drift (ESI source sensitivity, LC performance
# decay, etc.) using the QC pool injections as a reference: fits a smooth
# per-feature trend of QC pool intensity vs. injection order (from the
# manifest's acquisition_date, real for every sample), then rescales every
# sample's intensity by that trend, centered on the QC median for that
# feature (Dunn et al. 2011-style QC-RLSC). This must run BEFORE blank
# filtering and precision filtering below - both compare intensities across
# samples acquired at different times, so they're only as meaningful as the
# drift correction that (ideally) already made those intensities comparable.
qc_cols <- sample_map$area_col[sample_map$ms_role == "qc_pool"]
qc_cv_before <- NULL
qc_cv_after <- NULL
if (length(qc_cols) >= min_qc_for_drift) {
    acq_time <- as.POSIXct(sample_map$acquisition_date, format = "%m/%d/%Y %H:%M:%S")
    if (any(is.na(acq_time))) {
        stop("Failed to parse acquisition_date for ", sum(is.na(acq_time)),
             " manifest rows - expected format MM/DD/YYYY HH:MM:SS")
    }
    injection_order <- rank(acq_time, ties.method = "first")
    names(injection_order) <- sample_map$area_col
    qc_order <- injection_order[qc_cols]
    all_order <- injection_order[all_cols]

    qc_cv_before <- qc_cv_of(kept, qc_cols)

    area_matrix <- kept |> select(all_of(all_cols)) |> as.matrix()
    qc_matrix <- area_matrix[, qc_cols, drop = FALSE]

    n_corrected <- 0L
    for (i in seq_len(nrow(area_matrix))) {
        qc_vals <- qc_matrix[i, ]
        ok <- !is.na(qc_vals) & qc_vals > 0
        if (sum(ok) < min_qc_for_drift) next  # not enough QC signal for this feature - leave as-is

        fit <- tryCatch(
            loess(y ~ x, data.frame(x = qc_order[ok], y = qc_vals[ok]), span = 0.75, degree = 1),
            error = function(e) NULL
        )
        if (is.null(fit)) next

        ref_level <- median(qc_vals[ok])
        trend_at_sample <- predict(fit, newdata = data.frame(x = all_order))
        # Guard against non-positive/NA trend values (extrapolation beyond the
        # QC injection range, or a near-zero local fit) blowing up the ratio.
        bad <- is.na(trend_at_sample) | trend_at_sample <= 0
        trend_at_sample[bad] <- ref_level

        area_matrix[i, ] <- area_matrix[i, ] / trend_at_sample * ref_level
        n_corrected <- n_corrected + 1L
    }
    kept[all_cols] <- as.data.frame(area_matrix)
    qc_cv_after <- qc_cv_of(kept, qc_cols)

    message("Drift correction: fit a QC-based trend for ", n_corrected, " / ", nrow(kept),
            " features (>= ", min_qc_for_drift, " qc_pool detections required); ",
            "median QC pool CV ", sprintf("%.3f -> %.3f", median(qc_cv_before, na.rm = TRUE),
                                           median(qc_cv_after, na.rm = TRUE)))
} else {
    skipped_steps <- c(skipped_steps,
        paste0("drift correction (need >= ", min_qc_for_drift, " qc_pool runs, found ", length(qc_cols), ")"))
    message("[SKIP] drift correction: only ", length(qc_cols), " qc_pool runs (need >= ", min_qc_for_drift, ")")
}

# --- Step 2: blank filtering -------------------------------------------------
blank_cols <- sample_map$area_col[sample_map$ms_role == "blank"]
sample_role_cols <- sample_map$area_col[sample_map$ms_role == "sample"]
if (length(blank_cols) >= 1) {
    # A feature absent from a run reads as NA in MZmine's export, not 0;
    # treat "not detected" as zero abundance so blank/sample means and the
    # keep decision below are never NA (sum(!keep) must stay a plain count).
    blank_mean <- kept |>
        select(all_of(blank_cols)) |>
        mutate(across(everything(), zero_na)) |>
        rowMeans()
    sample_mean <- kept |>
        select(all_of(sample_role_cols)) |>
        mutate(across(everything(), zero_na)) |>
        rowMeans()
    keep <- blank_mean == 0 | (sample_mean / pmax(blank_mean, 1)) >= blank_fold_change
    message("Blank filtering: removing ", sum(!keep), " / ", nrow(kept), " features ",
            "(", length(blank_cols), " blank runs, fold-change >= ", blank_fold_change, ")")
    kept <- kept[keep, ]
} else {
    skipped_steps <- c(skipped_steps, "blank filtering (no blank runs in manifest for this group)")
    message("[SKIP] blank filtering: no blank runs")
}

# --- Step 3: precision (CV) filtering via QC pools --------------------------
# Uses the drift-corrected values from step 1 (if it ran) - a feature that's
# only noisy because of uncorrected run-order drift no longer looks noisy
# here, so this now measures analytical precision rather than a mix of
# precision and trend.
if (length(qc_cols) >= 3) {
    qc_cv <- qc_cv_of(kept, qc_cols)
    keep <- is.na(qc_cv) | qc_cv <= qc_max_cv
    message("Precision filtering: removing ", sum(!keep, na.rm = TRUE), " / ", nrow(kept),
            " features (CV across ", length(qc_cols), " QC pool runs > ", qc_max_cv, ")")
    kept <- kept[keep, ]
} else {
    skipped_steps <- c(skipped_steps,
        paste0("precision filtering (need >= 3 qc_pool runs, found ", length(qc_cols), ")"))
    message("[SKIP] precision filtering: only ", length(qc_cols), " qc_pool runs")
}

message("Final: ", nrow(kept), " / ", n_start, " features retained")

write_tsv(kept, args$out_table)

# --- Minimal HTML report -----------------------------------------------------
report_lines <- c(
    "<html><head><title>Metabolomics feature QC report</title></head><body>",
    str_glue("<h1>QC report: {args$group}</h1>"),
    str_glue("<p>Input features: {n_start}, retained: {nrow(kept)}</p>"),
    str_glue("<p>Roles present: {paste(names(role_counts), role_counts, sep='=', collapse=', ')}</p>")
)
if (length(skipped_steps) > 0) {
    report_lines <- c(report_lines, "<h2>Skipped steps</h2><ul>",
        paste0("<li>", skipped_steps, "</li>"), "</ul>")
}
if (!is.null(qc_cv_before)) {
    n_improved <- sum(qc_cv_after < qc_cv_before, na.rm = TRUE)
    n_comparable <- sum(!is.na(qc_cv_before) & !is.na(qc_cv_after))
    report_lines <- c(report_lines, "<h2>Drift correction</h2>",
        str_glue("<p>Median QC pool CV: {sprintf('%.3f', median(qc_cv_before, na.rm = TRUE))} -> ",
                 "{sprintf('%.3f', median(qc_cv_after, na.rm = TRUE))} ",
                 "({n_improved} / {n_comparable} features improved)</p>"))
    cv_df <- bind_rows(
        tibble(cv = qc_cv_before, when = "before"),
        tibble(cv = qc_cv_after, when = "after")
    ) |> filter(!is.na(cv))
    if (nrow(cv_df) > 0) {
        p <- ggplot(cv_df, aes(cv, fill = when)) +
            geom_histogram(position = "identity", alpha = 0.5, bins = 60) +
            coord_cartesian(xlim = c(0, quantile(cv_df$cv, 0.95))) +
            theme_minimal() + ggtitle("QC pool CV per feature, before vs after drift correction")
        cv_png_path <- str_replace(args$out_report, "\\.html$", "_drift_cv.png")
        ggsave(cv_png_path, p, width = 7, height = 4)
        report_lines <- c(report_lines,
            str_glue('<img src="{basename(cv_png_path)}">'))
    }
}
if (length(qc_cols) >= 3) {
    pca_matrix <- kept |> select(all_of(c(qc_cols, sample_role_cols))) |> as.matrix()
    pca_matrix[is.na(pca_matrix)] <- 0
    pca_input <- t(pca_matrix)
    pca_input <- pca_input[, apply(pca_input, 2, var) > 0, drop = FALSE]
    if (ncol(pca_input) >= 2 && nrow(pca_input) >= 3) {
        pca <- prcomp(log1p(pca_input), scale. = TRUE)
        pca_df <- as_tibble(pca$x[, 1:2]) |>
            mutate(area_col = rownames(pca_input),
                   role = if_else(area_col %in% qc_cols, "qc_pool", "sample"))
        p <- ggplot(pca_df, aes(PC1, PC2, color = role)) + geom_point(size = 2) +
            theme_minimal() + ggtitle("PCA: QC pools should cluster tightly")
        png_path <- str_replace(args$out_report, "\\.html$", "_pca.png")
        ggsave(png_path, p, width = 6, height = 5)
        report_lines <- c(report_lines,
            str_glue('<h2>PCA</h2><img src="{basename(png_path)}">'))
    }
}
report_lines <- c(report_lines, "</body></html>")
writeLines(report_lines, args$out_report)
message("Wrote ", args$out_table, " and ", args$out_report)
