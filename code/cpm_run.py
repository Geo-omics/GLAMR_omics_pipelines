#!/usr/bin/env python3
"""
Runs CPM (cpm.cpm.main, same CLI as the `cpm` command) with its figures
disabled. Several CPM plots size the canvas by the data (a legend entry per
file, half an inch per heatmap row or m/z bin) at 300 dpi; on a 624-file
dataset single figures reached ~300 MB PNGs, failed with std::bad_alloc, or
grew the process past 700 GB and were OOM-killed. The tables are unchanged.
Run inside the cpm conda env.
"""
import sys

import cpm.cpm as cpm_main
import cpm.adduct_pipeline as adduct_pipeline
import cpm.cyanometdb_match as cyanometdb_match


def skip_plot(*args, **kwargs):
    return None


for name in ('plot_rt_histograms', 'plot_indiv_counts', 'plot_precursor_rt',
             'plot_heatmaps', 'plot_indiv_scatter', 'plot_has_tilemap',
             'plot_matched_tiles'):
    setattr(cpm_main, name, skip_plot)
adduct_pipeline.draw_colored_graph = skip_plot
# write_outputs() only draws its unknown-features heatmap when seaborn imported
cyanometdb_match._HAS_SNS = False

sys.exit(cpm_main.main())
