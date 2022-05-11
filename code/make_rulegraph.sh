#! /bin/bash

snakemake binning sample_annotation --rulegraph  --dry-run | dot -Tpng > rulegraph.png
