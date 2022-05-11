#! /bin/bash

snakemake binning metaG_annotation --rulegraph  --dry-run | dot -Tpng > rulegraph.png
