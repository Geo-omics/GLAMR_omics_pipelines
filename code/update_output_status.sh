#!/bin/bash

singularity run --bind /geomicro:/geomicro,/nfs:/nfs docker://eandersk/r_microbiome ~/GLAMR/code/update_output_status.R