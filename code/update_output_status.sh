#!/bin/bash

# Pull and run the container directly from Docker Hub (rather than a local .sif),
# so this always uses the latest published image.
#
# APPTAINER_IGNORE_PROOT / PROOT_NO_SECCOMP are required (do not remove): pulling
# and converting an OCI image to SIF uses a proot-wrapped mksquashfs, and several
# geomicro hosts have kernel.yama.ptrace_scope hardened so proot's ptrace(TRACEME)
# is denied - which makes the pull fail outright without these set.
export APPTAINER_IGNORE_PROOT=1
export PROOT_NO_SECCOMP=1

# --dns is required (do not remove): this script connects to the cayman Postgres
# server and Google Sheets, both of which need working DNS in the container.
# 127.0.0.53 adapts to whatever network DNS is  available); 1.1.1.1/1.0.0.1 (Cloudflare) are fallbacks.
# Use `exec ... Rscript` rather than `run`: the container's runscript does not
# reliably execute a passed-in .R file to completion (it would exit 0 without
# actually running the script), so invoke Rscript on it explicitly.
singularity exec \
    --dns 127.0.0.53,1.1.1.1,1.0.0.1 \
    --bind /geomicro:/geomicro,/nfs:/nfs \
    docker://eandersk/r_microbiome \
    Rscript ~/GLAMR/code/update_output_status.R