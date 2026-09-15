#!/bin/bash
set -euo pipefail

function logto () {
    if [[ -z "$1" ]]; then
        echo >&2 "logto: Need the log file as argument"
        exit 1
    fi

    log=$1
    if [[ -n "$log" ]]; then
        exec > >(tee -i "$log" -p --output-error=exit)
        exec 2>&1
    fi
}

echo "Job executed on: ${HOSTNAME}"
echo "SLURM job id: ${SLURM_JOB_ID:-local job}"

export BASE=  # set by pypelib.shell_prep()
