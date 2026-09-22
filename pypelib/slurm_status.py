#!/usr/bin/env python3

# This script is taken from https://github.com/LUMC/slurm-cluster-status
#
# Copyright 2019 Leiden University Medical Center
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import argparse
import subprocess

STATE_MAP = {
    "BOOT_FAIL": "failed",
    "CANCELLED": "failed",
    "COMPLETED": "success",
    "CONFIGURING": "running",
    "COMPLETING": "running",
    "DEADLINE": "failed",
    "FAILED": "failed",
    "NODE_FAIL": "failed",
    "OUT_OF_MEMORY": "failed",
    "PENDING": "running",
    "PREEMPTED": "failed",
    "RUNNING": "running",
    "RESIZING": "running",
    "SUSPENDED": "running",
    "TIMEOUT": "failed",
    "UNKNOWN": "running",
    " ": "running",
    "": "running"
}


def fetch_status(batch_id):
    """fetch the status for the batch id"""
    sacct_args = ["sacct", "-j",  batch_id, "-o", "State", "--parsable2",
                  "--noheader"]

    try:
        output = subprocess.check_output(sacct_args).decode("utf-8").strip()
    except Exception:
        # If sacct fails for whatever reason, assume its temporary and return 'running'
        output = 'UNKNOWN'

    # The first output is the state of the overall job
    # See
    # https://stackoverflow.com/questions/52447602/slurm-sacct-shows-batch-and-extern-job-names
    # for details
    job_status = output.split("\n")[0]

    # If the job was cancelled manually, it will say by who, e.g "CANCELLED by 12345"
    # We only care that it was cancelled
    if job_status.startswith("CANCELLED by"):
        job_status = "CANCELLED"

    # Otherwise, return the status
    try:
        return STATE_MAP[job_status]
    except KeyError:
        raise NotImplementedError(f"Encountered unknown status {job_status} "
                                  f"when parsing output:\n{output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("batch_id", type=str)
    args = parser.parse_args()

    status = fetch_status(args.batch_id)
    print(status)
