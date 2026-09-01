#!/usr/bin/env python3
"""
Dev-only setup script: symlinks real GLAMR metabolomics data into this repo,
read-only, so the pipeline can be evaluated against real data without ever
writing into ~/GLAMR (canonically /nfs/turbo/lsa-kiledal/GLAMR).

This is NOT part of the pipeline and is not run by Snakemake. It exists only
to populate a local dev sandbox. Rerunning it is idempotent.

What it does:
  1. For each of the 352 already-imported set_41 samples (samp_11304-11655),
     creates data/omics/metabolomes/{sample}/spectra/raw_spectra.{vendor}.raw
     as a symlink to the real underlying .raw file (resolved, not the
     multi-hop chain GLAMR uses), and mirrors the data/projects/set_41/...
     symlink. The vendor tag (e.g. "thermo_raw") is embedded in the filename
     because Snakemake needs it as a wildcard to route conversion by vendor
     - Thermo .raw is a single file, but Waters .raw is a directory, so the
     bare extension alone can't disambiguate them; see pipelines/
     metabolomics.smk stage 1.
  2. For the 272 QC/blank/standard files that exist in GLAMR's staging area
     but were never assigned sample IDs or imported, assigns new synthetic
     samp_ IDs continuing the real sequence (starting at 11656) and creates
     the same local symlink structure, so the manifest/role-guessing/QC
     stages of the pipeline can be exercised end-to-end on real data.
  3. Writes dev/test_data_manifest.tsv recording every synthetic ID mapping,
     for reference (this is a dev fixture, not a pipeline input).

All writes are local to this repo. GLAMR is only ever read from.
"""
import csv
import re
import sys
from pathlib import Path

GLAMR = Path("/nfs/turbo/lsa-kiledal/GLAMR")
REPO = Path(__file__).resolve().parent.parent
STAGING = GLAMR / "import/staging/202407_metabolomics_from_usgs"

sys.path.insert(0, str(REPO))
from pypelib.metabolomics.detect_format import detect  # noqa: E402


def make_sample(sample_id: str, real_raw_path: Path, dest_dataset: str,
                 rows: list, year: str, filename: str, ms_role_hint: str):
    real_raw_path = real_raw_path.resolve()
    assert real_raw_path.is_file(), f"missing source file: {real_raw_path}"

    fmt = detect(real_raw_path)
    link_name = f"raw_spectra.{fmt.vendor}{real_raw_path.suffix}"

    sample_dir = REPO / "data/omics/metabolomes" / sample_id / "spectra"
    sample_dir.mkdir(parents=True, exist_ok=True)
    # clean up any stale symlink from a prior naming convention
    for stale in sample_dir.glob("raw_spectra.*"):
        stale.unlink()
    link = sample_dir / link_name
    link.symlink_to(real_raw_path)

    proj_dir = REPO / "data/projects" / dest_dataset / "metabolomes"
    proj_dir.mkdir(parents=True, exist_ok=True)
    proj_link = proj_dir / sample_id
    target = Path("../../../omics/metabolomes") / sample_id
    if proj_link.is_symlink() or proj_link.exists():
        proj_link.unlink()
    proj_link.symlink_to(target)

    rows.append({
        "SampleID": sample_id,
        "StudyID": dest_dataset,
        "year": year,
        "filename": filename,
        "ms_role_hint": ms_role_hint,
        "source_path": str(real_raw_path),
    })


def guess_role(filename: str) -> str:
    stem = filename.removesuffix(".raw")
    if re.match(r"^blank\d+$", stem, re.IGNORECASE):
        return "blank"
    if re.match(r"^QC_SA", stem, re.IGNORECASE):
        return "qc_spike"
    if re.match(r"^QC_A", stem, re.IGNORECASE):
        return "qc_pool"
    if re.match(r"^LCTXM", stem, re.IGNORECASE):
        return "cal_standard"
    return "sample"


def main():
    rows = []

    # 1. The 352 already-imported samples
    import_table = STAGING / "glamr_import_table.tsv"
    with open(import_table) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            src = GLAMR / r["path"]
            make_sample(
                r["SampleID"], src, r["StudyID"], rows,
                year=r["year"], filename=r["filename"], ms_role_hint="sample",
            )
    print(f"linked {len(rows)} previously-imported samples")

    # 2. The 272 un-imported QC/blank/standard files
    newdata = STAGING / "newdata11212024"
    next_id = 11656
    n_before = len(rows)
    for raw_file in sorted(newdata.glob("*/*.raw")):
        year = raw_file.parent.name
        filename = raw_file.name
        sample_id = f"samp_{next_id}"
        next_id += 1
        make_sample(
            sample_id, raw_file, "set_41", rows,
            year=year, filename=filename, ms_role_hint=guess_role(filename),
        )
    print(f"linked {len(rows) - n_before} QC/blank/standard files as {rows[n_before]['SampleID']}..{rows[-1]['SampleID']}")

    manifest_path = REPO / "dev/test_data_manifest.tsv"
    with open(manifest_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {manifest_path} ({len(rows)} rows total)")


if __name__ == "__main__":
    main()
