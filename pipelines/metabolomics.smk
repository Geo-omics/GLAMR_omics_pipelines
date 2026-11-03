"""
Metabolomics (LC-MS/MS) sub-pipeline.

Mirrors the amplicon sub-pipeline's structure (Snakefile:amplicon_* rules):
non-trivial logic lives in pypelib.metabolomics, invoked from `run:` blocks;
checkpoints handle the parts that can't be known until real data is
inspected (which run role a file has, which processing group it belongs to).

Pipeline shape, sample-level then dataset-level:

  data/omics/metabolomes/{sample}/spectra/raw_spectra.{vendor}.{ext}  (import)
    -> {sample}.mzML                          (stage 1: metab_convert_*)
    -> spectra/run_metadata.json              (stage 0: metab_harvest_metadata)

  data/projects/{dataset}/metabolomes/samp_*  (all samples in a dataset)
    -> metabolomics/run_manifest.tsv          (checkpoint metab_run_manifest)
    -> metabolomics/processing_groups.tsv     (same checkpoint)
    -> metabolomics/{group_id}/mzmine/...     (stage 2: metab_mzmine_batch)
    -> metabolomics/{group_id}/ms2query/...   (stage 3: metab_ms2query)
    -> metabolomics/{group_id}/sirius/...     (stage 3: metab_sirius, opt-in)
    -> metabolomics/{group_id}/metams/...     (stage 4: metab_metams)
    -> metabolomics/{group_id}/qc/...         (stage 5: metab_feature_qc)
    -> metabolomics_pipeline_done             (stage 6: top rule)

A "processing group" is (method_profile x polarity), the unit of MZmine
alignment - see pypelib/metabolomics/groups.py. For a single-method,
single-polarity project like set_41 this collapses to exactly one group.
"""
from pypelib.metabolomics.detect_format import detect as metab_detect_format
from pypelib.metabolomics.groups import get_all_group_ids as metab_get_group_ids
from pypelib.metabolomics.groups import get_group_samples as metab_get_group_samples
import pypelib.metabolomics.groups
import pypelib.metabolomics.harvest_metadata
import pypelib.metabolomics.run_manifest
import pypelib.metabolomics.rt_correction
import pypelib.metabolomics.rt_align
import shlex


wildcard_constraints:
    dataset = r'[^/]+',
    sample = r'samp_[A-Za-z0-9]+',
    group_id = r'[^/]+',


def get_method_profile(dataset):
    """
    Return the method_profile name assigned to a dataset.

    To add support for a new project's instrument/method combination, add an
    entry to config.yaml:

        metabolomics:
          method_profiles:
            qtof_rp_neg:
              mzbatch: config/mzmine/qtof_rp_neg.mzbatch
              corems_params: config/metams/qtof_rp_neg.toml
              scan_translator: config/metams/generic_no_filter.toml
          profile_assignment:
            my_new_study: qtof_rp_neg

    Falls back to config['metabolomics']['default_profile'] if the dataset
    has no explicit assignment (fine for single-method GLAMR projects; a
    dataset spanning multiple instruments/methods should be split into
    separate StudyIDs upstream instead, matching how GLAMR already scopes
    StudyID to one acquisition campaign).
    """
    metab_config = config.get('metabolomics', {})
    assignment = metab_config.get('profile_assignment', {})
    profile = assignment.get(dataset, metab_config.get('default_profile'))
    if profile is None:
        raise ValueError(
            f'no method_profile assigned to dataset {dataset!r} and no '
            f'metabolomics.default_profile configured - add an entry to '
            f"config.yaml's metabolomics.profile_assignment"
        )
    if profile not in metab_config.get('method_profiles', {}):
        raise ValueError(
            f'dataset {dataset!r} is assigned method_profile {profile!r}, '
            f'which has no entry in metabolomics.method_profiles'
        )
    return profile


def get_profile_param(dataset, key):
    profile = get_method_profile(dataset)
    return config['metabolomics']['method_profiles'][profile][key]


def get_group_profile_param(dataset, group_id, key):
    """
    Like get_profile_param, but a profile can override a key per polarity
    ('{key}_negative', e.g. mzbatch_negative) for datasets acquired in both
    modes; group_id is '{profile}_{polarity}'.
    """
    params = config['metabolomics']['method_profiles'][get_method_profile(dataset)]
    polarity = group_id.rsplit('_', 1)[-1]
    return params.get(f'{key}_{polarity}', params[key])


# ---------------------------------------------------------------------------
# Stage 0: metadata harvest + run manifest + processing groups
# ---------------------------------------------------------------------------

def get_sample_mzml_input(wc):
    """ Input function: the converted mzML for a sample, regardless of vendor """
    return f'data/omics/metabolomes/{wc.sample}/spectra/{wc.sample}.mzML'


def get_sample_thermo_raw(wc):
    """
    Params function: path to the original Thermo raw file for metadata
    enrichment, or None for other vendors (see harvest_metadata.py - the
    mzML stays authoritative either way, this only adds richer fields).
    """
    p = Path(f'data/omics/metabolomes/{wc.sample}/spectra/raw_spectra.thermo_raw.raw')
    return str(p) if p.exists() else None


rule metab_harvest_metadata:
    input:
        mzml = get_sample_mzml_input,
    output:
        meta = "data/omics/metabolomes/{sample}/spectra/run_metadata.json",
    params:
        thermo_raw = get_sample_thermo_raw,
    resources: mem_mb=2000, time_min=15
    benchmark: "benchmarks/metab_harvest_metadata/{sample}.txt"
    log: "logs/metab_harvest_metadata/{sample}.log"
    conda: "../config/conda_yaml/thermorawfileparser.yaml"
    run:
        pypelib.metabolomics.harvest_metadata.main(
            input.mzml,
            output.meta,
            thermo_raw=params.thermo_raw,
            conda_env=conda_env,
            log=log,
        )


def get_dataset_metadata_files(wc):
    """ All harvested metadata files for a dataset's metabolome samples """
    project_samples = Path("data/projects/") / wc.dataset / "metabolomes"
    if not project_samples.is_dir():
        return []
    return [
        f"data/omics/metabolomes/{i.name}/spectra/run_metadata.json"
        for i in sorted(project_samples.glob("samp_*"))
    ]


checkpoint metab_run_manifest:
    """
    Build the run manifest (auto metadata + auto-guessed + human-curated
    roles) and the processing-group assignment for a dataset.

    The manifest is update()-wrapped: rerunning after a human has edited
    ms_role_override/prep_variant/cal_level/BioSampleID/amount_basis
    preserves those edits (see pypelib.metabolomics.run_manifest).
    """
    input:
        get_dataset_metadata_files,
    output:
        manifest = update("data/projects/{dataset}/metabolomics/run_manifest.tsv"),
        groups = "data/projects/{dataset}/metabolomics/processing_groups.tsv",
    params:
        project_dir = "data/projects/{dataset}",
        method_profile = lambda wc: get_method_profile(wc.dataset),
    resources: mem_mb=500, time_min=10
    log: "logs/metab_run_manifest/{dataset}.log"
    run:
        import csv
        rows = pypelib.metabolomics.run_manifest.main(
            input, params.project_dir, output=output.manifest, log=log,
        )
        groups = pypelib.metabolomics.groups.assign(rows, params.method_profile)
        pypelib.metabolomics.groups.write_groups(groups, output.groups)


# ---------------------------------------------------------------------------
# Stage 1: conversion to mzML (format dispatch)
# ---------------------------------------------------------------------------
# Two rules producing the same output, disambiguated by which one's input
# actually exists - the same pattern the amplicon/metaT rules use elsewhere
# in this Snakefile (`ruleorder: prodigal_metaT > prodigal`), needed here
# because Snakemake directives (conda:/container:) are static per rule and
# Thermo (native, no Wine) vs. everything else (needs ProteoWizard/Wine)
# genuinely require different execution environments.

ruleorder: metab_convert_thermo > metab_mzml_fixup

rule metab_convert_thermo:
    """ Thermo .raw -> mzML via ThermoRawFileParser. No Wine, no container. """
    input:
        raw = "data/omics/metabolomes/{sample}/spectra/raw_spectra.thermo_raw.raw",
    output:
        mzml = "data/omics/metabolomes/{sample}/spectra/{sample}.mzML",
    params:
        base = lambda wc, output: str(Path(output.mzml).with_suffix("")),
    threads: 2
    resources: mem_mb=4000, time_min=30
    benchmark: "benchmarks/metab_convert_thermo/{sample}.txt"
    log: "logs/metab_convert_thermo/{sample}.log"
    conda: "../config/conda_yaml/thermorawfileparser.yaml"
    shell:
        """
        . code/shell_prelude {log}
        ThermoRawFileParser -i {input.raw} -b {params.base} -f 2 --logging 2
        """


def get_vendor_raw_file(wc):
    """
    Input function for metab_convert_vendor: glob for the vendor-tagged raw
    file/dir, excluding the vendors handled by the other two conversion
    rules (thermo_raw -> metab_convert_thermo, mzml -> metab_mzml_fixup directly).
    mzXML is converted here too: the pipeline expects real mzML (rt_align
    rewrites it, harvest_metadata parses its header), so it can't just be
    linked under an .mzML name.

    NOTE: {vendor} is deliberately NOT a rule wildcard here (an earlier
    version tried that and failed at DAG-build time - a wildcard that only
    appears in an input function, never in any rule's output pattern,
    cannot be bound by Snakemake's wildcard matching, since matching only
    happens against output paths). Instead this function alone determines
    the vendor via glob + exclusion, exactly like get_mzml_fixup_input does
    for its own case - {sample} is the only wildcard
    either rule needs.
    """
    all_raw = list(Path(f"data/omics/metabolomes/{wc.sample}/spectra").glob("raw_spectra.*"))
    candidates = [
        p for p in all_raw
        if not p.name.startswith("raw_spectra.thermo_raw.")
        and not p.name.startswith("raw_spectra.mzml.")
    ]
    if len(candidates) != 1:
        raise ValueError(
            f"expected exactly one non-Thermo, non-mzML raw_spectra.* for "
            f"sample {wc.sample}, found {len(candidates)} (of {len(all_raw)} total)"
        )
    return str(candidates[0])


rule metab_convert_vendor:
    """
    Agilent .d / Bruker .d / Sciex .wiff(2) / Shimadzu .lcd / Waters .raw /
    mzXML -> mzML via ProteoWizard msconvert (Wine).

    NOTE: the chambm/pwiz-skyline-i-agree-to-the-vendor-licenses image name
    is itself a license acceptance notice - see the container line below.
    ProteoWizard's Wine-based msconvert does not support ABI T2D, Bruker
    FID/YEP (bare, not in a .d dir), or Waters UNIFI; detect_format.py fails
    loudly on those rather than reaching this rule with an unconvertible
    input. Sciex .wiff conversion under Wine has a history of problems and
    should be validated against real data before production use.
    """
    input:
        raw = get_vendor_raw_file,
    output:
        mzml = temp("data/omics/metabolomes/{sample}/spectra/{sample}.msconvert.mzML"),
    params:
        outdir = lambda wc, output: str(Path(output.mzml).parent),
    threads: 2
    resources: mem_mb=8000, time_min=60
    benchmark: "benchmarks/metab_convert_vendor/{sample}.txt"
    log: "logs/metab_convert_vendor/{sample}.log"
    container: "docker://chambm/pwiz-skyline-i-agree-to-the-vendor-licenses"
    shell:
        """
        . code/shell_prelude {log}
        # no --filter: msconvert 3.0.26056 under Wine crashes mid-file with
        # "defaultArrayLength" (the empty scans it was meant to drop are
        # handled in metab_mzml_fixup instead)
        wine msconvert {input.raw:q} \
            --mzML --outdir {params.outdir:q} \
            --outfile {wildcards.sample}.msconvert.mzML
        """


def get_mzml_fixup_input(wc):
    """ Imported mzML as is, otherwise the msconvert output """
    imported = list(Path(f"data/omics/metabolomes/{wc.sample}/spectra").glob("raw_spectra.mzml.*"))
    if imported:
        return str(imported[0])
    return f"data/omics/metabolomes/{wc.sample}/spectra/{wc.sample}.msconvert.mzML"


rule metab_mzml_fixup:
    """
    msconvert output or imported mzML -> final mzML: adds the <run>
    defaultInstrumentConfigurationRef that msconvert leaves out for some
    inputs (seen with mzXML) and MZmine requires, and empty m/z / intensity
    arrays to scans with no peaks (written with no arrays at all, which
    MassQL/pyteomics can't read), rebuilding the index
    (pypelib.metabolomics.mzml_rewrite). Imported mzML that needs neither is
    copied unchanged.
    """
    input:
        mzml = get_mzml_fixup_input,
    output:
        mzml = "data/omics/metabolomes/{sample}/spectra/{sample}.mzML",
    resources: mem_mb=1000, time_min=15
    log: "logs/metab_mzml_fixup/{sample}.log"
    run:
        from pypelib.metabolomics.mzml_rewrite import fixup_msconvert_mzml
        fixup_msconvert_mzml(input.mzml, output.mzml)


# ---------------------------------------------------------------------------
# Reference data
# ---------------------------------------------------------------------------

rule metab_download_gnps_msp:
    """ GNPS-curated MS2 library, filtered to high-quality/high-res spectra by NMDC """
    output:
        msp = "data/reference/metabolomics/gnps_curated.msp",
    resources: mem_mb=500, time_min=30, heavy_network=1
    log: "logs/metab_download_gnps_msp/download.log"
    shell:
        """
        . code/shell_prelude {log}
        mkdir -p $(dirname {output.msp:q})
        curl -fL -o {output.msp:q} \
            "https://nmdcdemo.emsl.pnnl.gov/metabolomics/databases/20250407_gnps_curated.msp"
        """


rule metab_download_ms2query_models:
    """ MS2Query's pretrained Spec2Vec/MS2DeepScore/random-forest models, once per polarity """
    output:
        marker = touch("data/reference/metabolomics/ms2query_models/{polarity}/.models_downloaded"),
    wildcard_constraints:
        polarity = "positive|negative",
    params:
        outdir = lambda wc, output: str(Path(output.marker).parent),
    resources: mem_mb=4000, time_min=60, heavy_network=1
    log: "logs/metab_download_ms2query_models/{polarity}.log"
    conda: "../config/conda_yaml/ms2query.yaml"
    shell:
        """
        . code/shell_prelude {log}
        python code/ms2query_download_models.py \
            --ion-mode {wildcards.polarity} \
            --outdir {params.outdir:q} \
            --marker {output.marker:q}
        """


# ---------------------------------------------------------------------------
# Stage 2: MZmine - aligned feature matrix (per processing group)
# ---------------------------------------------------------------------------

localrules: metab_mzmine_check_user

rule metab_mzmine_check_user:
    """
    Preflight: fail fast on the head node if the mzmine login token is
    missing, rather than after hours of queued cluster compute. This is
    best-effort - it only checks the .mzuser file exists and is non-empty;
    mzmine itself validates the actual login server-side during the batch
    run (confirmed by testing: mzmine's CLI parses/loads a full batch file
    without any login check, and only fails with "Invalid user: You need to
    be logged in to continue" once it starts real data processing).
    """
    output:
        touch("data/.mzmine_user_checked"),
    run:
        mzuser = config.get('metabolomics', {}).get('mzuser')
        if not mzuser:
            raise ValueError(
                "config['metabolomics']['mzuser'] is not set - generate a "
                "token with 'mzmine -login-console' or 'mzmine -login' and "
                "point this at the resulting .mzuser file"
            )
        p = Path(mzuser)
        if not p.is_file() or p.stat().st_size == 0:
            raise ValueError(f"mzuser file missing or empty: {p}")


def get_rt_method():
    """
    config metabolomics.rt_correction.method:
      "landmark" (default) - pypelib.metabolomics.rt_align rewrites each run's
                 RTs before mzmine (rules metab_rt_landmarks / metab_rt_align),
                 mzmine then runs the plain batch file
      "mzmine"   - mzmine's own scan RT correction (the *_rtcorr batch file),
                 when rt_correction.eligible() says so
      "none"     - no RT correction
    """
    method = config.get('metabolomics', {}).get('rt_correction', {}).get('method', 'landmark')
    if method not in ('landmark', 'mzmine', 'none'):
        raise ValueError(f'unknown metabolomics.rt_correction.method {method!r}')
    return method


def get_rt_landmark_param(dataset, key, default):
    """
    Landmark detection settings (method profile `rt_landmarks:` block), for
    instruments whose intensity scale or mass accuracy differs from the
    Orbitrap defaults in pypelib.metabolomics.rt_align
    """
    profile = config['metabolomics']['method_profiles'][get_method_profile(dataset)]
    return profile.get('rt_landmarks', {}).get(key, default)


def use_mzmine_rtcorr(wc):
    return get_rt_method() == 'mzmine' and get_rt_correction_decision(wc)[0]


def get_group_mzml_files(wc):
    """ mzML inputs for mzmine: RT-aligned copies in "landmark" mode, the originals otherwise """
    manifest = checkpoints.metab_run_manifest.get(dataset=wc.dataset).output.groups
    samples = metab_get_group_samples(manifest, wc.group_id)
    if get_rt_method() == 'landmark':
        return [f"data/projects/{wc.dataset}/metabolomics/{wc.group_id}/rt_align/mzml/{s}.mzML"
                for s in samples]
    return [f"data/omics/metabolomes/{s}/spectra/{s}.mzML" for s in samples]


def _read_group_manifest_rows(manifest_path, groups_path, group_id):
    """ This processing group's rows from run_manifest.tsv, as dicts """
    samples = set(metab_get_group_samples(groups_path, group_id))
    import csv
    with open(manifest_path) as ifile:
        return [r for r in csv.DictReader(ifile, delimiter='\t') if r['sample'] in samples]


def get_group_manifest_rows(wc):
    """
    params-lambda version of _read_group_manifest_rows: resolves the
    checkpoint itself from wildcards. Rules that already have the manifest/
    groups paths as resolved `input.*` (e.g. metab_mzmine_metadata's `run:`
    block) should call _read_group_manifest_rows directly instead, to avoid
    hitting the checkpoint a second time.
    """
    manifest_path = checkpoints.metab_run_manifest.get(dataset=wc.dataset).output.manifest
    groups_path = checkpoints.metab_run_manifest.get(dataset=wc.dataset).output.groups
    return _read_group_manifest_rows(manifest_path, groups_path, wc.group_id)


def get_rt_reference_config():
    """ (reference_mode, references_per_batch) from config, see rt_correction.py """
    rt_config = config.get('metabolomics', {}).get('rt_correction', {})
    return (
        rt_config.get('reference', pypelib.metabolomics.rt_correction.DEFAULT_REFERENCE_MODE),
        rt_config.get('references_per_batch',
                      pypelib.metabolomics.rt_correction.DEFAULT_REFERENCES_PER_BATCH),
    )


def get_rt_correction_decision(wc):
    """
    Returns (eligible: bool, reason: str) - see pypelib.metabolomics.
    rt_correction for the real numbers behind min_batches/min_qc_per_batch
    and why they're configured this way.
    """
    rt_config = config.get('metabolomics', {}).get('rt_correction', {})
    min_batches = rt_config.get(
        'min_batches', pypelib.metabolomics.rt_correction.DEFAULT_MIN_BATCHES)
    min_qc_per_batch = rt_config.get(
        'min_qc_per_batch', pypelib.metabolomics.rt_correction.DEFAULT_MIN_QC_PER_BATCH)
    rows = get_group_manifest_rows(wc)
    return pypelib.metabolomics.rt_correction.eligible(
        rows, min_batches, min_qc_per_batch, get_rt_reference_config()[0])


def get_group_run_metadata(wc):
    """ run_metadata.json for every sample in a processing group (RT-reference quality) """
    groups = checkpoints.metab_run_manifest.get(dataset=wc.dataset).output.groups
    return [f"data/omics/metabolomes/{s}/spectra/run_metadata.json"
            for s in metab_get_group_samples(groups, wc.group_id)]


def _group_run_quality(run_meta_files):
    """ sample -> MS2 scan count from run_metadata.json (used to avoid failed injections) """
    import json
    quality = {}
    for f in run_meta_files:
        meta = json.load(open(f))
        if meta.get('thermo_ms2_scans'):
            quality[Path(f).parents[1].name] = int(meta['thermo_ms2_scans'])
    return quality


def get_rt_reference_mzmls(wc):
    """ Original mzMLs of the runs used to build RT landmarks (see rt_correction.select_reference_samples) """
    rows = get_group_manifest_rows(wc)
    quality = _group_run_quality(get_group_run_metadata(wc))
    refs = pypelib.metabolomics.rt_correction.select_reference_samples(
        rows, get_rt_reference_config()[1], quality)
    return [f"data/omics/metabolomes/{s}/spectra/{s}.mzML" for s in sorted(refs)]


rule metab_rt_landmarks:
    """
    Landmark ions for RT alignment of one processing group: intense,
    unambiguous ions found in >= 80% of the reference runs (a few real
    samples per acquisition batch, spread over injection order), with their
    median RT across those runs as the alignment target. See
    pypelib.metabolomics.rt_align for why this replaces mzmine's own RT
    correction by default.
    """
    input:
        refs = get_rt_reference_mzmls,
        run_meta = get_group_run_metadata,
    output:
        landmarks = "data/projects/{dataset}/metabolomics/{group_id}/rt_align/landmarks.tsv",
    resources: mem_mb=16000, time_min=60
    benchmark: "benchmarks/metab_rt_landmarks/{dataset}_{group_id}.txt"
    log: "logs/metab_rt_landmarks/{dataset}_{group_id}.log"
    params:
        min_height = lambda wc: get_rt_landmark_param(wc.dataset, 'min_height', pypelib.metabolomics.rt_align.DEFAULT_MIN_HEIGHT),
        ppm = lambda wc: get_rt_landmark_param(wc.dataset, 'ppm', pypelib.metabolomics.rt_align.DEFAULT_PPM),
    run:
        pypelib.metabolomics.rt_align.landmarks_main(
            input.refs, output.landmarks, min_height=params.min_height, ppm=params.ppm, log=log)


rule metab_rt_align:
    """
    Fit one run's RT offset vs RT against the group's landmarks (robust
    LOESS) and write an RT-corrected copy of its mzML (index and checksum
    rebuilt) plus the fit statistics / correction curve as JSON.
    """
    input:
        mzml = get_sample_mzml_input,
        landmarks = rules.metab_rt_landmarks.output.landmarks,
    output:
        mzml = "data/projects/{dataset}/metabolomics/{group_id}/rt_align/mzml/{sample}.mzML",
        stats = "data/projects/{dataset}/metabolomics/{group_id}/rt_align/stats/{sample}.json",
    resources: mem_mb=4000, time_min=20
    benchmark: "benchmarks/metab_rt_align/{dataset}_{group_id}_{sample}.txt"
    log: "logs/metab_rt_align/{dataset}_{group_id}_{sample}.log"
    params:
        min_height = lambda wc: get_rt_landmark_param(wc.dataset, 'min_height', pypelib.metabolomics.rt_align.DEFAULT_MIN_HEIGHT),
        ppm = lambda wc: get_rt_landmark_param(wc.dataset, 'ppm', pypelib.metabolomics.rt_align.DEFAULT_PPM),
    run:
        pypelib.metabolomics.rt_align.align_main(
            input.mzml, input.landmarks, output.mzml, output.stats,
            min_height=params.min_height, ppm=params.ppm, log=log)


rule metab_mzmine_metadata:
    """
    Write the mzmine project-metadata TSV (sample type + run date per file)
    for one processing group. Always generated (cheap - one row per sample),
    even for groups that turn out not to be RT-correction eligible; whether
    metab_mzmine_batch actually passes it to mzmine via --metadata is a
    separate, real per-run decision (get_rt_correction_decision) made at
    that rule's params time so a group's eligibility can change (e.g. more
    samples/batches added later) without needing to touch this rule.

    The sample type column also decides which runs mzmine uses as RT
    references (config metabolomics.rt_correction.reference: "qc" or
    "samples"; see rt_correction.py). In "samples" mode the MS2 scan count
    from each run's run_metadata.json is used to avoid failed injections.
    """
    input:
        manifest = lambda wc: checkpoints.metab_run_manifest.get(dataset=wc.dataset).output.manifest,
        groups = lambda wc: checkpoints.metab_run_manifest.get(dataset=wc.dataset).output.groups,
        run_meta = get_group_run_metadata,
    output:
        tsv = "data/projects/{dataset}/metabolomics/{group_id}/mzmine/metadata.tsv",
    params:
        reference_mode = lambda wc: get_rt_reference_config()[0],
        references_per_batch = lambda wc: get_rt_reference_config()[1],
    resources: mem_mb=500, time_min=10
    log: "logs/metab_mzmine_metadata/{dataset}_{group_id}.log"
    run:
        rows = _read_group_manifest_rows(input.manifest, input.groups, wildcards.group_id)
        quality = _group_run_quality(input.run_meta)
        with open(log[0], 'w') as logfile:
            n, refs = pypelib.metabolomics.rt_correction.write_mzmine_metadata(
                rows, output.tsv, params.reference_mode, params.references_per_batch, quality)
            print(f'wrote {n} rows to {output.tsv}', file=logfile)
            print(f'RT reference mode {params.reference_mode!r}: {len(refs)} reference runs: '
                  f'{", ".join(sorted(refs))}', file=logfile)


rule metab_mzmine_batch:
    """
    Align all mzML files in a processing group with a single mzmine JVM
    (alignment is inherently cross-sample, so this is one job per group,
    not one job per sample). The .mzbatch file is generated once via the
    mzmine wizard for the group's method_profile and committed to
    config/mzmine/ - see get_profile_param(dataset, 'mzbatch').

    mzTab-M export is NOT a required output. Confirmed by real testing: mzmine's
    MZTabmExportModule can throw a jmztab2 MZTabErrorOverflowException
    ("number of Data line items" != header count) from its own internal
    write-then-validate step on unidentified features - a real upstream
    export-module bug (still present as of 4.10.6), not something wrong with
    our batch parameters (the same run's CSV and SIRIUS MGF exports, which
    carry the same feature data, both succeed). The batch runs mzTab-M export
    LAST, after CSV and MGF are already safely on disk, and the shell block
    tolerates mzmine exiting non-zero as long as those two required outputs
    exist.

    Runs mzmine via the eandersk/r_microbiome container (mzmine 4.10.6),
    not conda - bioconda's mzmine package tops out at 4.7.29 (confirmed via
    `conda search -c bioconda mzmine`; quay.io/biocontainers, built from the
    same recipe, is capped identically). mzmine is on the container's PATH
    (confirmed via testing - `which mzmine` resolves).
    Re-validated on the same 94-sample subset used throughout this pipeline's
    development: identical feature count (5478), identical column set (1713
    columns, reordered but name-identical - harmless since every consumer
    reads this CSV by header name), identical MS2 spectra count in the MGF
    (7502) versus the prior 4.7.x/conda run. CLI flag names are unchanged
    between versions (-user/-batch/-input/-output/-temp/-memory/-threads all
    still accepted).

    NOTE: this container's `latest` tag is mutable and Snakemake's local
    singularity/apptainer image cache (.snakemake/singularity/*.simg) is
    keyed by an md5 hash of the docker:// URL string, not the image digest -
    confirmed by real testing that a container update does NOT invalidate
    this cache. After any update to eandersk/r_microbiome, delete the stale
    cached file (md5("docker://eandersk/r_microbiome") identifies it) before
    relying on new behavior, or every rule using this container - not just
    this one - will silently keep running the old image.

    RT correction: by default (metabolomics.rt_correction.method:
    "landmark") the mzML inputs are already RT-aligned by rule
    metab_rt_align and the plain batch file is used. The rest of this
    paragraph applies to method "mzmine" only. There, whether this group's alignment includes mzmine's
    scan-based RT correction step is decided per-run by
    get_rt_correction_decision (pypelib.metabolomics.rt_correction), which
    counts distinct batch_id values (see run_manifest.py's
    batch_id_guess/batch_id_override) and each batch's qc_pool coverage. Real
    measurement on set_41 (mzmine 4.10.6, drift-split diagnostic on the
    aligned feature table): 1-2 batches makes correction slightly worse
    (candidate drift-split rate +15%/+5%), 3 batches is roughly neutral
    (-4%), 4-5 batches is a real win (-26%/-39%) - see rt_correction.py's
    module docstring for the full numbers. metabolomics.rt_correction.
    min_batches in config.yaml controls the cutoff (currently 3, the point
    where it stops being harmful - this is a real config knob, not something
    to hardcode further). Uses config['metabolomics']['method_profiles']
    [profile]['mzbatch_rtcorr'] when eligible, ['mzbatch'] otherwise; the
    metadata TSV (metab_mzmine_metadata) is generated unconditionally
    either way since it's cheap, but only passed to mzmine via --metadata
    when RT correction is actually applied.
    """
    input:
        mzml = get_group_mzml_files,
        user_check = "data/.mzmine_user_checked",
        metadata = "data/projects/{dataset}/metabolomics/{group_id}/mzmine/metadata.tsv",
    output:
        # Filenames confirmed by real testing: mzmine's -output flag sets a
        # base name, but each export module appends its own configured
        # suffix before the extension (CSVExportModularParameters defaults
        # to "_full_feature_table", SiriusExportParameters to "_sirius") -
        # it does NOT produce exactly "{base}.csv"/"{base}.mgf" as the CLI
        # docs alone would suggest.
        csv = "data/projects/{dataset}/metabolomics/{group_id}/mzmine/features_full_feature_table.csv",
        mgf = "data/projects/{dataset}/metabolomics/{group_id}/mzmine/features_sirius.mgf",
    params:
        mzbatch = lambda wc: get_group_profile_param(
            wc.dataset, wc.group_id,
            'mzbatch_rtcorr' if use_mzmine_rtcorr(wc) else 'mzbatch',
        ),
        metadata_flag = lambda wc, input: (
            # NOTE: :q is a Snakemake shell-string-only quoting feature, not
            # valid in a plain f-string (confirmed by testing - raises
            # ValueError) - use shlex.quote() directly instead.
            f"--metadata {shlex.quote(input.metadata)}"
            if use_mzmine_rtcorr(wc) else ""
        ),
        mzuser = lambda wc: config['metabolomics']['mzuser'],
        outbase = lambda wc, output: str(Path(output.csv).parent / "features"),
        input_list = lambda wc, output: str(Path(output.csv).parent / "mzml_input_list.txt"),
        tmpdir = lambda wc: config['metabolomics'].get(
            'mzmine_temp_dir', f"/tmp/mzmine_{wc.dataset}_{wc.group_id}"
        ),
    threads: 8
    resources: mem_mb=64000, time_min=720
    benchmark: "benchmarks/metab_mzmine_batch/{dataset}_{group_id}.txt"
    log: "logs/metab_mzmine_batch/{dataset}_{group_id}.log"
    container: "docker://eandersk/r_microbiome"
    shell:
        """
        . code/shell_prelude {log}
        mkdir -p {params.tmpdir:q} $(dirname {params.outbase:q})
        printf '%s\n' {input.mzml:q} > {params.input_list:q}
        mzmine \
            -user {params.mzuser:q} \
            -batch {params.mzbatch:q} \
            -input {params.input_list:q} \
            {params.metadata_flag} \
            -output {params.outbase:q} \
            -temp {params.tmpdir:q} \
            -memory none \
            -threads {threads} \
        || echo "[WARNING] mzmine exited non-zero - tolerated if it was only the known mzTab-M export bug (see rule docstring); checking required outputs below"

        test -s {output.csv:q} && test -s {output.mgf:q}
        """


# ---------------------------------------------------------------------------
# Stage 3: annotation
# ---------------------------------------------------------------------------

rule metab_ms2query:
    """
    Offline MS2 analogue/exact-match search - the default annotation engine
    (no account, deterministic; see metab_sirius for the opt-in alternative
    that needs a SIRIUS account and live internet).
    """
    input:
        mgf = rules.metab_mzmine_batch.output.mgf,
        models = lambda wc: (
            f"data/reference/metabolomics/ms2query_models/"
            f"{wc.group_id.rsplit('_', 1)[-1]}/.models_downloaded"
        ),
    output:
        results = "data/projects/{dataset}/metabolomics/{group_id}/ms2query/results.csv",
    params:
        models_dir = lambda wc, input: str(Path(input.models).parent),
        outdir = lambda wc, output: str(Path(output.results).parent),
    threads: 4
    resources: mem_mb=16000, time_min=240
    benchmark: "benchmarks/metab_ms2query/{dataset}_{group_id}.txt"
    log: "logs/metab_ms2query/{dataset}_{group_id}.log"
    conda: "../config/conda_yaml/ms2query.yaml"
    shell:
        """
        . code/shell_prelude {log}
        python code/ms2query_annotate.py \
            --mgf {input.mgf:q} \
            --models-dir {params.models_dir:q} \
            --outdir {params.outdir:q}
        """


rule metab_sirius:
    """
    Opt-in (config['metabolomics']['run_sirius']) - formula/ZODIAC/
    CSI:FingerID/CANOPUS. Needs a SIRIUS account and live internet from the
    compute node; not gated behind metab_mzmine_check_user-style preflight
    because the account is separate from the mzmine one and this rule is
    already opt-in.

    Verified against real data (5 real features, real CSI:FingerID/CANOPUS
    webservice calls, real login): write-summaries does NOT produce a file
    named "compound_identifications.tsv" as an earlier, untested version of
    this rule assumed - real output filenames are formula_identifications.tsv,
    structure_identifications.tsv (CSI:FingerID structure search - has InChI/
    InChIKey, the actual compound-level annotation), denovo_structure_
    identifications.tsv, canopus_formula_summary.tsv, canopus_structure_
    summary.tsv, spectral_matches.tsv, spectral_matches_analog.tsv. This rule
    declares structure_identifications.tsv as the required output since it's
    the one with real structure predictions; the others are still written
    alongside it in the same summaries/ directory.
    """
    input:
        mgf = rules.metab_mzmine_batch.output.mgf,
    output:
        summary = "data/projects/{dataset}/metabolomics/{group_id}/sirius/summaries/structure_identifications.tsv",
    params:
        sirius_bin = lambda wc: config['metabolomics']['sirius_bin'],
        project_dir = lambda wc, output: str(Path(output.summary).parents[1] / "project"),
        summaries_dir = lambda wc, output: str(Path(output.summary).parent),
    threads: 8
    resources: mem_mb=16000, time_min=480
    benchmark: "benchmarks/metab_sirius/{dataset}_{group_id}.txt"
    log: "logs/metab_sirius/{dataset}_{group_id}.log"
    conda: "../config/conda_yaml/sirius.yaml"
    shell:
        """
        . code/shell_prelude {log}
        # SIRIUS appends to an existing project, so a rerun after a killed
        # run would import every feature a second time
        rm -rf {params.project_dir:q} {params.project_dir:q}.sirius
        {params.sirius_bin:q} \
            --input {input.mgf:q} --project {params.project_dir:q} \
            formulas fingerprints classes structures --database pubchem \
            write-summaries --output {params.summaries_dir:q}
        """


# ---------------------------------------------------------------------------
# Stage 4: metaMS - parallel per-file NMDC-comparable annotated tables
# ---------------------------------------------------------------------------

def get_sample_thermo_or_mzml_file(wc):
    """
    metaMS accepts Thermo .raw natively (no conversion needed) and mzML for
    everything else - prefer the raw file when available to skip the extra
    conversion round-trip.
    """
    thermo = Path(f"data/omics/metabolomes/{wc.sample}/spectra/raw_spectra.thermo_raw.raw")
    return str(thermo) if thermo.exists() else f"data/omics/metabolomes/{wc.sample}/spectra/{wc.sample}.mzML"


rule metab_metams_one:
    """
    Runs microbiomedata/metaMS on a single sample (not aligned - see mzmine
    for that). Bypasses miniwdl and calls the same CLI command it wraps
    directly.

    One job per sample, not per group: confirmed by reading
    metaMS/lcms_metabolomics_workflow.py upstream that metaMS does no
    cross-sample step whatsoever (each file gets its own independent
    process_complete_workflow() call and output path - only a one-time,
    read-only MSP-library parse is shared). metaMS also hard-codes cores=1
    whenever any input is a Thermo .raw file (its own in-process
    multiprocessing.Pool is not safe with the .NET raw-file reader) -
    confirmed by real testing (a 624-file single-invocation run pegged one
    core the whole time). That restriction is per-process, not per-machine,
    so one job per sample gets real OS-level parallelism via Snakemake's own
    scheduler instead, without touching metaMS's internals.

    Isolated working directory: metaMS's prepare_metadata()/process_ms2()
    write/read "metadata.pkl" relative to the CURRENT WORKING DIRECTORY, not
    a configurable path. Parallel jobs launched from a shared cwd would
    silently clobber each other's metadata.pkl; each job here gets its own
    scratch dir via `cd`.

    metaMS names its per-file output directory after the INPUT file's own
    basename stem, not any GLAMR sample ID - confirmed by real testing:
    QC/blank/calibration-standard files reuse the same generic vendor
    filename every acquisition batch (e.g. "blank10", "QC_A_1"), which
    silently overwrote earlier batches' output when they were all run
    together in one invocation (205/624 real samples lost this way in the
    first production run). Fixed by staging the file under its unique
    samp_ID name before invoking metaMS, so its own naming scheme becomes
    collision-proof.
    """
    input:
        file = get_sample_thermo_or_mzml_file,
        msp = "data/reference/metabolomics/gnps_curated.msp",
    output:
        csv = "data/projects/{dataset}/metabolomics/{group_id}/metams/{sample}.corems/{sample}.csv",
    params:
        # All absolute: confirmed by real testing that this rule's own `cd
        # {workdir}` (needed to isolate metaMS's cwd-relative metadata.pkl -
        # see docstring) breaks every relative path passed to metaMS,
        # including its OWN inputs - it checks corems_toml_path first and
        # fails there ("Corems toml file not found") without ever reaching
        # the sample file, but every one of these would have broken the
        # same way. Path(...).resolve() runs at rule-parse time, while cwd
        # is still the repo root, so it's unaffected by the shell's cd.
        corems_params = lambda wc: str(Path(get_group_profile_param(wc.dataset, wc.group_id, 'corems_params')).resolve()),
        scan_translator = lambda wc: str(Path(get_profile_param(wc.dataset, 'scan_translator')).resolve()),
        outdir = lambda wc, output: str(Path(output.csv).resolve().parents[1]),
        workdir = lambda wc: f"/tmp/metams_{wc.dataset}_{wc.group_id}_{wc.sample}",
        msp_abs = lambda wc, input: str(Path(input.msp).resolve()),
        file_abs = lambda wc, input: str(Path(input.file).resolve()),
        # Thermo's .NET raw-file reader (via pythonnet in this container)
        # fails to open the file through a symlink ("Error opening
        # (ThermoFisher.CommonCore.RawFileReader.Facade.RawFileLoader)"),
        # confirmed by real testing - copy instead (~54MB/file for this
        # dataset, well within scratch disk space).
        staged_file = lambda wc, input: f"/tmp/metams_{wc.dataset}_{wc.group_id}_{wc.sample}/{wc.sample}{Path(input.file).suffix}",
        stage_script = str(Path("code/metams_stage_mzml.sh").resolve()),
    threads: 1
    resources: cpus=1, mem_mb=4000, time_min=60
    benchmark: "benchmarks/metab_metams/{dataset}_{group_id}_{sample}.txt"
    log: "logs/metab_metams/{dataset}_{group_id}_{sample}.log"
    container: "docker://microbiomedata/metams:3.5.2"
    shell:
        """
        . code/shell_prelude {log}
        mkdir -p {params.workdir:q} {params.outdir:q}
        cd {params.workdir:q}
        # copies Thermo input; gives other vendors' mzML the Thermo-style ids and
        # filter strings CoreMS expects (see the script)
        bash {params.stage_script:q} {params.file_abs:q} {params.staged_file:q}

        metaMS run-lcms-metabolomics-workflow \
            -i {params.staged_file:q} \
            -o {params.outdir:q} \
            -c {params.corems_params:q} \
            -m {params.msp_abs:q} \
            -s {params.scan_translator:q} \
            -j 1

        cd /
        rm -rf {params.workdir:q}
        """


def get_group_metams_outputs(wc):
    manifest = checkpoints.metab_run_manifest.get(dataset=wc.dataset).output.groups
    samples = metab_get_group_samples(manifest, wc.group_id)
    return [
        f"data/projects/{wc.dataset}/metabolomics/{wc.group_id}/metams/{s}.corems/{s}.csv"
        for s in samples
    ]


rule metab_metams:
    """Aggregates the per-sample metab_metams_one outputs for a group."""
    input:
        get_group_metams_outputs,
    output:
        done = touch("data/projects/{dataset}/metabolomics/{group_id}/metams/done"),


# ---------------------------------------------------------------------------
# Stage 5: QC, filtering, normalization
# ---------------------------------------------------------------------------

rule metab_feature_qc:
    """
    Blank filtering, drift correction, precision (CV) filtering, amount
    normalization - each conditional on the required ms_role rows actually
    existing in the manifest, so a project with no blanks or no QC pools
    still completes (with a report noting what was skipped) rather than
    failing. Not every future project will have set_41's QC design.
    """
    input:
        feature_table = rules.metab_mzmine_batch.output.csv,
        manifest = "data/projects/{dataset}/metabolomics/run_manifest.tsv",
    output:
        filtered = "data/projects/{dataset}/metabolomics/{group_id}/qc/feature_table_filtered.tsv",
        report = "data/projects/{dataset}/metabolomics/{group_id}/qc/qc_report.html",
    params:
        blank_fold_change = config.get('metabolomics', {}).get('blank_fold_change', 5),
        qc_max_cv = config.get('metabolomics', {}).get('qc_max_cv', 0.3),
    threads: 4
    resources: mem_mb=8000, time_min=60
    benchmark: "benchmarks/metab_feature_qc/{dataset}_{group_id}.txt"
    log: "logs/metab_feature_qc/{dataset}_{group_id}.log"
    container: "docker://eandersk/r_microbiome"
    shell:
        """
        . code/shell_prelude {log}
        ./code/metabolomics_feature_qc.R \
            --features {input.feature_table} \
            --manifest {input.manifest} \
            --group {wildcards.group_id} \
            --blank-fold-change {params.blank_fold_change} \
            --qc-max-cv {params.qc_max_cv} \
            --out-table {output.filtered} \
            --out-report {output.report}
        """


# ---------------------------------------------------------------------------
# Stage 5b: CyanopeptideMatching (CPM) - targeted cyanopeptide screen (opt-in)
# ---------------------------------------------------------------------------
# https://github.com/sheffera01/CyanopeptideMatchingPythonVersion: MassQL search
# of MS2 spectra for class-specific diagnostic product ions, precursor merging,
# MS1 XIC areas, and m/z matching to CyanoMetDB. Runs on the same (RT-aligned
# when method "landmark") mzML files as mzmine. Enabled with
# config metabolomics.run_cpm. Uses the glamr-fixes branch of CPM (see
# config/conda_yaml/cpm.yaml), which fixes precursor merging, cross-file scan
# joins, MS1 areas (every file, not only MS2-triggered ones) and charge-aware
# ppm CyanoMetDB matching; see presentations/metabolomics_pipeline/cpm_evaluation.md.

CPM_CLASSES = ['MC', 'MP', 'AB', 'AR', 'MG']


def get_cpm_param(dataset, key, default):
    """ metabolomics.cpm setting, overridden by the dataset's method profile `cpm:` block """
    profile_cpm = config['metabolomics']['method_profiles'][get_method_profile(dataset)].get('cpm', {})
    return profile_cpm.get(key, config.get('metabolomics', {}).get('cpm', {}).get(key, default))


rule metab_cpm_metadata:
    """ CPM metadata CSV (filename, sample_type, batch) from the run manifest """
    input:
        manifest = lambda wc: checkpoints.metab_run_manifest.get(dataset=wc.dataset).output.manifest,
        groups = lambda wc: checkpoints.metab_run_manifest.get(dataset=wc.dataset).output.groups,
    output:
        csv = "data/projects/{dataset}/metabolomics/{group_id}/cpm/metadata.csv",
    resources: mem_mb=500, time_min=10
    log: "logs/metab_cpm_metadata/{dataset}_{group_id}.log"
    run:
        import csv
        rows = _read_group_manifest_rows(input.manifest, input.groups, wildcards.group_id)
        with open(output.csv, 'w', newline='') as ofile:
            w = csv.writer(ofile)
            w.writerow(['filename', 'sample_type', 'batch'])
            for r in rows:
                # CPM only treats sample_type == "blank" specially; keep other roles as-is
                w.writerow([f"{r['sample']}.mzML", r['ms_role'], r['batch_id']])


rule metab_cpm_ms1_points:
    """
    CPM's MS1 point table, built once with CPM's own extractor and shared by
    all class jobs (the CLI would rebuild it per class). intensity_min keeps
    it tractable: unfiltered it is ~170 MB of CSV per set_41 file (~100 GB
    for the dataset, loaded whole into pandas); 1e5 keeps ~9 % of points and
    only clips the low tails of the XICs CPM integrates.
    """
    input:
        mzml = get_group_mzml_files,
    output:
        csv = temp("data/projects/{dataset}/metabolomics/{group_id}/cpm/ms1_points.csv"),
    params:
        outdir = lambda wc, output: str(Path(output.csv).parent / "ms1_points_tmp"),
        intensity_min = lambda wc: get_cpm_param(wc.dataset, 'ms1_intensity_min', 1e5),
        rt_min = lambda wc: get_cpm_param(wc.dataset, 'rt_min', 2.0),
        rt_max = lambda wc: get_cpm_param(wc.dataset, 'rt_max', 25.0),
    resources: mem_mb=60000, time_min=240
    benchmark: "benchmarks/metab_cpm_ms1_points/{dataset}_{group_id}.txt"
    log: "logs/metab_cpm_ms1_points/{dataset}_{group_id}.log"
    conda: "../config/conda_yaml/cpm.yaml"
    shell:
        """
        . code/shell_prelude {log}
        rm -rf {params.outdir:q}
        printf '%s\n' {input.mzml:q} > {params.outdir:q}.list
        python - {params.outdir:q} {params.outdir:q}.list {params.intensity_min} {params.rt_min} {params.rt_max} <<'PY'
import sys
from pathlib import Path
from cpm.cpm import build_ms1_points
outdir, listfile, imin, rtmin, rtmax = sys.argv[1:]
files = Path(listfile).read_text().split()
out = build_ms1_points(files, Path(outdir), mz_round=None, intensity_min=float(imin),
                       assume_time_unit="min", rt_window=(float(rtmin), float(rtmax)))
print(out)
PY
        mv {params.outdir:q}/ms1_points_*.csv {output.csv:q}
        rm -rf {params.outdir:q} {params.outdir:q}.list
        """


rule metab_cpm:
    """
    One CPM class per job (MC microcystins, MP cyanopeptolins, AB
    anabaenopeptins, AR aeruginosins, MG microginins), with CPM's blank filter
    and batch correction. CPM writes a timestamped run directory; it is moved
    to a stable path. Each job gets its own TMPDIR, removed after (upstream
    CPM leaves a temporary copy of every mzML there). CPM is run through code/cpm_run.py, which disables its figures: several
    scale the canvas with the data and, at 624 files, failed with
    std::bad_alloc or grew past 700 GB and were OOM-killed.
    """
    input:
        mzml = get_group_mzml_files,
        metadata = rules.metab_cpm_metadata.output.csv,
        ms1_points = rules.metab_cpm_ms1_points.output.csv,
    output:
        outdir = directory("data/projects/{dataset}/metabolomics/{group_id}/cpm/{cpm_class}"),
    wildcard_constraints:
        cpm_class = "|".join(CPM_CLASSES),
    params:
        workdir = lambda wc, output: output.outdir + "_work",
        tol = lambda wc: get_cpm_param(wc.dataset, 'tol', 0.01),
        tol_da = lambda wc: get_cpm_param(wc.dataset, 'tol_da', 0.05),
        tol_ppm = lambda wc: get_cpm_param(wc.dataset, 'tol_ppm', 5),
        rt_min = lambda wc: get_cpm_param(wc.dataset, 'rt_min', 2.0),
        rt_max = lambda wc: get_cpm_param(wc.dataset, 'rt_max', 25.0),
    threads: 1
    # set_41 (624 files, 12.8 GB MS1 table): ~2-2.3 h per class
    resources: mem_mb=45000, time_min=720
    benchmark: "benchmarks/metab_cpm/{dataset}_{group_id}_{cpm_class}.txt"
    log: "logs/metab_cpm/{dataset}_{group_id}_{cpm_class}.log"
    conda: "../config/conda_yaml/cpm.yaml"
    shell:
        """
        . code/shell_prelude {log}
        rm -rf {params.workdir:q}
        mkdir -p {params.workdir:q}/tmp
        export TMPDIR=$(realpath {params.workdir:q}/tmp)
        python code/cpm_run.py --class-tag {wildcards.cpm_class} \
            --files {input.mzml:q} \
            --metadata {input.metadata:q} \
            --ms1-points-file {input.ms1_points:q} \
            --output-root {params.workdir:q} \
            --tol {params.tol} --tol-da {params.tol_da} --tol-ppm {params.tol_ppm} \
            --rt-min {params.rt_min} --rt-max {params.rt_max} \
            --blank-filter --batch-correct
        rm -rf {params.workdir:q}/tmp
        mv {params.workdir:q}/{wildcards.cpm_class}_run_* {output.outdir:q}
        mv {params.workdir:q}/pipeline_log {output.outdir:q}/ 2>/dev/null || true
        rm -rf {params.workdir:q}
        """


# ---------------------------------------------------------------------------
# Stage 6: top rule
# ---------------------------------------------------------------------------

def get_metabolomics_pipeline_output(wc):
    groups_file = checkpoints.metab_run_manifest.get(dataset=wc.dataset).output.groups
    group_ids = metab_get_group_ids(groups_file)
    files = []
    for g in group_ids:
        files.append(f"data/projects/{wc.dataset}/metabolomics/{g}/qc/feature_table_filtered.tsv")
        files.append(f"data/projects/{wc.dataset}/metabolomics/{g}/ms2query/results.csv")
        files.append(f"data/projects/{wc.dataset}/metabolomics/{g}/metams/done")
        if config.get('metabolomics', {}).get('run_sirius', False):
            files.append(f"data/projects/{wc.dataset}/metabolomics/{g}/sirius/summaries/structure_identifications.tsv")
        # CPM matches [M+H]+ precursors, positive mode only
        if config.get('metabolomics', {}).get('run_cpm', False) and g.endswith('_positive'):
            files.extend(f"data/projects/{wc.dataset}/metabolomics/{g}/cpm/{c}" for c in CPM_CLASSES)
    return files


rule metabolomics_pipeline_dataset:
    """ The top rule for the metabolomics pipeline """
    input: get_metabolomics_pipeline_output
    output: "data/projects/{dataset}/metabolomics_pipeline_done"
    resources: mem_mb=100, time_min=2
    shell:
        """
        for i in {input:q}; do echo "$i" >> {output}; done
        """
