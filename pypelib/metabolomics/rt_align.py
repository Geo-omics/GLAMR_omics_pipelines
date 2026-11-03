"""
Landmark-based retention-time alignment of mzML files, applied before MZmine.

Why this exists instead of MZmine's own scan RT correction: on set_41 the
between-batch RT shift is a hump-shaped function of RT (raw apex RTs,
2024-06-20 vs 06-25 batches: ~0 at 2-5 min, +0.28 min at 7-8 min, ~0 again at
12-13 min). MZmine 4.10.6's only correction function (MultiLinearCorrection)
smooths each file's landmark offsets and then fits an asymmetric least-squares
baseline (asymmetry 0.01, i.e. the lower envelope) to them; for a hump that
over-corrects the flat regions (-0.14 min at 2-5 min on set_41), which created
new split features there. This module does the same job the standard way
(XCMS "peak groups" style): measure landmark ions in every run, fit a robust
LOESS of RT offset vs RT per run, and rewrite the scan start times.

Steps (all numpy, no scipy - runs in the pipeline's snakemake environment):

1. build_landmarks(): pick candidate ions from the reference runs (intense,
   chromatographically unambiguous), measure their apex RT in every
   reference run, keep those found in most references; the landmark's
   reference RT is the median apex RT across reference runs.
2. fit_run(): measure the landmarks in one run, fit offset(RT) with a robust
   local-linear LOESS, and return a monotonic mapping original -> corrected RT.
3. write_aligned_mzml(): copy the mzML with corrected "scan start time"
   values, rebuilding the byte-offset index and SHA-1 checksum.
"""
import argparse
import base64
import csv
import json
import re
import xml.etree.ElementTree as ET
import zlib
from pathlib import Path

import numpy as np

from ..utils import logme
from .mzml_rewrite import rewrite_indexed_mzml

NS = '{http://psi.hupo.org/ms/mzml}'

DEFAULT_PPM = 5.0
DEFAULT_MIN_HEIGHT = 2e6       # landmark apex intensity (counts)
DEFAULT_WINDOW = 0.8           # min, search window around a landmark's reference RT
DEFAULT_SPAN = 0.25            # LOESS span (fraction of a run's landmarks)
DEFAULT_MIN_REF_FRACTION = 0.8
DEFAULT_MAX_CANDIDATES = 3000


# ---------------------------------------------------------------------------
# Reading MS1 data
# ---------------------------------------------------------------------------

def _decode(bda):
    params = {c.get('accession') for c in bda.iter(NS + 'cvParam')}
    raw = base64.b64decode(bda.find(NS + 'binary').text or '')
    if 'MS:1000574' in params:
        raw = zlib.decompress(raw)
    arr = np.frombuffer(raw, dtype=np.float64 if 'MS:1000523' in params else np.float32)
    kind = 'mz' if 'MS:1000514' in params else 'int' if 'MS:1000515' in params else None
    return kind, arr


def load_ms1(path):
    """ Returns (rts [min], list of m/z arrays, list of intensity arrays) for MS1 scans """
    rts, mzs, ints = [], [], []
    for _, el in ET.iterparse(str(path), events=('end',)):
        if el.tag != NS + 'spectrum':
            continue
        cv = {c.get('accession'): c for c in el.iter(NS + 'cvParam')}
        if cv.get('MS:1000511') is not None and cv['MS:1000511'].get('value') == '1':
            rt = float(cv['MS:1000016'].get('value'))
            if cv['MS:1000016'].get('unitName') == 'second':
                rt /= 60
            arrays = dict(_decode(b) for b in el.iter(NS + 'binaryDataArray'))
            # empty scans (seen in Agilent data via msconvert) carry no binary arrays at all
            arrays.setdefault('mz', np.empty(0)); arrays.setdefault('int', np.empty(0))
            order = np.argsort(arrays['mz'])
            rts.append(rt)
            mzs.append(arrays['mz'][order].astype(np.float64))
            ints.append(arrays['int'][order].astype(np.float64))
        el.clear()
    return np.array(rts), mzs, ints


def xic_matrix(run, target_mz, ppm=DEFAULT_PPM):
    """ (n_targets x n_scans) summed intensity within +-ppm of each target m/z """
    rts, mzs, ints = run
    target_mz = np.asarray(target_mz, dtype=float)
    lo_mz, hi_mz = target_mz * (1 - ppm * 1e-6), target_mz * (1 + ppm * 1e-6)
    out = np.zeros((len(target_mz), len(rts)))
    for j, (mz, it) in enumerate(zip(mzs, ints)):
        cs = np.concatenate([[0.0], np.cumsum(it)])
        out[:, j] = cs[np.searchsorted(mz, hi_mz, 'right')] - cs[np.searchsorted(mz, lo_mz, 'left')]
    return out


def apex_rts(run, target_mz, expected_rt, window=DEFAULT_WINDOW, min_height=DEFAULT_MIN_HEIGHT,
             ppm=DEFAULT_PPM, max_other=0.3):
    """
    Apex RT of each target within expected_rt +- window, or NaN when the
    peak is too weak or ambiguous (another local signal > max_other x apex
    more than 0.1 min away). Apex refined by a 3-point parabola.
    """
    rts = run[0]
    x = xic_matrix(run, target_mz, ppm)
    expected_rt = np.asarray(expected_rt, dtype=float)
    out = np.full(len(target_mz), np.nan)
    for i in range(len(target_mz)):
        sel = np.where(np.abs(rts - expected_rt[i]) <= window)[0]
        if len(sel) < 5:
            continue
        y = x[i, sel]
        j = int(np.argmax(y))
        if y[j] < min_height or j == 0 or j == len(y) - 1:
            continue
        far = np.abs(rts[sel] - rts[sel][j]) > 0.1
        if far.any() and y[far].max() > max_other * y[j]:
            continue
        y0, y1, y2 = y[j - 1], y[j], y[j + 1]
        denom = y0 - 2 * y1 + y2
        shift = 0.5 * (y0 - y2) / denom if denom < 0 else 0.0
        step = (rts[sel][j + 1] - rts[sel][j - 1]) / 2
        out[i] = rts[sel][j] + np.clip(shift, -0.5, 0.5) * step
    return out


# ---------------------------------------------------------------------------
# Landmarks
# ---------------------------------------------------------------------------

def candidate_ions(run, min_height=DEFAULT_MIN_HEIGHT, ppm=DEFAULT_PPM, max_n=DEFAULT_MAX_CANDIDATES):
    """ (mz, apex_rt) of intense, unambiguous ions in one run """
    rts, mzs, ints = run
    peaks = np.concatenate([mz[it >= min_height] for mz, it in zip(mzs, ints)])
    if len(peaks) == 0:
        return np.empty(0), np.empty(0)
    peaks.sort()
    # cluster consecutive m/z values closer than ppm into one candidate ion
    breaks = np.where(np.diff(peaks) > peaks[:-1] * ppm * 1e-6)[0] + 1
    groups = np.split(peaks, breaks)
    cand_mz = np.array([np.median(g) for g in groups if len(g) >= 3])
    x = xic_matrix(run, cand_mz, ppm)
    best = x.max(axis=1)
    keep = np.argsort(best)[::-1][:max_n]
    cand_mz, apex_guess = cand_mz[keep], rts[np.argmax(x[keep], axis=1)]
    apex = apex_rts(run, cand_mz, apex_guess, window=1.0, min_height=min_height, ppm=ppm)
    ok = np.isfinite(apex)
    return cand_mz[ok], apex[ok]


@logme()
def build_landmarks(reference_runs, min_height=DEFAULT_MIN_HEIGHT, ppm=DEFAULT_PPM,
                    window=DEFAULT_WINDOW, min_fraction=DEFAULT_MIN_REF_FRACTION, log=None):
    """
    reference_runs: dict name -> load_ms1() result.
    Candidates come from every reference run (so batch-specific RTs are all
    represented), merged by m/z/RT, then measured in every reference run.
    Returns a list of dicts: mz, ref_rt, n_found.
    """
    names = list(reference_runs)
    cands = []
    for n in names:
        mz, rt = candidate_ions(reference_runs[n], min_height, ppm)
        cands.extend(zip(mz, rt))
        print(f'{n}: {len(mz)} candidate ions')
    # merge candidates seen in several runs (same m/z, RT within window)
    cands.sort()
    merged = []
    for mz, rt in cands:
        if merged and abs(mz - merged[-1][0]) <= mz * ppm * 1e-6 and abs(rt - merged[-1][1]) <= window:
            merged[-1][2].append(rt)
        else:
            merged.append([mz, rt, [rt]])
    if not merged:
        print(f'no candidate ions >= {min_height:g} in {len(names)} reference runs: no landmarks, '
              f'runs will be left unaligned')
        return []
    lm_mz = np.array([m[0] for m in merged])
    lm_rt = np.array([np.median(m[2]) for m in merged])
    found = np.array([apex_rts(reference_runs[n], lm_mz, lm_rt, window, min_height / 4, ppm) for n in names])
    frac = np.isfinite(found).mean(axis=0)
    keep = frac >= min_fraction
    ref_rt = np.nanmedian(found[:, keep], axis=0)
    # drop landmarks whose reference RTs collide with another landmark of the same m/z
    out = [dict(mz=float(m), ref_rt=float(r), n_found=int(np.isfinite(found[:, k]).sum()))
           for m, r, k in zip(lm_mz[keep], ref_rt, np.where(keep)[0])]
    print(f'{len(merged)} merged candidates, {len(out)} landmarks found in >= {min_fraction:.0%} '
          f'of {len(names)} reference runs')
    return out


# ---------------------------------------------------------------------------
# Per-run correction
# ---------------------------------------------------------------------------

def robust_loess(x, y, x_eval, span=DEFAULT_SPAN, iterations=3):
    """ Local-linear LOESS with tricube weights and bisquare robustness iterations """
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    k = max(5, int(np.ceil(span * n)))
    robust_w = np.ones(n)

    def fit_at(x0, w_extra):
        d = np.abs(x - x0)
        h = np.partition(d, k - 1)[k - 1] or 1e-9
        w = np.clip(1 - (d / h) ** 3, 0, None) ** 3 * w_extra
        sw = w.sum()
        if sw <= 0:
            return np.nan
        xm, ym = (w * x).sum() / sw, (w * y).sum() / sw
        sxx = (w * (x - xm) ** 2).sum()
        b = (w * (x - xm) * (y - ym)).sum() / sxx if sxx > 0 else 0.0
        return ym + b * (x0 - xm)

    for _ in range(iterations):
        fitted = np.array([fit_at(xi, robust_w) for xi in x])
        resid = y - fitted
        s = np.median(np.abs(resid)) or 1e-9
        robust_w = np.clip(1 - (resid / (6 * s)) ** 2, 0, None) ** 2
    return np.array([fit_at(xe, robust_w) for xe in np.asarray(x_eval, float)])


def fit_run(run, landmarks, window=DEFAULT_WINDOW, min_height=DEFAULT_MIN_HEIGHT, ppm=DEFAULT_PPM,
            span=DEFAULT_SPAN, min_landmarks=20):
    """
    Returns (grid_rt, corrected_grid_rt, stats). The mapping is monotonic and
    is evaluated on a 0.01-min grid over the run; offsets outside the
    landmark RT range are held at the nearest fitted value.
    """
    lm_mz = np.array([l['mz'] for l in landmarks])
    lm_ref = np.array([l['ref_rt'] for l in landmarks])
    rts = run[0]
    if not landmarks:
        grid = np.arange(0, np.ceil(rts.max() * 100) / 100 + 0.01, 0.01)
        return grid, grid.copy(), dict(n_landmarks=0, applied=False)
    obs = apex_rts(run, lm_mz, lm_ref, window, min_height / 4, ppm)
    ok = np.isfinite(obs)
    grid = np.arange(0, np.ceil(rts.max() * 100) / 100 + 0.01, 0.01)
    if ok.sum() < min_landmarks:
        return grid, grid.copy(), dict(n_landmarks=int(ok.sum()), applied=False)
    x, d = obs[ok], obs[ok] - lm_ref[ok]
    inside = (grid >= x.min()) & (grid <= x.max())
    offset = np.empty_like(grid)
    offset[inside] = robust_loess(x, d, grid[inside], span)
    offset[grid < x.min()] = offset[inside][0]
    offset[grid > x.max()] = offset[inside][-1]
    corrected = grid - offset
    corrected = np.maximum.accumulate(corrected)
    corrected = np.maximum(corrected, 0) + np.arange(len(corrected)) * 1e-7   # strictly increasing
    resid = d - np.interp(x, grid, offset)
    stats = dict(n_landmarks=int(ok.sum()), applied=True,
                 median_abs_offset=float(np.median(np.abs(d))),
                 median_abs_residual=float(np.median(np.abs(resid))),
                 max_abs_fitted_offset=float(np.abs(offset).max()))
    return grid, corrected, stats


# ---------------------------------------------------------------------------
# Writing the aligned mzML
# ---------------------------------------------------------------------------

_RT_RE = re.compile(r'(accession="MS:1000016"[^>]*?value=")([^"]+)(")')
_UNIT_RE = re.compile(r'unitName="(second|minute)"')


def write_aligned_mzml(src, dst, grid, corrected):
    """
    Copy src -> dst replacing every "scan start time" with the corrected RT
    (same unit as the input), then rebuild the indexedmzML spectrum /
    chromatogram offsets and the SHA-1 file checksum.
    """
    def fix(line):
        m = _RT_RE.search(line)
        if not m:
            return line
        unit = _UNIT_RE.search(line)
        to_min = 1 / 60 if unit and unit.group(1) == 'second' else 1.0
        new = float(np.interp(float(m.group(2)) * to_min, grid, corrected)) / to_min
        return line[:m.start(2)] + f'{new:.7f}' + line[m.end(2):]

    rewrite_indexed_mzml(src, dst, fix)


# ---------------------------------------------------------------------------
# Snakemake entry points
# ---------------------------------------------------------------------------

@logme()
def landmarks_main(reference_mzmls, output, min_height=DEFAULT_MIN_HEIGHT, ppm=DEFAULT_PPM, log=None):
    runs = {Path(p).stem: load_ms1(p) for p in reference_mzmls}
    lms = build_landmarks(runs, min_height=min_height, ppm=ppm)
    with open(output, 'w', newline='') as ofile:
        w = csv.DictWriter(ofile, fieldnames=['mz', 'ref_rt', 'n_found'], delimiter='\t')
        w.writeheader()
        w.writerows(lms)


@logme()
def align_main(mzml, landmarks_tsv, output_mzml, output_json, min_height=DEFAULT_MIN_HEIGHT,
               ppm=DEFAULT_PPM, log=None):
    with open(landmarks_tsv) as ifile:
        lms = [dict(mz=float(r['mz']), ref_rt=float(r['ref_rt']))
               for r in csv.DictReader(ifile, delimiter='\t')]
    run = load_ms1(mzml)
    grid, corrected, stats = fit_run(run, lms, min_height=min_height, ppm=ppm)
    print(json.dumps(stats))
    write_aligned_mzml(mzml, output_mzml, grid, corrected)
    sample_grid = np.arange(0, grid.max(), 0.25)
    stats['correction_min'] = {f'{g:.2f}': round(float(np.interp(g, grid, corrected) - g), 4)
                               for g in sample_grid}
    with open(output_json, 'w') as ofile:
        json.dump(stats, ofile, indent=2)


def cli():
    argp = argparse.ArgumentParser(description=__doc__)
    sub = argp.add_subparsers(dest='cmd', required=True)
    a = sub.add_parser('landmarks')
    a.add_argument('output')
    a.add_argument('reference_mzmls', nargs='+')
    b = sub.add_parser('align')
    b.add_argument('mzml')
    b.add_argument('landmarks_tsv')
    b.add_argument('output_mzml')
    b.add_argument('output_json')
    args = argp.parse_args()
    if args.cmd == 'landmarks':
        landmarks_main(args.reference_mzmls, args.output)
    else:
        align_main(args.mzml, args.landmarks_tsv, args.output_mzml, args.output_json)


if __name__ == '__main__':
    cli()
