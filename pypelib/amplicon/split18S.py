"""
Dispatch the 18S reads of primer pairs that amplify both 16S and 18S

Some primer pairs, e.g. 515F-Y/926R (Parada et al. 2016), amplify 18S along
with 16S.  Such samples get assigned to their 16S target as usual.  Their reads
are additionally split into 16S and 18S with bbsplit (rule
amplicon_split_16S_18S) and the 18S reads are processed as a separate dada2
target on the 18S model.  This module makes the dispatch files (target
assignments, sample info, and target info) for that 18S target, analogous to
those made by the dispatch module for the regular targets.
"""
import csv
from pathlib import Path
from statistics import median

from .dispatch import SKIP, get_assignments, load_target_table
from .hmm import HMM


def get_split_samples(assignments_file, split_targets):
    """ Get IDs of the samples assigned to one of the given (16S) targets """
    return sorted(
        sample for sample, target
        in get_assignments(path=assignments_file).items()
        if target in split_targets
    )


def read_pair_count(refstats_file, name='18S'):
    """ Get number of read pairs assigned to given reference from bbsplit refstats """
    with open(refstats_file) as ifile:
        for row in csv.DictReader(ifile, delimiter='\t'):
            if row['#name'] == name:
                return int(row['assignedReads']) // 2
    return 0


def make(
    assignments_file,
    target_tab,
    split_targets,
    target_18S,
    reads_files,
    refstats_files,
    min_pairs,
    out_assignments,
    out_samples,
    out_target_tab,
):
    """
    Write the dispatch files for the 18S target

    reads_files:
        dict sample -> (fwd fastq, rev fastq) of the split 18S reads, paths
        relative to the GLAMR root
    refstats_files:
        dict sample -> bbsplit refstats file
    min_pairs:
        Samples with fewer 18S read pairs are set to SKIP.
    """
    hmm, fwd_primer, rev_primer = HMM.parse_target(target_18S)
    samples = get_split_samples(assignments_file, split_targets)
    target_info = {row['sample']: row for row in load_target_table(target_tab)}

    # read extents (on the 16S model) per sample, samples with errors in target
    # detection may lack them, those get the median
    def extent(row, direc):
        try:
            return int(row[f'{direc}_hmmto']) - int(row[f'{direc}_hmmfrom'])
        except (KeyError, ValueError):
            return None
    extents = {}
    for direc in ['fwd', 'rev']:
        known = [e for i in samples if (e := extent(target_info[i], direc)) is not None]
        default = round(median(known))
        extents[direc] = {
            i: e if (e := extent(target_info[i], direc)) is not None else default
            for i in samples
        }

    assignments = []
    sample_info = []
    target_rows = []
    for sample in samples:
        num_pairs = read_pair_count(refstats_files[sample])
        override = '' if num_pairs >= min_pairs else SKIP
        assignments.append({'sample': sample, 'target': target_18S, 'override': override})
        print(f'{sample}: {num_pairs} 18S read pairs{" -> SKIP" if override else ""}')

        fwd_fastq, rev_fastq = reads_files[sample]
        sample_info.append({
            'sample': sample,
            'sample_dir': Path(fwd_fastq).parent,
            'fwd_fastq': fwd_fastq,
            'rev_fastq': rev_fastq,
        })

        # The 18S reads come from the already primer-trimmed and (if needed)
        # swapped 16S target reads.  Coordinates are on the 18S model, taking
        # the read extents from the 16S hits (dada2 only uses them for the
        # read lengths when concatenating).
        row = target_info[sample]
        fwd_len = extents['fwd'][sample]
        rev_len = extents['rev'][sample]
        start, end = hmm.get_target_interval(fwd_primer, rev_primer)
        target_rows.append({
            'sample': sample,
            'layout': row['layout'],
            'model_name': hmm.name,
            'swapped': False,
            'fwd_primer': fwd_primer.name,
            'rev_primer': rev_primer.name,
            'fwd_clean': True,
            'rev_clean': True,
            'fwd_rev_clean': True,
            'rev_fwd_clean': True,
            'fwd_hmmfrom': start,
            'fwd_hmmto': start + fwd_len,
            'rev_hmmfrom': end - rev_len,
            'rev_hmmto': end,
            'overlap_start': end - rev_len,
            'overlap_end': start + fwd_len,
            'errors': '',
        })

    write_tsv(assignments, out_assignments)
    write_tsv(sample_info, out_samples)
    write_tsv(target_rows, out_target_tab)


def write_tsv(rows, output_file):
    print(f'Writing {output_file} ... ', end='', flush=True)
    with open(output_file, 'w') as ofile:
        cols = list(rows[0].keys())
        ofile.write('\t'.join(cols) + '\n')
        for row in rows:
            ofile.write('\t'.join(str(row[i]) for i in cols) + '\n')
    print('[Done]')
