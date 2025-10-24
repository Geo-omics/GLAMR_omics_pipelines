def load_stats(path):
    """
    Load and parse data from a "seqkit stats -T" created reads statistics file

    Returns a dict of dicts.  Top-level keys are the reads fastq file names.
    """
    TEXT_COLS = ['file', 'format', 'type']

    rows = {}
    with open(path) as ifile:
        line = ifile.readline().rstrip('\n')
        header = line.split('\t')
        for line in ifile:
            row0 = line.rstrip('\n').split('\t')
            row = {}
            for col, value in zip(header, row0, strict=True):
                if col not in TEXT_COLS:
                    try:
                        value = int(value)
                    except ValueError as e1:
                        try:
                            value = float(value)
                        except ValueError as e2:
                            raise ValueError(
                                f'Expected a numerical value: {e1} / {e2}'
                            ) from e2

                if col == 'file':
                    if value in rows:
                        raise RuntimeError('duplicate file: {value}')
                    rows[value] = row
                else:
                    row[col] = value
    if not rows:
        raise RuntimeError('no data in file')
    return rows
