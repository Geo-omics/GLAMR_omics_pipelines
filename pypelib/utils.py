from contextlib import contextmanager, ExitStack
import json


@contextmanager
def save_error_file(path):
    """
    Run a block of code and save exception info (the args) to a file.

    Will save the args of an exception to a file (json formatted) and re-raise
    the exception but with all but the first arg (usually the message) dropped.
    Use this if you want to record certain errors for later inspection.

    If path is None, then this is a no-op.
    """
    try:
        yield
    except Exception as exc:
        if path is None:
            raise
        data = {
            'Exception': exc.__class__.__name__,
            'args': exc.args,
        }
        try:
            with open(path, 'w') as ofile:
                json.dump(data, ofile, indent=4)
            print(f'Error info saved as: {path}')
        except Exception as e:
            print(f'failed saving error: {e}')

        if len(exc.args) > 1:
            exc.args = (exc.args[0], f'{len(exc.args) - 1} args suppressed')
        raise  # raise original exception as usual


def load_stats(file_or_path):
    """
    Load and parse data from a "seqkit stats -T" created reads statistics file

    Returns a dict of dicts.  Top-level keys are the reads fastq file names.
    """
    TEXT_COLS = ['file', 'format', 'type']

    rows = {}
    with ExitStack() as estack:
        try:
            ifile = open(file_or_path)
        except TypeError:
            # assume we're given an open filehandle
            ifile = file_or_path
        else:
            estack.enter_context(ifile)

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


class UsageError(Exception):
    pass
