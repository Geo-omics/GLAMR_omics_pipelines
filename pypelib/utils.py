from contextlib import (contextmanager, ExitStack, redirect_stdout,
                        redirect_stderr)
from functools import wraps
from inspect import signature
import json
from pathlib import Path
import sys
import traceback


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
                        raise RuntimeError(f'duplicate file: {value}')
                    rows[value] = row
                else:
                    row[col] = value
    if not rows:
        raise RuntimeError('no data in file')
    return rows


class UsageError(Exception):
    pass


class logme:
    """
    3-in-1 context manager to to write stdout/stderr/tracebacks to a log file.

    For use in with statements in Snakefile rules.  Can also be used as a
    function decorator, which is the preferred usage for rules that make a
    single function call in their run section.

    After the code block or wrapped function is done, the log file will be read
    and stdout and stderr will be also be printed to the regular stdout.  If an
    exception occurred inside the block the tracback will not be printed to
    stdout as snakemake is expected to print it always.

    The rule's log section must define a single log file, which is then either
    passed to logme() in the with statement or is passed via log=log kwarg to
    the decorated function.  The passed log can also be a str or PathLike.  The
    wrapper will consume the log kwarg unless the decocrated function has also
    a log parameter.  This allows the decorated function to write additional
    data to the log file without also printing it to stdout.  The decorated
    function must xpect log to be file-like, open for writing.
    """
    def __init__(self, log=None):
        # log is the rule's log but can also be a str to PathLike
        self.log = log

    def __enter__(self):
        if self.log is None:
            raise ValueError('need path to log file')

        self.logfile = open(str(self.log), 'w+')
        self.redir_out_cm = redirect_stdout(self.logfile)
        self.redir_out_cm.__enter__()
        self.redir_err_cm = redirect_stderr(self.logfile)
        self.redir_err_cm.__enter__()
        return self.logfile

    def __exit__(self, exc_type, exc, exc_tb):
        self.redir_err_cm.__exit__(exc_type, exc, exc_tb)
        self.redir_out_cm.__exit__(exc_type, exc, exc_tb)
        self.logfile.seek(0)
        print(self.logfile.read())  # print all output
        if exc is not None:
            # print traceback to log file
            traceback.print_tb(exc_tb, file=self.logfile)
        self.logfile.close()

    def _recreate_cm(self):
        """ cf. contextlib.ContextDecorator """
        return self

    def __call__(self, func):
        """ cf. contextlib.ContextDecorator """
        @wraps(func)
        def wrapper(*args, log=None, **kwargs):
            keep_logkw = signature(func).parameters.get('log') is not None
            # NOTE: too lazy to check the parameter kind, func() call will fail
            # if log is positional-only

            if log is None:
                # Allow funct to be called without the log kwargs, in which
                # case the function runs as-is
                if keep_logkw:
                    kwargs['log'] = None
                return func(*args, **kwargs)
            else:
                self.log = log
                with self._recreate_cm() as logfile:
                    if keep_logkw:
                        kwargs['log'] = logfile
                    return func(*args, **kwargs)
        return wrapper


@logme()
def test_logme(good='good', bad='bad'):
    print(good)
    print(bad, file=sys.stderr)
    raise RuntimeError('the ugly')
>>>>>>> 258a535 (pypelib/utils: add logging helper for snakemake rules)
