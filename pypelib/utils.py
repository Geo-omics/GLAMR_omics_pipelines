from contextlib import (chdir, contextmanager, ExitStack, redirect_stdout,
                        redirect_stderr)
from functools import total_ordering, wraps
import gzip
from inspect import signature
import json
from itertools import batched
from os import PathLike
from pathlib import Path
import random
import re
import shutil
from subprocess import PIPE, run
import sys
import tarfile
from tempfile import TemporaryDirectory
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


def conda_run(prefix, cmd, *args, **kwargs):
    """
    Wrapper around subprocess.run to run command inside conda environment

    prefix:
        Path to the conda env, the prefix in conda lingo.  If this is None,
        then the command is run as-is.

    cmd:
        The first argument to subprocess.run/Popen, , what the subprocess docs
        call 'args', can be a list of str or a PathLike.

    *args, **kargs:
        Any other positional and keyword arguments to be passed to run()

    If conda is not installed, the command in run as usual.
    """
    if prefix is None or shutil.which('conda') is None:
        return run(cmd, *args, **kwargs)
    else:
        conda_run_cmd = ['conda', 'run', '--prefix', str(prefix)]
        if isinstance(cmd, PathLike):
            cmd = conda_run_cmd.append(cmd)
        else:
            # a list
            cmd = conda_run_cmd + cmd
        return run(cmd, *args, **kwargs)


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
            traceback.print_exception(exc, file=self.logfile)
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


@total_ordering
class PipelineVersion:
    """
    Models pipeline versions

    Assumes we're run out of a git repository and that there is at least one
    annotated tag with a version name, e.g. "v1.2.3" in the history.
    """
    git_describe_long_pat = re.compile(
        r'^v(?P<base>[0-9.]+)-(?P<commit_count>[0-9]+)-g(?P<hash>[0-9a-f]+)'
        r'(-(?P<suffix>\w+))?$'
    )

    def __init__(self, version_string):
        """
        Get an instance from output of 'git describe --long'
        """
        m = self.git_describe_long_pat.match(version_string)
        if m is None:
            # e.g. no version tag was set before and git describe only returns
            # an abbreviated hash
            raise ValueError(f'invalid version string: {version_string}')
        self.base = tuple(int(i) for i in m['base'].split('.'))
        self.count = int(m['commit_count'])
        self.hash = m['hash']
        self.suffix = m['suffix']

    def __repr__(self):
        return f"{type(self).__name__}('{self}')"

    def __str__(self):
        base = '.'.join(map(str, self.base))
        val = f'v{base}-{self.count}-g{self.hash}'
        if self.suffix:
            val += f'-{self.suffix}'
        return val

    def __hash__(self):
        return hash((self.base, self.count, self.hash, self.suffix))

    def __eq__(self, other):
        if other is None:
            return False
        return all((
            self.base == other.base,
            self.count == other.count,
            self.hash == other.hash,
            self.suffix == other.suffix,
        ))

    def __lt__(self, other):
        if self.base == self.base:
            if self.count == other.count:
                if self.hash == other.hash:
                    # non-suffix sorts before suffix, otherwise
                    # suffices are incomparable
                    if self.suffix is None and other.suffix:
                        return True
                    elif self.suffix and other.suffix and self.suffix != other.suffix:  # noqa:E501

                        raise ValueError(
                            'versions differing by suffix are non-comparable'
                        )
                    else:
                        False  # they're equal
                else:
                    raise ValueError(
                        'versions differing only by hash are non-comparable'
                    )
            else:
                return self.count < other.count
        else:
            return self.base < other.base

    LOG_MSG_PREFIX = 'omics pipeline version: '

    @classmethod
    def write_to_log(cls, path):
        """ Append current version message to a file """
        with open(path, 'a') as ofile:
            try:
                pl_version = PipelineVersion.current()
            except Exception as e:
                ofile.write(
                    f'[WARNING] failed getting omics pipeline version: '
                    f'{e.__class__.__name__}: {e}\n'
                )
            ofile.write(f'{cls.LOG_MSG_PREFIX}{pl_version}\n')

    @classmethod
    def current(cls):
        """ get the current version """
        return cls.from_commit(None)

    @classmethod
    def from_commit(cls, commit_hash):
        """
        Get version for given git commit hash

        commit_hash [str]: This can also be None to get the current version.
        """
        cmd = ['git', 'describe', '--always', '--long']
        if commit_hash is None:
            cmd += ['--dirty', '--broken']
        else:
            cmd.append(commit_hash)

        p = run(cmd, check=True, stdout=PIPE)
        version_txt = p.stdout.decode().splitlines()
        if len(version_txt) == 1:
            return cls(version_txt[0].strip())
        else:
            raise RuntimeError(
                f'parsing output of {cmd} failed {version_txt=}'
            )

    @classmethod
    def from_timestamp(cls, timestamp):
        if timestamp is None:
            commit_hash = None
        else:
            commit_hash = cls.get_commit_before(timestamp)
        return cls.from_commit(commit_hash)

    @staticmethod
    def get_commit_before(timestamp):
        """
        Get ID of most recent commit before given datetime

        timestamp:
            ISO-formatted str of timestamp or similar formatting that git
            understands.

        This assumes running from within a git repository.
        """
        cmd = ['git', 'rev-list', '--max-count', '1', '--before', timestamp]
        p = run(cmd, check=True, stdout=PIPE)
        lines = p.stdout.decode().splitlines()
        if len(lines) == 1:
            return lines[0].strip()
        else:
            raise RuntimeError(
                f'expected single line of stdout but got {lines=}'
            )


def make_test_dataset(
    parent,
    num_samples=3,
    num_reads=1500,
    sample_type='amplicons',
    random_seed=67,
):
    """ compile a small test dataset from a given parent data set """
    random.seed(random_seed)
    parent_dir = Path('data/projects/') / parent
    parent_samp_dirs = random.sample(
        sorted((parent_dir / sample_type).glob('samp_*')),
        k=num_samples,
    )
    parent_raws = [
        (i.name, list((i / 'reads').glob('raw_*_reads.fastq.gz')))
        for i in parent_samp_dirs
    ]

    min_num_reads = int(num_reads * 0.8)
    max_num_reads = int(num_reads * 1.2)
    sra_accn_pat = re.compile(rb'@SRR[0-9]+')  # expected start of fastq record

    pref, _, number = parent.partition('_')
    name = f'{pref}_T{number}'  # name of test data set
    with TemporaryDirectory() as tmpd:
        tmpd = Path(tmpd)
        samp_dir = tmpd / 'data' / 'omics' / sample_type
        samp_dir.mkdir(parents=True)
        link_dir = tmpd / 'data' / 'projects' / name / sample_type
        link_dir.mkdir(parents=True)

        for parent_samp, parent_raw_files in parent_raws:
            pref, _, samp_number = parent_samp.partition('_')
            samp_name = f'{pref}_T{samp_number}'
            test_accn = b'@TEST' + str(samp_number).encode()

            reads_dir = samp_dir / samp_name / 'reads'
            reads_dir.mkdir(parents=True)
            link_target = Path('..') / '..' / '..' / 'omics' / sample_type / samp_name  # noqa:E501
            (link_dir / samp_name).symlink_to(link_target)

            read_count = random.randrange(min_num_reads, max_num_reads)
            for p in parent_raw_files:
                test_raw = reads_dir / p.name
                with gzip.open(p) as ifile, gzip.open(test_raw, 'wb') as ofile:
                    # 4 lines per read in fastq file
                    reads = batched(ifile, 4, strict=True)
                    for num, (head, seq, plus, qual) in enumerate(reads):
                        if num >= read_count:
                            break

                        # also replacing these to avoid any confusion w/read
                        # data
                        head = sra_accn_pat.sub(test_accn, head)
                        ofile.write(head)
                        ofile.write(seq)
                        ofile.write(plus)
                        ofile.write(qual)

        def reset(tarinfo):
            tarinfo.uid = tarinfo.gid = 0
            tarinfo.uname = tarinfo.gname = 'root'
            return tarinfo

        tarpath = Path(f'{name}.tar.gz').absolute()
        with chdir(tmpd):
            print(f'Writing {tarpath.name} ...', end='', flush=True)
            with tarfile.open(tarpath, 'w:gz') as tarf:
                tarf.add('data', filter=reset)
    print('[OK]')
    print('To unpack test data in the real data location, run:')
    if Path('data').is_symlink():
        # Have to unpack at the real place as tar -xf would overwrite the
        # symlink, creating completely new "data" directory
        print('  $ cd', Path('data').resolve().parent)
    print('  $ tar -xf', tarpath)
