import bisect
from contextlib import (contextmanager, ExitStack, redirect_stdout,
                        redirect_stderr)
from datetime import datetime, UTC
from functools import wraps
from inspect import signature
from io import StringIO
import json
from pathlib import Path
import re
import subprocess
import sys
import traceback

from snakemake.deployment.conda import Env
import yaml


DEFAULT_OVERRIDE_FILENAME = 'override'


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


def override(path, **kwargs):
    for key in kwargs:
        if not key.isidentifier():
            raise ValueError('kwargs keys must be identifiers: {key}')

    path = Path(path)
    if path.is_dir():
        path = path / DEFAULT_OVERRIDE_FILENAME
    if path.is_file():
        with open(path) as ifile:
            lines = ifile.readlines()
    else:
        lines = []

    changed = False
    for lnum, line in enumerate(lines, start=1):
        if not line.strip() or line.strip().startswith('#'):
            continue

        key0, eq, old_value = line.partition('=')
        if not eq:
            raise RuntimeError(
                f'not a key=value assignment on line {lnum}: {path}'
            )
        key = key0.strip()
        if not key.isidentifier():
            raise RuntimeError(
                f'key must be an identifier on line {lnum}: {path}'
            )
        if key in kwargs:
            value = str(kwargs.pop(key)).strip()
            if value is None:
                value = ''
            if value == old_value.strip():
                # keep old line as-is
                continue
            # add space if needed
            spacer = '' if key0[-1].isspace() else ' '
            lines.pop(lnum - 1)
            lines.insert(lnum - 1, f'{key0}{spacer}= {value}\n')
            changed = True

    # add new lines
    for key, value in kwargs.items():
        value = str(value).strip()
        lines.append(f'{key} = {value}\n')
        changed = True

    if changed:
        with open(path, 'w') as ofile:
            ofile.writelines(lines)


def get_override(file, key):
    """
    Return an override

    If the override file exists and the key is overridden this returns a
    non-empty string, otherwise None is returned.
    """
    try:
        with open(file) as ifile:
            for line in ifile:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                key0, _, value = line.partition('=')
                if key == key0.strip():
                    if value := value.strip():
                        return value
    except FileNotFoundError:
        return None


_deps = {}


def dependencies(env_dir):
    pkgname_pat = re.compile(r'^([\w0-9]+)(.*)')
    deps = set()
    for i in Path(env_dir).glob('*.yaml'):
        with open(i) as ifile:
            spec = yaml.load(ifile, yaml.CLoader)
            for j in spec['dependencies']:
                m = pkgname_pat.match(j)
                pkg, version = m.groups()
                deps.add(pkg)
    return deps


class JobIdReused(Exception):
    pass


def get_info_from_log(log):
    """
    Collect certain information from a snakemake run's log file

    Returns list of output files of completed jobs and set of names of rules
    that were run to completion at least once.
    """
    rule_keywords = ['checkpoint', 'rule', 'localrule', 'localcheckpoint']
    rule_pat = re.compile(
        r'^(' + r'|'.join(rule_keywords) + r') (?P<rule_name>\w+):\n'
        r'(?P<keyvals>(\s+\w+: .*\n)+)',
        re.MULTILINE
    )
    finished_pat = re.compile(
        r'^Finished jobid: (?P<jobid>[0-9]+)',
        re.MULTILINE,
    )

    linepos = [0]  # index to get line numbers for error messages
    with open(log) as ifile:
        log_txt = ifile.read()
        ifile.seek(0)
        for line in ifile:
            linepos.append(linepos[-1] + len(line))

    def linenum(pos_or_match):
        if isinstance(pos_or_match, int):
            pos = pos_or_match
        else:
            # assume re.Match object, use start position
            pos = pos_or_match.start()
        return bisect.bisect_right(linepos, pos)

    no_output_jobs = set()
    jobs = {}
    for m in rule_pat.finditer(log_txt):
        job = {}
        job['name'] = m['rule_name']
        jobid = None

        for line in m['keyvals'].splitlines():
            key, _, value = line.strip().partition(': ')
            match key:
                case 'output':
                    if 'output' in job:
                        raise RuntimeError('duplicate output line')
                    job['output'] = [
                        i.removesuffix(' (update)')
                        for i in value.split(', ')
                    ]
                case 'jobid':
                    if jobid is not None:
                        raise RuntimeError('duplicate jobid in job block')
                    jobid = int(value)
                case _:
                    pass

        if jobid is None:
            raise RuntimeError(
                f'jobid missing near line {linenum(m)} {m=}'
            )

        if jobid in jobs or jobid in no_output_jobs:
            raise JobIdReused(
                f'duplicate jobid {jobid} near line {linenum(m)}'
            )

        if 'output' not in job:
            # e.g. target rule
            no_output_jobs.add(jobid)
            continue

        job['lines'] = (linenum(m.start()), linenum(m.end()))  # for debugging
        jobs[jobid] = job

    output = []
    rules = set()
    for m in finished_pat.finditer(log_txt):
        jobid = int(m['jobid'])
        if jobid in no_output_jobs:
            continue
        try:
            job = jobs[jobid]
        except KeyError:
            print(f'for log file {log}')
            print(f'unknown jobid {jobid} near line {linenum(m)}')
            raise
        output += jobs[jobid]['output']
        rules.add(jobs[jobid]['name'])

    return output, rules


def read_omics_checkout(checkout_file):
    """
    Helper to read and parse the checkout file and return dict mapping files to
    most recent mtime.
    """
    MIN_TIME = datetime.min.replace(tzinfo=UTC)
    most_recent = {}
    with open(checkout_file) as ifile:
        for line in ifile:
            mtime, _, path = line.rstrip('\n').partition('\t')
            mtime = datetime.fromisoformat(mtime)
            if most_recent.get(path, MIN_TIME) < mtime:
                most_recent[path] = mtime
    return most_recent


def update_omics_checkout(
    output_files,
    data_root,
    checkout_file=None,
    dry_run=False,
):
    """
    Update checkout file by appending one line for each new or newer file.

    Returns dict of statistics.
    """
    # For historical reasons the checkout file lists files relative to
    # data/omics (the base), but we'll have to allow files under at least
    # data/projects as well, those paths will start with ../projects instead.
    # Below, we resolve() to prevent path traversal escape from data root with
    # Snakemake output files.
    data_root = Path(data_root).resolve(strict=True)
    base = data_root / 'omics'

    old_data = {}
    if checkout_file:
        checkout_file = Path(checkout_file)
        if checkout_file.is_file():
            old_data = read_omics_checkout(checkout_file)

    missing = no_change = ignored = errors = 0
    new_data = []
    for file00 in sorted(output_files):
        file0 = Path(file00)
        if file0.is_absolute() or file0.parts[0] != 'data':
            # TODO: some files get saved outside of data_root, unsure how to
            # handle these
            ignored += 1
            continue

        try:
            file = (data_root.parent / file0).resolve(strict=True)
        except (FileNotFoundError, NotADirectoryError):
            # e.g. output file marked temp()
            # or replay of old log file, data deleted long time ago
            missing += 1
            continue
        except PermissionError:
            # TODO: needs reporting
            errors += 1
            continue

        if not file.is_relative_to(data_root):
            # TODO: some files get saved outside of data_root, unsure how to
            # handle these
            ignored += 1
            continue

        mtime = file.stat().st_mtime
        mtime = datetime.fromtimestamp(mtime).astimezone()
        relpath = Path(file).relative_to(base, walk_up=True)

        if old_mtime := old_data.get(relpath):
            if mtime == old_mtime:
                # no change, nothing to do
                no_change += 1
                continue
            elif mtime < old_mtime:
                # ? maybe data loss
                raise RuntimeError('old time more recent than current')
            else:
                # file got updated, append below
                pass

        new_data.append((mtime, relpath))

    with ExitStack() as estack:
        if checkout_file is None or dry_run:
            ofile = sys.stdout
        else:
            ofile = estack.enter_context(open(checkout_file, 'a'))

        for mtime, relpath in new_data:
            ofile.write(f'{mtime}\t{relpath}\n')

    return {
        'new': len(new_data),
        'missing': missing,
        'no_change': no_change,
        'ignored': ignored,
        'errors': errors,
    }


def get_conda_env_package_versions(env):
    """
    Get installed versions of important packages in given conda environment

    env:
        A snakemake.deployment.conda.Env object.

    Returns dict of package names mapping to version.
    """
    pkgname_pat = re.compile(r'^([\w0-9]+)(.*)')
    packages = []

    with StringIO(env.content.decode()) as ifile:
        spec = yaml.safe_load(ifile)

    for j in spec['dependencies']:
        m = pkgname_pat.match(j)
        pkg, version = m.groups()
        packages.append(pkg)

    if not packages:
        raise RuntimeError('parsing env yaml failed, no packaged')
    regex = '(' + '|'.join(packages) + ')'
    cmd = ['conda', 'list', '--prefix', env.address, '--export', regex]
    p = subprocess.run(cmd, check=True, stdout=subprocess.PIPE)
    out = p.stdout.decode()
    versions = {}
    for line in out.splitlines():
        line = line.strip()
        if line.startswith('#'):
            continue
        pkg_name, version, *_ = line.split('=')
        if not pkg_name or not version:
            raise RuntimeError(f'failed parsing export list: {line}')
        versions[pkg_name] = version

    return versions


def update_versions_file(envs, versions_file=None, dry_run=False):
    # get git description
    cmd = ['git', 'describe', '--dirty', '--broken', '--always']
    p = subprocess.run(cmd, check=True, stdout=subprocess.PIPE)
    git_descr = p.stdout.decode().splitlines()
    if len(git_descr) == 1:
        git_descr = git_descr[0]
    else:
        raise RuntimeError(f'parsing output of {cmd} failed {git_descr=}')

    # get actually installed package versions
    new_data = {}
    for env in envs:
        env_name = Path(env.file.path).stem
        for pkg, version in get_conda_env_package_versions(env).items():
            if pkg not in new_data:
                new_data[pkg] = []
            new_data[pkg].append((env_name, version))

    if versions_file and versions_file.is_file():
        with open(versions_file) as ifile:
            all_data = json.load(ifile)
        data = all_data.get(git_descr, {})  # data for current git checkout
        if not data:
            print(f'new git name: {git_descr}')
    else:
        all_data = {}
        data = {}

    # update existing package versions, deduplicate as needed
    for pkg in data:
        versions = data[pkg] + new_data.pop(pkg, [])
        versions = sorted(set(versions))
        if len(version) == 1:
            data[pkg] = versions[0]
        else:
            data[pkg] = versions

    # add previously unseen packages
    for pkg, versions in new_data.items():
        data[pkg] = versions[0] if len(versions) == 1 else versions

    all_data[git_descr] = data
    output = json.dumps(all_data, indent=4)

    if versions_file and not dry_run:
        print(f'Writing {versions_file} ...', end='', flush=True)
        with open(versions_file, 'w') as ofile:
            ofile.write(output)
            ofile.write('\n')
        print('[OK]')
    else:
        # for testing
        print('[dry_run or no versions file configured]')
        print(output)


def post_production(log, config, rules=None, data_root=None, dry_run=False):
    """
    Post-snakemake run processing

      * Update the checkout file
      *

    Intended to be called via onsuccess and onerror.

    log:
        str or PathLike to snakemake log file.
    config:
        Snakefile's config.
    rules:
        The snakemake rules variable, passed from the Snakefile.

    checkout_file:
        str or PathLike to two-column tab-separated text file listing output
        files and their modification times.  This can be None if no such file
        is configured.
    data_root:
        The data directory.  If None, this will be derived from path to log
        file.
    dry_run [bool]:
        Set to True for testing if no suitable checkout file is available.
        This will print output to stdout.
    """
    checkout_file = config.get('checkout_file')
    versions_file = config.get('versions_file')

    log = Path(log)
    if data_root is None:
        # assuming normal OMICS pipeline conventions
        data_root = log.parent.parent.parent / 'data'

    if not data_root.is_dir():
        raise FileNotFoundError(f'no such directory: {data_root}')

    outputs, rules_run = get_info_from_log(log)
    if outputs and checkout_file or dry_run:
        stats = update_omics_checkout(
            outputs,
            data_root,
            checkout_file=checkout_file,
            dry_run=dry_run,
        )
        if checkout_file:
            print(f'[{checkout_file} OK]', end='', flush=True)
        if any(stats.values()):
            print(f' <-- same={stats["no_change"]} new={stats["new"]} '
                  f'ignored={stats["ignored"]} missing={stats["missing"]}',
                  end='', flush=True)
            if stats['errors']:
                print(f' *** errors={stats["errors"]} ***')
            else:
                print()
        else:
            print(' (no output files)')
    else:
        print('[post production] no output files per log')

    # FIXME TODO NOTE DEV ONLY
    if not rules_run:
        rules_run = ['get_reads_prep', 'get_reads_single']

    if versions_file:
        if rules is None:
            print('[WARNING] rules parameter not supplied, versions file will '
                  'not be updated')
        else:
            envs = []
            for name, rule in rules._rules.items():
                rule = rule.rule  # 1st rule is proxy obj
                if name not in rules_run:
                    continue

                if env_file := rule.conda_env:
                    envs.append(Env(rule.workflow, env_file=env_file))
            if envs:
                update_versions_file(
                    envs,
                    Path(versions_file),
                    dry_run=dry_run,
                )


def post_production_replay(log_dir, data_root=None, checkout_file=None,
                           dry_run=False, skip_on_error=False):
    """
    Process all logs in given directory.

    log_dir:
        Directory containing snakemake log files, like .snakemake/log
    data_root:
        Corresponding output data directory for the workflow.
    checkout_file:
        The checkout file to be written (will be updated if it exists)
    skip_on_error [bool]:
        Skip a log file if certain errors are encountered.
    """
    log_dir = Path(log_dir)
    if not log_dir.is_dir():
        # as glob() doesn't mind non-existing dirs
        raise FileNotFoundError(f'no such directory: {log_dir}')

    if data_root is None:
        data_root = log_dir.parent.parent / 'data'

    for i in sorted(Path(log_dir).glob('*.snakemake.log')):
        print(f'Processing {i.name} ... ', end='', flush=True)
        try:
            post_production(i, checkout_file=checkout_file, dry_run=dry_run)
        except JobIdReused as e:
            if skip_on_error:
                print(f'{e.__class__.__name__}: {e} [SKIP]')
            else:
                raise
