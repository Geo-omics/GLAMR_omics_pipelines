"""
Post-workflow code
"""
import argparse
import bisect
from contextlib import ExitStack
from datetime import datetime, UTC
from io import BytesIO, StringIO
import json
from pathlib import Path
import re
import shutil
from subprocess import run, PIPE
import sys

import yaml

from .utils import PipelineVersion


BASE_SNAKEMAKE_CONDA_ENV = 'config/conda_yaml/snakemake.yaml'
""" conda env that snakemake runs under """


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

    Returns a dict mapping rule names to list of output files.
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
    version_pat = re.compile(
        rf'^{PipelineVersion.LOG_MSG_PREFIX}(?P<version>.*)$',
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

    if m := version_pat.search(log_txt):
        pl_version = m['version']
    else:
        pl_version = None

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

    data = {}
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
        rulename = jobs[jobid]['name']
        if rulename not in data:
            data[rulename] = []
        data[rulename] += jobs[jobid]['output']

    return data, pl_version


def read_omics_checkout(checkout_file):
    """
    Helper to read and parse the checkout file and return dict mapping files to
    most recent mtime.
    """
    MIN_TIME = datetime.min.replace(tzinfo=UTC)
    most_recent = {}
    with open(checkout_file) as ifile:
        for line in ifile:
            mtime, version, rule, path = line.rstrip('\n').split('\t')
            mtime = datetime.fromisoformat(mtime)
            if most_recent.get(path, MIN_TIME) < mtime:
                most_recent[path] = mtime
    return most_recent


def update_omics_checkout(
    output_files,
    data_root,
    version_before_mtime=False,
    checkout_file=None,
    dry_run=False,
):
    """
    Update checkout file by appending one line for each new or newer file.

    output_files:
        Dict mapping rule name to list of files.

    checkout_file:
        Path to new or existing checkout file.  If this is None, then the
        output file will be written to stdout instead.

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

    output_files1 = []
    for rule, files in output_files.items():
        output_files1 += [(rule, i) for i in sorted(files)]

    if version_before_mtime:
        # Get list of pairs (commit time, commit hash)
        cmd = ['git', 'log', '--format=%cI %h']
        p = run(cmd, check=True, stdout=PIPE)
        git_commits = [
            line.strip().split(' ', maxsplit=1)
            # git log has most recents first, reverse for bisect's benefit
            for line in reversed(p.stdout.decode().splitlines())
        ]
        git_commits = [
            (datetime.fromisoformat(t), h)
            for t, h in git_commits
        ]
        version_str_cache = {}
    else:
        version = PipelineVersion.current()

    missing = no_change = ignored = errors = 0
    new_data = []
    for rule, file00 in output_files1:
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

        if version_before_mtime:
            idx = bisect.bisect(git_commits, mtime, key=lambda x: x[0])
            commit = git_commits[idx - 1][1]
            if commit not in version_str_cache:
                try:
                    version_str_cache[commit] = PipelineVersion.from_commit(commit)  # noqa:E501
                except Exception:
                    print(f'[DEBUG] getting pipeline version for {commit} for '
                          f'{mtime=} {rule=} {file00=}')
                    raise
            version = version_str_cache[commit]

        new_data.append((mtime, version, rule, relpath))

    with ExitStack() as estack:
        if checkout_file is None or dry_run:
            ofile = sys.stdout
        else:
            ofile = estack.enter_context(open(checkout_file, 'a'))
            print(f'Writing {checkout_file} ...', end='', flush=True)

        for mtime, version, rule, relpath in new_data:
            ofile.write(f'{mtime}\t{version or ""}\t{rule or ""}\t{relpath}\n')

    if checkout_file and not dry_run:
        print('[OK] ', end=' ')  # expect stats printed after this

    return {
        'new': len(new_data),
        'missing': missing,
        'no_change': no_change,
        'ignored': ignored,
        'errors': errors,
    }


class NoConda(Exception):
    pass


def get_conda_env_package_versions(env):
    """
    Get installed versions of important packages in given conda environment

    env:
        A snakemake.deployment.conda.Env object or path to a yaml formatted
        conda spec file.  If a spec file was passed, we assume that it
        corresponds to the currently active conda environment.  There is no
        good way to check this, though it seems the information is availale in
        conda-meta/history.

    Returns dict of package names mapping to version. Raises NoConda if
    we're not running inside an active conda environment.  Not being inside an
    env is fine if env is a snakemake Env.
    """
    pkgname_pat = re.compile(r'^([\w0-9]+)(.*)')
    packages = []

    if isinstance(env, Path):
        spec_ifile = open(env)
        prefix = None  # don't know environment location
    else:
        # assume snakemake Env
        spec_ifile = StringIO(env.content.decode())
        prefix = env.address  # this fails if conda not installed

    with spec_ifile:
        spec = yaml.safe_load(spec_ifile)

    for j in spec['dependencies']:
        m = pkgname_pat.match(j)
        pkg, version = m.groups()
        packages.append(pkg)

    if not packages:
        raise RuntimeError('parsing env yaml failed, no packaged')

    if prefix:
        prefix_args = ['--prefix', prefix]
    else:
        prefix_args = []
    regex = '^(' + '|'.join(packages) + ')$'

    if shutil.which('conda') is None:
        # avoid raising potentially non-specific error in run()
        raise NoConda('conda is not in PATH')

    if not prefix:
        info_proc = run(['conda', 'info', '--json'], check=True, stdout=PIPE)
        with BytesIO(info_proc.stdout) as buf:
            info = json.load(buf)
        if not info.get('active_prefix'):
            # avoid misleading message to stderr by "conda list"
            raise NoConda('not in an active conda env')

    cmd = ['conda', 'list'] + prefix_args + ['--export', regex]
    p = run(cmd, check=True, stdout=PIPE)
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


class EnvPackageInfo:
    """
    Information on a conda environment

    Installed package versions and which snakemake rules used the env.
    """
    def __init__(self, env_name, rules=None, packages=None):
        if not env_name:
            raise ValueError('env name must be provided')
        self.name = env_name
        if rules:
            self.rules = tuple()
            self.add_rules(rules)
        else:
            if env_name == Path(BASE_SNAKEMAKE_CONDA_ENV).stem:
                self.rules = None
            else:
                raise ValueError('at least one rule must be provided')
        if not packages:
            raise ValueError('some package versions must be provided')
        self.packages = dict(sorted(packages.items()))

    def to_dict(self):
        """ helper for json encoding """
        return dict(rules=self.rules, packages=self.packages)

    def add_rules(self, rules):
        """
        add more rules

        Ensures order and no duplicates
        """
        if rules:
            self.rules = tuple(sorted(list(self.rules) + list(rules)))

    def update(self, other):
        """
        Update from another instance

        Rules are added, the package versions get replaced.
        """
        if self.name != other.name:
            raise ValueError('other name must be the same')
        self.rules = tuple(sorted(set(self.rules + other.rules)))
        self.packages = other.packages


class VersionInfoFile:
    """
    Conda env versioning history

    The history of installed package versions of conda envs used in Snakefile
    rules.
    """
    def __init__(self, data=None):
        if data is None:
            data = {}
        if list(data.keys()) != sorted(data.keys()):
            raise RuntimeError('pipeline versions are not well-ordered')
        self.data = data
        self.changed = False

    @classmethod
    def from_file(cls, version_file_path):
        with open(version_file_path) as ifile:
            file_data = json.load(ifile)

        data = {}
        for pl_version, data_part in file_data.items():
            pl_version = PipelineVersion(pl_version)
            data[pl_version] = {}
            for env_name, env_data in data_part.items():
                data[pl_version][env_name] = \
                    EnvPackageInfo(env_name, **env_data)

        return cls(data)

    def get_current(self, env_name):
        """
        look up most recent pkg info for given env name

        Returns None if the env is not yet listed.
        """
        rules = None
        packages = None
        for plver, data in self.data.items():
            if env_data := data.get(env_name):
                if env_data.rules is None:
                    # usually a base env
                    rules = None
                else:
                    # normal env for rules
                    if rules is None:
                        rules = set()
                    rules.update(env_data.rules)
                packages = env_data.packages

            elif rules is not None:
                # Remove any rules that are now used by other env.  There can
                # be only one (env per rule.)
                rules.difference_update([i.rules for i in data.values()])

        if packages is None:
            # env not listed yet
            return None

        if rules or rules is None:
            return EnvPackageInfo(env_name, rules=rules, packages=packages)
        else:
            # env not used by any rules (anymore)
            return None

    def update(self, pipeline_version, env_pkg_info):
        """
        Incorporate given information on a single conda env.
        """
        if self.data:
            last_version = list(self.data)[-1]
        else:
            last_version = None
        if last_version and pipeline_version < last_version:
            # rebasing will get us here, too
            raise ValueError(
                f'can\'t update data for old pipeline version: '
                f'{last_version=} {pipeline_version=}'
            )

        if current := self.get_current(env_pkg_info.name):
            if env_pkg_info.packages == current.packages:
                if env_pkg_info.rules is None and current.rules is None:
                    # nothing to be done
                    return
                if set(env_pkg_info.rules).issubset(current.rules):
                    # also nothing to be done
                    return

        # something changed!
        if pipeline_version == last_version:
            # get old data for update
            data = self.data[pipeline_version]
            if env_pkg_info.name in data:
                # keeping those rules even if packages/versions change
                env_pkg_info.add_rules(*data[env_pkg_info.name].rules)
        else:
            # data for new pipline version to be added
            data = {}

        data[env_pkg_info.name] = env_pkg_info
        # re-sort in case env was new to existing pipeline version
        data = dict(sorted(data.items()))
        self.data[pipeline_version] = data
        self.changed = True

    @staticmethod
    def json_default(obj):
        """ helper for json encoding with dump() """
        if isinstance(obj, EnvPackageInfo):
            return obj.to_dict()
        else:
            raise TypeError('don\'t know how to handle this type')

    def dump(self):
        """ return json text from instance """
        # convert the keys (PipelineVersion instance) to str, dumps() can't do
        # that currently, default is not passed on to process dict keys (maybe
        # in future python version)
        data = {str(k): v for k, v in self.data.items()}
        return json.dumps(data, indent=4, default=self.json_default)

    def save(self, output_file):
        """ save data to file """
        json_text = self.dump()
        with open(output_file, 'w') as ofile:
            ofile.write(json_text)
            ofile.write('\n')


def update_versions_file(workflow, dry_run=False):
    """
    Make new or update an existing version file

    This will only work from within a proper snakemake run.
    """
    versions_file = workflow.config.get('version_file')
    if versions_file is None and not dry_run:
        # not configured, nothing to be done
        return

    try:
        from snakemake.deployment.conda import Env
        from snakemake.settings.types import DeploymentMethod
    except ImportError:
        print('[WARNING] Snakemake not installed, versions file not updated')
        return

    if shutil.which('conda') is None:
        print('[WARNING] conda is not installed, version file not updated')
        # return

    if DeploymentMethod.CONDA not in workflow.deployment_settings.deployment_method:  # noqa:E501
        print('[WARNING] conda is not a deployment method, versions file '
              'not updated')
        return

    # get pipeline version
    pipeline_version = PipelineVersion.current()

    if versions_file is not None:
        versions_file = Path(versions_file)

    if versions_file and versions_file.is_file():
        # get existing data from file
        vinfo = VersionInfoFile.from_file(versions_file)
    else:
        # start from scratch
        vinfo = VersionInfoFile()

    # distinct conda-using rules used for finished jobs
    rules = set(
        job.rule for job
        in workflow.dag.finished_jobs
        if job.rule.conda_env
    )

    envs = [Path(BASE_SNAKEMAKE_CONDA_ENV)]
    env_rules_map = {}  # tracks which rules use which conda env
    for rule in rules:
        if rule.conda_env in env_rules_map:
            env_rules_map[rule.conda_env].append(rule.name)
        else:
            env_rules_map[rule.conda_env] = [rule.name]
            envs.append(Env(rule.workflow, env_file=rule.conda_env))

    # get currently installed package versions
    for env in envs:
        if isinstance(env, Path):
            env_name = env.stem
            env_rules = None
        else:
            # snakemake Env
            env_name = Path(env.file.path).stem
            env_rules = env_rules_map[env.file.path]

        try:
            packages = get_conda_env_package_versions(env)
        except NoConda as e:
            print(f'[WARNING] {e}: not updating versions file for env '
                  f'"{env_name}"')
            continue

        if not packages:
            continue

        vinfo.update(
            pipeline_version,
            EnvPackageInfo(env_name, rules=env_rules, packages=packages)
        )

    if vinfo.changed:
        print(f'Writing {versions_file} ...', end='', flush=True)
        if dry_run or not versions_file:
            print()
            print(vinfo.dump())
            print('[dry run OK]')
        else:
            vinfo.save(versions_file)
            print('[OK]')
    else:
        print('[OK] no changes to versions file')


def load_benchmark(path):
    with open(path) as ifile:
        head = ifile.readline().split()
        row = ifile.readline().split()
        return dict(zip(head, row, strict=True))


def collect_benchmarks(workflow, dry_run=False):
    """
    Collect benchmarks from jobs into single table with added job data

    This does nothing unless 'collect_benchmark' set in the config.
    """
    if outpath := workflow.config.get('collect_benchmarks'):
        outpath = Path(outpath)
    elif dry_run:
        outpath = None
    else:
        return

    # columns we're adding
    base_cols = ['rule', 'wildcards', 'timestamp', 'threads', 'jobsize']
    # columns from benchmarks
    bm_cols = ['s', 'h:m:s', 'max_rss', 'max_vms', 'max_uss', 'max_pss',
               'io_in', 'io_out', 'mean_load', 'cpu_time']
    rows = []
    for j in filter(lambda x: bool(x.benchmark), workflow.dag.finished_jobs):
        bm_file = Path(j.benchmark)
        jobsize = sum(
            Path(i).stat().st_size
            for i in j.input
        )
        jobstats = load_benchmark(bm_file)
        bm_mtime = datetime.fromtimestamp(bm_file.stat().st_mtime)
        rows.append([
            j.name,
            ','.join(j.wildcards),
            bm_mtime.isoformat(),
            str(j.threads),
            str(jobsize),
            *(jobstats[col] for col in bm_cols)
        ])

    with ExitStack() as estack:
        if dry_run:
            ofile = None
        else:
            ofile = estack.enter_context(outpath.open('a'))

        if dry_run or ofile.tell() == 0:
            # write header if it's a new file
            print(*base_cols, *bm_cols, sep='\t', file=ofile)

        for row in rows:
            print(*row, sep='\t', file=ofile)

    if not dry_run:
        print(f'[OK] {len(rows)} benchmarks written to {outpath}')


def post_production(log, workflow=None, *, data_root=None, checkout_file=None,
                    dry_run=False):
    """
    Post-snakemake run processing

      * Update the output-file-checkout file
      * Update the package version tracking file

    Intended to be called via onsuccess and onerror.  Beware, makes use of
    snakemake internals.

    log:
        str or PathLike to snakemake log file.
    workflow:
        Workflow instance passed from Snakefile.  If this is None, then some
        functionality will be disabled.
    data_root:
        The data directory.  If None, this will be derived from path to log
        file.
    checkout_file:
        str or PathLike to two-column tab-separated text file listing output
        files and their modification times.  This can be None if no such file
        is configured.  This will take precedence if a workflow is also passed
        and has the checkout_file config setting set.
    dry_run [bool]:
        Set to True for testing if no suitable checkout file is available.
        This will print output to stdout.
    """
    if workflow:
        if not any((
            workflow.config.get('checkout_file'),
            workflow.config.get('versions_file'),
            workflow.config.get('collect_benchmarks'),
        )):
            # nothing to do
            return

    if workflow:
        checkout_file = workflow.config.get('checkout_file')

    if checkout_file or dry_run:
        log = Path(log)
        if data_root is None:
            # assuming normal OMICS pipeline conventions
            data_root = log.parent.parent.parent / 'data'

        if not data_root.is_dir():
            raise FileNotFoundError(f'no such directory: {data_root}')

        if workflow:
            outputs = {j.name: j.output for j in workflow.dag.finished_jobs}
            pl_version = PipelineVersion.current()
        else:
            # replay from log file
            outputs, pl_version = get_info_from_log(log)

        if outputs and checkout_file or dry_run:
            stats = update_omics_checkout(
                outputs,
                data_root,
                version_before_mtime=pl_version is None,
                checkout_file=checkout_file,
                dry_run=dry_run,
            )
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

    if workflow and workflow.config.get('version_file'):
        update_versions_file(workflow, dry_run=dry_run)

    if workflow:
        collect_benchmarks(workflow, dry_run=dry_run)


def replay(log_dir, data_root=None, checkout_file=None, dry_run=False,
           skip_on_error=False):
    """
    Process all logs in given directory.  Compile a checkout file.

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
            post_production(
                i,
                workflow=None,
                checkout_file=checkout_file,
                dry_run=dry_run,
            )
        except JobIdReused as e:
            if skip_on_error:
                print(f'{e.__class__.__name__}: {e} [SKIP]')
            else:
                raise


def test_conda_envs(path):
    from snakemake.api import (
        OutputSettings,
        ResourceSettings,
        SnakemakeApi,
    )
    from snakemake.deployment.conda import Env

    with SnakemakeApi(
        OutputSettings(
            verbose=True,
            show_failed_logs=True,
        ),
    ) as snakemake_api:
        workflow_api = snakemake_api.workflow(
            resource_settings=ResourceSettings(),
            snakefile=path,
        )
        envs = {}
        for i in workflow_api._workflow.rules:
            if i.conda_env and i.conda_env not in envs:
                envs[i.conda_env] = \
                    Env(i.workflow, i.conda_env, envs_dir='.snakemake/conda')

        for path, env in sorted(envs.items()):
            name = Path(path).stem
            if not Path(path).is_file():
                print(f'{name}: no such spec file: {path}')
            elif not Path(env.address).is_dir():
                print(f'{name}: not such directory: {env.address}')
            else:
                vers = get_conda_env_package_versions(env)
                print(f'{name}: {vers}')


if __name__ == '__main__':
    """ command-line interface for replay function and dev/testing """
    argp = argparse.ArgumentParser(
        description='Run the post-production replay function'
    )
    subps = argp.add_subparsers(dest='cmd', required=True)

    subp1 = subps.add_parser(
        'replay',
        help='Replay post production for multiple snakemake log files',
    )
    subp1.add_argument(
        'log_directory',
        nargs='?',
        default='./.snamemake/log',
        help='Snakemake log directory, defaults to "./snakemake/log"',
    )
    subp1.add_argument(
        '--data-root', default='./data',
        help='The data directory',
    )
    subp1.add_argument(
        '--checkout-file',
        help='Write file listing to this file. By default everything is '
             'written to stdout',
    )
    subp1.add_argument(
        '--dry-run', '-n', action='store_true',
        help='Do not change any files'
    )

    subp2 = subps.add_parser(
        'single',
        help='Test post-production on a single log file',
    )
    subp2.add_argument(
        'log',
        help='A snakemake log file',
    )
    subp2.add_argument(
        '--data-root', default='./data',
        help='The data directory',
    )
    subp2.add_argument(
        '--checkout-file',
        help='Write file listing to this file. By default everything is '
             'written to stdout',
    )
    subp2.add_argument(
        '--dry-run', '-n', action='store_true',
        help='Do not change any files'
    )

    subp3 = subps.add_parser(
        'conda',
        help='Test getting installed package versions from conda environments',
    )
    subp3.add_argument(
        'snakefile', default='./Snakefile',
        nargs='?',
        help='A Snakefile, defaults to "./Snakefile"',
    )

    args = argp.parse_args()
    match args.cmd:
        case 'replay':
            path = Path(args.log_directory)
            if not path.is_dir():
                argp.error(f'no such directory: {path}')
            replay(
                path,
                args.data_root,
                checkout_file=args.checkout_file,
                dry_run=args.dry_run,
                skip_on_error=True,
            )
        case 'single':
            log = Path(args.log)
            if not log.is_file():
                argp.error(f'no such file: {log}')
            config = {}
            post_production(
                log,
                workflow=None,
                data_root=Path(args.data_root),
                checkout_file=args.checkout_file,
                dry_run=args.dry_run,
            )
        case 'conda':
            snakefile = Path(args.snakefile)
            if not snakefile.is_file():
                argp.error(f'no such file: {snakefile}')
            test_conda_envs(snakefile)
        case _:
            argp.error('invalid subcommand')
