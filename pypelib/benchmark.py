"""
Plot benchmark data

This can do

 (1) Plot the benchmark data collected during post-production.  These plots
     will include indication of job sizes.  There will be one plot per rule.

 (2) Plot the data straight from the benchmark directory without job sizes.
     There will be one plot per subdirectory.
"""
import argparse
import math
from pathlib import Path

from matplotlib.backends.backend_pdf import PdfPages
import pandas


def cli():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument(
        'benchmarks',
        help='Your post-production benchmark collection file.',
    )
    argp.add_argument(
        '-o', '--output',
        help='Name of output PDF file.',
    )
    args = argp.parse_args()
    if args.output:
        output = args.output
    else:
        output = Path(args.benchmarks).stem + '.pdf'

    inpath = Path(args.benchmarks)
    if inpath.is_file():
        # post-production-collected benchmark table
        plot_all_table(inpath, output)
    elif inpath.is_dir():
        # assume this is the benchmark directory
        plot_all_plain(inpath, output)
    else:
        argp.error(f'no such file or directory: {inpath}')


def load_plain_table(path, colums):
    with open(path) as ifile:
        head = ifile.readline().rstrip().split('\t')
        row = ifile.readline().rstrip().split('\t')
        if ifile.read(1):
            # rule was run with benchmark=repeat(..., >1)
            raise RuntimeError('only expecting two lines')

        try:
            data = {key: row[head.index(key)] for key in colums}
        except ValueError as e:
            raise RuntimeError(
                f'Failed parsing {path}: missing column? {e}'
            ) from e
        except IndexError as e:
            raise RuntimeError(
                f'Failed parsing {path}: missing data values? {e}'
            ) from e

        try:
            data = {
                k: float('nan') if v == 'NA' else float(v)
                for k, v in data.items()
            }
        except ValueError as e:
            raise RuntimeError(
                f'Failed parsing {path}: expect a (FP) number: {e}'
            ) from e

        return data


def load_extended_table(path):
    df = pandas.read_table(
        path,
        parse_dates=['timestamp'],
        dtype=dict(rule='category'),
    )
    return df


def plot_rule(df, rule_name):
    if 'rule' in df:
        df = df.loc[df['rule'] == rule_name]

    do_size_legend = 'jobsize_log' in df

    ax = df.plot.scatter(
        's', 'max_rss', c='b',
        s='jobsize_log' if do_size_legend else None,
        alpha=0.6,
        # needed to enable legend stuff below
        label='jobsize' if do_size_legend else None
    )
    ax.set_title(rule_name)

    if do_size_legend:
        handles, labels = ax.get_legend_handles_labels()
        # handles[0] is the PathCollection instance that ax.scatter() would
        # return
        handles, labels = handles[0].legend_elements(
            prop='sizes',
            num=4,
            alpha=0.6,  # same as for plot
            color='b',  # same as for plot
        )
        ax.legend(
            list(reversed(handles)),
            list(reversed(labels)),
            title='Job sizes',
        )
    return ax


def plot_all_table(inpath, outpath):
    """ plot data from post-production benchmark table """
    df = load_extended_table(inpath)
    df['jobsize_log'] = df['jobsize'].apply(math.log2) * 10

    with PdfPages(outpath) as pp:
        for rule in df.rule.cat.categories:
            print(f'Plotting for rule {rule} ...')
            ax = plot_rule(df, rule)
            pp.savefig(ax.figure)
            ax.cla()
        breakpoint()
    print(f'Saved as {outpath}')


def get_subdirs(parent):
    for path in parent.iterdir():
        if path.is_dir():
            yield path
            yield from get_subdirs(path)


def plot_all_plain(inpath, outpath):
    """ plot data from (per-rule?) benchmark directories """
    bm_columns = ['s', 'max_rss']

    with PdfPages(outpath) as pp:
        for subdir in get_subdirs(inpath):

            name = '_'.join(subdir.relative_to(inpath).parts)

            data = {key: [] for key in bm_columns}
            for bm_file in subdir.iterdir():
                if not bm_file.is_file():
                    continue
                for k, v in load_plain_table(bm_file, bm_columns).items():
                    data[k].append(v)

            df = pandas.DataFrame(data)
            print(f'Plotting {name} ...')
            ax = plot_rule(df, name)
            pp.savefig(ax.figure)
            ax.cla()
    print(f'Saved as {outpath}')


if __name__ == '__main__':
    cli()
