""" Module for dealings with the SRA """
import argparse
from contextlib import contextmanager
from io import BytesIO
import json
from os import environ
from random import uniform
from time import sleep
from xml.etree.ElementTree import indent as xml_indent

from Bio import Entrez

import defusedxml.ElementTree as ET


Entrez.email = 'GLAMR-omics@umich.edu'
Entrez.tool = 'glamr'
if api_key := environ.get('NCBI_API_KEY'):
    Entrez.api_key = api_key


def set_api_key(key):
    # allow the env to override this
    if key and not Entrez.api_key:
        Entrez.api_key = key


def el2dict(elem):
    """ Convert ElementTree into a dict """
    if (text := elem.text) is not None:
        text = text.strip()  # fix for single newlines

    body = dict(**elem.attrib)
    lists = {}
    for child in (el2dict(i) for i in elem):
        child_tag, child_body = list(child.items())[0]
        if child_tag in body:
            # first turn key/value into into key/list
            lists[child_tag] = [body.pop(child_tag)]
            lists[child_tag].append(child_body)
        elif child_tag in lists:
            lists[child_tag].append(child_body)
        else:
            body[child_tag] = child_body

    for tag, children in lists.items():
        body[tag + '_list'] = children
    del lists

    if body:
        if text is None:
            if len(body) == 1:
                tag, value = list(body.items())[0]
                if value is None:
                    # remove the useless None value
                    body = tag
        else:
            # Empty text is preserved but not None
            body['text'] = text
    elif text is not None:
        body = text
    else:
        # don't return empty dict
        body = None

    return {elem.tag: body}


@contextmanager
def try_three_times(fn, *args, **kwargs):
    """
    Slow down Entrez access if "too many requests," try up to three times.

    Use as this:
        with try_three_times(Entrez.esearch, db='sra', ...) as handle:
            data = Entrez.read(handle)
    """
    attempts_left = 3
    while attempts_left:
        try:
            resource = fn(*args, **kwargs)
        except Exception as e:
            attempts_left -= 1
            if 'Too Many Requests' in str(e):
                # Try to catch urllib.error.HttpError 429? but can't seem to
                # import that here (maybe circular?)
                if attempts_left:
                    sleep(uniform(10.0, 30.0))
                continue
            # give up immediately on other errors
            raise
        else:
            with resource as resource:
                yield resource
            break
    else:
        raise RuntimeError('tried maximum number of attempts')


def search(accession, with_errors=False, quiet=False):
    """ Run SRA search and get list of SRA-internal IDs """
    kwargs = dict(db='sra', term=f'{accession}[accn]')
    with try_three_times(Entrez.esearch, **kwargs) as h:
        res = Entrez.read(h)

    errs = []
    for key, value in res.get('ErrorList', {}).items():
        if value:
            if not isinstance(value, list):
                value = [value]
            for item in value:
                errs.append(f'[ERROR] (Entrez.esearch) {key}: {item}')
    for key, value in res.get('WarningList', {}).items():
        if key == 'OutputMessage':
            for item in value:
                errs.append(f'[NOTICE] (Entrez.esearch) {item}')
        elif value:
            if not isinstance(value, list):
                value = [value]
            for item in value:
                errs.append(f'[WARNING] (Entrez.esearch) {key}: {item}')

    if not quiet:
        for line in errs:
            print(line)

    if with_errors:
        return res['IdList'], errs
    else:
        return res['IdList']


def get_entry_raw(accession, multi=False, slow=False):
    """ retrieve SRA entry by accession, return byte file handle """
    if slow:
        sleep(uniform(0.0, 10.0))
    ids, search_errors = search(accession, with_errors=True)
    if len(ids) == 0:
        raise RuntimeError('search failed?', *search_errors)
    elif not multi and len(ids) > 1:
        raise RuntimeError(
            f'search for {accession=} returned multiple hits: {ids}',
            *search_errors,
        )

    if slow:
        sleep(uniform(0.0, 10.0))
    with try_three_times(Entrez.efetch, db='sra', id=','.join(ids)) as h:
        data = BytesIO()
        data.write(h.read())
        data.seek(0)
        return data


def get_entry(accession=None, file=None, multi=False, auto_parse=False,
              slow=False):
    """ retrieve SRA entry by accession or from saved file """
    if file is not None and accession is None:
        data = open(file, 'rb')
    elif accession is not None and file is None:
        data = get_entry_raw(accession, multi=multi, slow=slow)
    else:
        raise TypeError('either accession if file must be given, but not both')

    with data as data:
        if not auto_parse:
            try:
                # biopython parser won't work with the EXPERIMENT_PACKAGE_SET
                # that we're expecting here
                entry = el2dict(ET.fromstring(data.read().decode()))
            except Exception as e:
                print(f'ERROR: e2dict failed! {e.__class__.__name__}: {e}')
                raise
            else:
                if isinstance(entry, dict):
                    if 'EXPERIMENT_PACKAGE_SET' in entry:
                        return entry
                # fall-back to Biopython's parsing
                print(
                    f'NOTICE: {accession=} did not return a exp pkg set, but:'
                    f'\n', str(entry)[:1000]
                )
                data.seek(0)

        try:
            return Entrez.read(data)
        except ValueError as e:
            # no biopython support for the document type or so?
            print(f'auto-decoding failed for {accession}: '
                  f'{e.__class__.__name__}: {e}')
            raise


def get_experiment(accn, data, sample_type=None):
    """
    Get the experiment record for a given accession

    accn: The SRA accession used to retrieve the data.
    data: return of get_entry()
    """
    strategy2sample_type = {
        'AMPLICON': 'amplicons',
        'WGS': 'metagenome',
    }
    strict_mode = environ.get('SRA_ESEARCH_STRICT', False)
    epset = data['EXPERIMENT_PACKAGE_SET']
    errargs = (accn, sample_type, epset)

    if 'EXPERIMENT_PACKAGE_list' in epset:
        # multiple experiments, need to find the right one
        sample_type_matches = []
        for expack in epset['EXPERIMENT_PACKAGE_list']:
            if accn.startswith('SRR'):
                if 'RUN' in expack['RUN_SET']:
                    if accn == expack['RUN_SET']['RUN']['accession']:
                        break
                else:
                    for run in expack['RUN_SET']['RUN_list']:
                        if accn == run['accession']:
                            break  # inner
                    else:
                        continue  # next expack
                    break  # outer
            elif accn.startswith('SRX'):
                if expack['EXPERIMENT']['accession'] == accn:
                    break
            elif sample_type is not None:
                strategy = expack['EXPERIMENT']['DESIGN']['LIBRARY_DESCRIPTOR']['LIBRARY_STRATEGY']  # noqa:E501
                try:
                    if sample_type == strategy2sample_type[strategy]:
                        sample_type_matches.append(expack)
                except KeyError:
                    if strict_mode:
                        # TODO: add these to map
                        raise NotImplementedError(
                            f'unknown library strategy: {strategy}', *errargs
                        )
            else:
                # TODO
                raise NotImplementedError(
                    'need a way to discriminate these', *errargs,
                )
        else:
            match len(sample_type_matches):
                case 0:
                    raise RuntimeError(
                        f'{accn}: no matching experiment package found',
                        *errargs
                    )
                case 1:
                    expack = sample_type_matches[0]
                case _:
                    raise RuntimeError(
                        f'{accn}/{sample_type}: multiple matching experiment '
                        f'packages',
                        *errargs,
                    )
    else:
        expack = epset['EXPERIMENT_PACKAGE']

    if sample_type is not None:
        strategy = expack['EXPERIMENT']['DESIGN']['LIBRARY_DESCRIPTOR']['LIBRARY_STRATEGY']  # noqa:E501
        try:
            if sample_type != strategy2sample_type[strategy]:
                if strict_mode:
                    raise RuntimeError(
                        f'{accn}: {sample_type=} does not match {strategy=}',
                        *errargs,
                    )
        except KeyError:
            if strict_mode:
                # TODO: add these to map, see above
                raise NotImplementedError(strategy, *errargs)

    if 'RUN_list' in expack['RUN_SET']:
        raise NotImplementedError(f'{accn}: multiple runs', *errargs)

    return expack


def stringify(value):
    """ Helper to make a dict into a line of text """
    return ' '.join(i.strip("'{}: ") for i in str(value).split())


def compile_srr_info(expack):
    info = {}
    exp = expack['EXPERIMENT']
    info['experiment'] = exp['accession']
    info['study'] = exp['STUDY_REF']['accession']
    info['sample'] = exp['DESIGN']['SAMPLE_DESCRIPTOR']['accession']
    lib = exp['DESIGN']['LIBRARY_DESCRIPTOR']
    info['strategy'] = lib['LIBRARY_STRATEGY']
    info['source'] = lib['LIBRARY_SOURCE']
    info['selection'] = lib['LIBRARY_SELECTION']
    info['layout'] = lib['LIBRARY_LAYOUT']
    info['platform'] = stringify(exp['PLATFORM'])
    info['run'] = expack['RUN_SET']['RUN']['accession']
    info['spots'] = expack['RUN_SET']['RUN']['total_spots']
    info['bases'] = expack['RUN_SET']['RUN']['total_bases']
    return info


def main():
    argp = argparse.ArgumentParser(description=__doc__)
    subs = argp.add_subparsers(dest='cmd', required=True)
    subs.add_parser(
        'xml',
        help='Print raw xml search result to stdout',
    )
    subp1 = subs.add_parser(
        'srr',
        help='Look up SRRxxx accession from other SRA accesions',
    )
    subp1.add_argument('--sample_type', help='a sample type')
    subp2 = subs.add_parser(
        'runinfo',
        help='Get info for a run, writes json data to stdout.',
    )
    subp2.add_argument('--sample_type', help='a sample type')
    argp.add_argument('accession', help='An SRA accession')
    args = argp.parse_args()

    match args.cmd:
        case 'xml':
            xml = get_entry_raw(args.accession, multi=True, slow=False)
            tree = ET.fromstring(xml.read().decode())
            xml_indent(tree)
            print(ET.tostring(tree).decode())
        case 'srr':
            data = get_entry(accession=args.accession, multi=True,
                             slow=False)
            expack = get_experiment(args.accession, data, args.sample_type)
            print(expack['RUN_SET']['RUN']['accession'])
        case 'runinfo':
            data = get_entry(accession=args.accession, multi=True,
                             slow=False)
            expack = get_experiment(args.accession, data, args.sample_type)
            print(json.dumps(compile_srr_info(expack), indent=4))
        case _:
            argp.error('invalid subcommand')


if __name__ == '__main__':
    main()
