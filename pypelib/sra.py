""" Module for dealings with the SRA """
import argparse
from contextlib import contextmanager
from io import BytesIO
import json
from os import environ
from random import uniform
from time import sleep

from Bio import Entrez

import defusedxml.ElementTree as ET

from .utils import save_error_file


Entrez.email = 'GLAMR-omics@umich.edu'
Entrez.tool = 'glamr'
if api_key := environ.get('NCBI_API_KEY'):
    Entrez.api_key = api_key


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


def get_entry_xml(accession, multi=False):
    """
    retrieve SRA entry by accession, return xml Element object

    Also saves XML as file.  This is for debugging and testing.
    """
    data = get_entry_raw(accession, multi=multi)
    xml_txt = data.read().decode()
    with open(f'{accession}.efetch.xml', 'w') as ofile:
        ofile.write(xml_txt.replace('><', '>\n<'))  # line break between tags
    return ET.fromstring(xml_txt)


def get_entry(accession, multi=False, auto_parse=True, slow=False):
    """ retrieve SRA entry by accession """
    with get_entry_raw(accession, multi=multi, slow=slow) as data:
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


def get_srr(accn, sample_type=None, slow=False):
    """
    Get the SRR for a given other accession

    accn: some SRA accession

    Even if an SRR is given, some sanity checks will be run before the same
    accession is returned.
    """
    strategy2sample_type = {
        'AMPLICON': 'amplicons',
        'WGS': 'metagenome',
    }
    entry = get_entry(accn, slow=slow, auto_parse=False, multi=True)
    epset = entry['EXPERIMENT_PACKAGE_SET']
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
                    # TODO: add these to map
                    raise NotImplementedError(f'missing: {strategy}', *errargs)
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
                        f'{accn}/{sample_type}: multiple experiment packages',
                        *errargs,
                    )
    else:
        expack = epset['EXPERIMENT_PACKAGE']

    if sample_type is not None:
        try:
            strategy = expack['EXPERIMENT']['DESIGN']['LIBRARY_DESCRIPTOR']['LIBRARY_STRATEGY']  # noqa:E501
        except KeyError as e:
            e.args = (*e.args, *errargs)
            raise
        try:
            if sample_type != strategy2sample_type[strategy]:
                raise RuntimeError(
                    f'{accn}: {sample_type=} does not match {strategy=}',
                    *errargs,
                )
        except KeyError:
            # TODO: add these to map, see above
            raise NotImplementedError(strategy, *errargs)

    if 'RUN_list' in expack['RUN_SET']:
        raise NotImplementedError(f'{accn}: multiple runs', *errargs)

    return expack['RUN_SET']['RUN']['accession']


def stringify(value):
    """ Helper to make a dict into a line of text """
    return ' '.join(i.strip("'{}: ") for i in str(value).split())


def compile_srr_info(accn):
    info = {}
    entry = get_entry(accn, auto_parse=False, multi=True)
    # FIXME: fails in case of multiple exp packages
    pkg = entry['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']
    exp = pkg['EXPERIMENT']
    info['experiment'] = exp['accession']
    info['study'] = exp['STUDY_REF']['accession']
    info['sample'] = exp['DESIGN']['SAMPLE_DESCRIPTOR']['accession']
    lib = exp['DESIGN']['LIBRARY_DESCRIPTOR']
    info['strategy'] = lib['LIBRARY_STRATEGY']
    info['source'] = lib['LIBRARY_SOURCE']
    info['selection'] = lib['LIBRARY_SELECTION']
    info['layout'] = lib['LIBRARY_LAYOUT']
    info['platform'] = stringify(exp['PLATFORM'])
    info['run'] = pkg['RUN_SET']['RUN']['accession']
    info['spots'] = pkg['RUN_SET']['RUN']['total_spots']
    info['bases'] = pkg['RUN_SET']['RUN']['total_bases']

    print(json.dumps(info, indent=4))


def main():
    argp = argparse.ArgumentParser(description=__doc__)
    subs = argp.add_subparsers(dest='cmd', required=True)
    subp1 = subs.add_parser(
        'srr',
        help='Look up SRRxxx accession from other SRA accesions',
    )
    subp1.add_argument('any_accession')
    subp1.add_argument('--error-file', help='path where to save error infos')
    subp2 = subs.add_parser(
        'runinfo',
        help='Get info for a run, writes json data to stdout.',
    )
    subp2.add_argument(
        'accession',
        help='SRR or SRS accession',
    )
    args = argp.parse_args()

    match args.cmd:
        case 'srr':
            with save_error_file(args.error_file):
                print(get_srr(args.any_accession))
        case 'runinfo':
            compile_srr_info(args.accession)
        case _:
            argp.error('invalid subcommand')


if __name__ == '__main__':
    main()
