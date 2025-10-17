""" Module for dealings with the SRA """
import argparse
from io import BytesIO
import json
from os import environ

from Bio import Entrez

import defusedxml.ElementTree as ET


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
            # turn key/value into into key/list
            lists[child_tag] = [body.pop(child_tag)]
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


def search(accession):
    """ Run SRA search and get list of SRA-internal IDs """
    with Entrez.esearch(db='sra', term=accession, idtype='acc') as h:
        res = Entrez.read(h)
        return res['IdList']


def get_entry_raw(accession, multi=False):
    """ retrieve SRA entry by accession, return byte file handle """
    ids = search(accession)
    if len(ids) == 0:
        raise RuntimeError('search failed?')
    elif not multi and len(ids) > 1:
        raise RuntimeError(f'search returned multiple hits: {ids}')

    with Entrez.efetch(db='sra', id=','.join(ids)) as h:
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


def get_entry(accession, multi=False):
    """ retrieve SRA entry by accession """
    pref = accession[:3]
    with get_entry_raw(accession, multi=multi) as data:
        if pref in ['SRS', 'SRX', 'SRR', 'SRP']:
            # biopython parser won't work with the EXPERIMENT_PACKAGE_SET
            # that we're expecting here
            return el2dict(ET.fromstring(data.read().decode()))

        try:
            return Entrez.read(data)
        except ValueError as e:
            # no biopython support for the document type or so?
            print(f'auto-decoding failed for {accession}: '
                  f'{e.__class__.__name__}: {e}')


def srs2srr(srs_accession):
    """ Look up SRRxxxxxxx from SRSxxxxxxxx accessions """
    if not srs_accession.startswith('SRS'):
        raise ValueError(f'expected SRSxxxxx accession, got {srs_accession}')
    entry = get_entry(srs_accession)
    try:
        pkg = entry['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']
        run = pkg['RUN_SET']['RUN']
        return run['accession']
    except KeyError as e:
        raise LookupError(
            f'Failed parsing SRA record: {e}\n{str(entry)[:300]}...'
        ) from e


def stringify(value):
    """ Helper to make a dict into a line of text """
    return ' '.join(i.strip("'{}: ") for i in str(value).split())


def compile_srr_info(accn):
    if not accn.startswith('SRS') and not accn.startswith('SRR'):
        raise ValueError('expecting SRR or SRS accession')
    info = {}
    entry = get_entry(accn)
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
        'srs2srr',
        help='Look up SRRxxx accession from SRSxxx',
    )
    subp1.add_argument('SRS_accession')
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
        case 'srs2srr':
            print(srs2srr(args.SRS_accession))
        case 'runinfo':
            compile_srr_info(args.accession)
        case _:
            argp.error('invalid subcommand')


if __name__ == '__main__':
    main()
