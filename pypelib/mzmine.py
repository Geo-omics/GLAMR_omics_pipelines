""" mzmine-related things """
from pathlib import Path
from xml.etree.ElementTree import ElementTree, indent, SubElement

import defusedxml.ElementTree


def get_el(root, part, parameter_name):
    """ Get a parameter element """
    for part_el in root.findall('wiz_part'):
        if part_el.attrib.get('part') == part:
            # found <wiz_part part="<part>" ...>
            break
    else:
        raise RuntimeError(f'<wiz_part> tag with path="{part}" not found')

    for i in part_el.findall('parameter'):
        if i.attrib.get('name') == parameter_name:
            return i  # <parameter name="<parameter_name>">

    raise RuntimeError(f'<parameter> with name="{parameter_name}" not found')


def make_presets(template, input_data, output=None):
    """
    Compile the presets template

    template:
        Path to the template file
    input_data:
        Path to the mzml formatted mzmine input data.
    output:
        Path to output file.  Print to stdout if None.
    """
    with open('mzmlfiles.txt', 'w') as ofile:
        for i in input_data:
            ofile.write(f'{Path(i).resolve()}\n')

    # root expected to be a <wizard> tag
    root = defusedxml.ElementTree.fromstring(Path(template).read_text())

    # 1. Insert file listing into <paramter name='File names'> tag
    file_names = get_el(root, 'DATA_IMPORT', 'File names')

    # get clean tag, need to iterate over copy of list so remove() won't miss
    # anything
    for i in file_names[:]:
        file_names.remove(i)

    for i in input_data:
        file_tag = SubElement(file_names, 'file')
        file_tag.text = str(i)

    # 2. insert export path
    el = get_el(root, 'WORKFLOW', 'Export path')
    for i in el[:]:
        el.remove(i)

    if output is None:
        export_path = './mzmine.out/X'
    else:
        export_path = Path(output).parent / 'mzmine.out' / 'X'
    current_file = SubElement(el, 'current_file')
    current_file.text = str(export_path)
    last_file = SubElement(el, 'last_file')
    last_file.text = str(export_path)

    # Write out presets file
    if output is None:
        print(ElementTree.tostring(root))
    else:
        tree = ElementTree(root)
        indent(tree, space=' ' * 4, level=0)
        with open(output, 'w') as ofile:
            tree.write(ofile, encoding='unicode', xml_declaration=True)
            ofile.write('\n')
