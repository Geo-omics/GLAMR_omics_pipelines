""" mzmine-related things """
from pathlib import Path


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
    template = Path(template).read_text()
    
