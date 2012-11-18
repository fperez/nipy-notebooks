#!/usr/bin/env python
from __future__ import print_function

DESCRIP = 'Modified diff between two IPython notebook files'
EPILOG = \
"""
Strips code output and code line numbers from notebooks `filename1` and
`filename2` before running unified diff on the resulting JSON strings.
"""
import io

from argparse import ArgumentParser, RawDescriptionHelpFormatter

import difflib

from IPython.nbformat import current


def remove_outputs(nb):
    """remove the outputs from a notebook"""
    for ws in nb.worksheets:
        for cell in ws.cells:
            if cell.cell_type == 'code':
                cell.outputs = []


def remove_prompt_numbers(nb):
    """remove the prompt numbers from a notebook"""
    for ws in nb.worksheets:
        for cell in ws.cells:
            if cell.cell_type == 'code':
                cell.prompt_number = 0


def main():
    parser = ArgumentParser(description=DESCRIP,
                            epilog=EPILOG,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('filename1', type=str,
                        help='notebook filename')
    parser.add_argument('filename2', type=str,
                        help='notebook filename')
    args = parser.parse_args()
    nb_strs = []
    for fname in (args.filename1, args.filename2):
        with io.open(fname, 'r') as f:
            nb = current.read(f, 'json')
        remove_outputs(nb)
        remove_prompt_numbers(nb)
        nb_strs.append(current.writes(nb, 'json').split('\n'))
    diff_gen = difflib.unified_diff(nb_strs[0], nb_strs[1])
    print('\n'.join(diff_gen))


if __name__ == '__main__':
    main()
