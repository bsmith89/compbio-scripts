#! /usr/bin/env python
"""From http://stackoverflow.com/a/20844506/1951857

To use::

    # Add the script to your path
    chmod +x path/to/this/ipynb_output_filter.py
    echo "*.ipynb    filter=dropoutput_ipynb" >> ~/.gitattributes
    git config --global core.attributesfile ~/.gitattributes
    git config --global filter.dropoutput_ipynb.clean ipynb_output_filter.py
    git config --global filter.dropoutput_ipynb.smudge cat

When you run ``git status`` you'll see changes not yet staged, but
diff-ing, committing, etc. should all ignore the output/prompt number
portions of the notebook.

You may find that ``git add *.ipynb`` cleans up your status output without
changing the content of the staging area.

"""

import sys
from IPython.nbformat.current import read, write

def clean(in_handle, out_handle):
    json_in = read(in_handle, 'json')

    if "signature" in json_in.metadata:
        json_in.metadata.pop("signature")
    for sheet in json_in.worksheets:
        for cell in sheet.cells:
            if "outputs" in cell:
                cell.outputs = []
            if "prompt_number" in cell:
                cell.pop("prompt_number")

    write(json_in, out_handle, 'json')


if __name__ == '__main__':
    clean(sys.stdin, sys.stdout)
