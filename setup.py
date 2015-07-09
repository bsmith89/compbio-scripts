"""A setuptools based setup module.

"""

from setuptools import setup

setup(
        name="sequtils",
        version="0.1.0",
        description=("Utility scripts and libraries for working with "
                     "sequence data"),
        url="https://github.com/bsmith89/compbio-scripts",
        author="Byron J. Smith",
        author_email="me@byronjsmith.com",
        license="MIT",
        packages=["sequtils"],
        install_requires=["biopython", "ipython[notebook]"],
        entry_points={
            'console_scripts': [
                'codonalign=codonalign:main',
                'convert=convert:main',
                'complement=complement:main',
                'convert=convert:main',
                'drop_seqs=drop_seqs:main',
                'fetch_seqs=fetch_seqs:main',
                'ipynb_output_filter=ipynb_output_filer:main',
                'ls_ids=ls_ids:main',
                'rename_seqs=rename_seqs:main',
                'translate=translate:main',
                'unalign=unalign:main',
                ],
            },
        )
