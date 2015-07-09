"""A setuptools based setup module.

"""

from setuptools import setup, find_packages

setup(
        name="sequtils",
        version="0.1.0",
        description=("Utility scripts and libraries for working with "
                     "sequence data"),
        url="https://github.com/bsmith89/compbio-scripts",
        author="Byron J. Smith",
        author_email="me@byronjsmith.com",
        license="MIT",
        packages=find_packages(),
        install_requires=["biopython", "ipython[notebook]"],
        entry_points={
            'console_scripts': [
                'codonalign=sequtils.scripts.codonalign:main',
                'convert=sequtils.scripts.convert:main',
                'complement=sequtils.scripts.complement:main',
                'drop_seqs=sequtils.scripts.drop_seqs:main',
                'fetch_seqs=sequtils.scripts.fetch_seqs:main',
                'ipynb_output_filter=sequtils.scripts.ipynb_output_filter:main',
                'ls_ids=sequtils.scripts.ls_ids:main',
                'rename_seqs=sequtils.scripts.rename_seqs:main',
                'translate=sequtils.scripts.translate:main',
                'unalign=sequtils.scripts.unalign:main',
                ],
            },
        )
