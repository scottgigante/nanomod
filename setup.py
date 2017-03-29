import os
from setuptools import setup

version_py = os.path.join(os.path.dirname(__file__), 'nanomod', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"','').strip()
print version

def read(fname):
    """
    Utility function to read the README file.
    Used for the long_description.  It's nice, because now 1) we have a top level
    README file and 2) it's easier to type in the README file than to put a raw
    string in below ...

    :param fname: string Path to file

    :returns: string File contents
    """
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
  name = 'nanomod',
  packages = ['nanomod'],
  package_dir={'nanomod': "nanomod"},
  version=version,
  install_requires=['numpy>=1.7', 'h5py>2.2.0', 'biopython', 'pysam', 'poretools'],
  requires=['python (>=2.7, <3.0)'],
  description = 'A tool designed to call base modifications on Oxford Nanopore Technologies\' MinION sequencing data.',
  long_description=read('README.md'),
  author = 'Scott Gigante',
  author_email = 'scottgigante@gmail.com',
  url = 'https://github.com/scottgigante/nanomod', # use the URL to the github repo
  keywords = ['nanopore', 'methylation', 'epigenomics', 'hdf5', 'fast5', 'oxford', 'minion', 'basecalling'], # arbitrary keywords
  classifiers = [
        "Development Status :: 3 - Alpha",
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ],
  entry_points = {
        'console_scripts': ['nanomod = nanomod.__main__:run'],
    },
)
