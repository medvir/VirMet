# Installation Guide

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/virmet/README.html)

VirMet is available through [Bioconda](https://bioconda.github.io), a channel
for the [conda](http://conda.pydata.org/docs/intro.html) package manager. Once
conda is [installed](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) 
and the channels are set up, you can install the package with all its dependencies by doing:

`conda install virmet`

## Dependencies

If you don't want to install VirMet the easy way with `conda install virmet` 
(which is highly recommended), try the classic:

`pip install virmet` 

This should also work.
However, VirMet relies on several third-party tools used to access databases, trim,
convert, filter and map reads. These will need to be installed manually.

Python and python dependencies:

- python >=3.9
- biopython (already installed with the `pip` command)
- numpy (already installed with the `pip` command)
- pandas (already installed with the `pip` command)
- requests (already installed with the `pip` command)

Third-party tools:

- blast >=2.3
- bwa
- fastp
- htslib
- kraken2 >=2.1.6
- ncbi-datasets-cli
- r-ggplot2
- samtools >=1.3
- seqkit
- seqtk

Please, note that these tools are all automatically available if VirMet is
installed with `conda install virmet`.