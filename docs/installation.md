### Installation
VirMet relies on a number of third-party tools used to access databases, trim,
convert, filter and map reads. One can refer to the files
[`.travis.yml`](../.travis.yml)
and [`install-dependencies.sh`](../install-dependencies.sh) for details or
further down in this README.
The dependencies are:

- bwa
- samtools 1.3
- tabix
- seqtk
- prinseq-lite
- edirect command line tools
- blast+ 2.3.0
- python (3.x, it's 2016...) with pandas and Biopython
