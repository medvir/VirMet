VirMet
------
A set of tools for viral metagenomics. This relies on locally installed tools
and resources.

[`hunter.sh`](https://github.com/ozagordi/VirMet/blob/master/hunter.sh) is used to

- trim and clean with [`seqtk`](https://github.com/lh3/seqtk);
- filter with [`prinseq`](http://prinseq.sourceforge.net);
- decontaminate with [`STAR`](https://code.google.com/p/rna-star/);
- hit viral database with [`blast`](http://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/);
- list the matching organisms with [`tax_orgs.py`](https://github.com/ozagordi/VirMet/blob/master/tax_orgs.py).

