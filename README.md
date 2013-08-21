VirMet
------
A set of tools for viral metagenomics. This relies on locally installed tools
and resources.

[`prindefr.sh`](https://github.com/ozagordi/VirMet/blob/master/prindefr.sh) is used to

- trim and clean with [`seqtk`](https://github.com/lh3/seqtk);
- filter with [`prinseq`](http://prinseq.sourceforge.net);
- decontaminate with [`deconseq`](http://deconseq.sourceforge.net);
- hit viral database with [`blast`](http://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/);
- list the matching organisms with [`tax_orgs.py`](https://github.com/ozagordi/VirMet/blob/master/tax_orgs.py).
