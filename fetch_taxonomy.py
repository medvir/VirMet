#!/usr/bin/env python3.4

'''
Taxonomy file was obtained with the following query
esearch -db protein -query "refseq [filter] AND txid10239 [orgn] NOT txid[131567]" |
efetch -format docsum |
xtract -pattern DocumentSummary -element Gi TaxId Organism > taxonomy.tsv

What follows takes this taxonomy file and a blast result file (tab separated,
outfmt 6) and returns the list of organisms
'''

import sys
import pandas as pd
import argparse

# default action is 'store'
parser = argparse.ArgumentParser(description='Lists organisms in a tsv blast output file',
                                 epilog='Input are mandatory')
parser.add_argument('-t', '--taxonomy', dest='taxonomy',
                    help='tsv taxonomy file with: gi taxid name',
                    default='/data/databases/viral_protein/taxonomy.tsv')
parser.add_argument('-r', '--results', dest='results',
                    help='blast output file with outfmt 6 options')
args_l = parser.parse_args()

taxonomy = pd.read_csv(args_l.taxonomy, names=['gi', 'taxid', 'name'],
                       delimiter='\t', index_col='gi')

# blastn hits to the viral database
hits = pd.read_csv(args_l.results, index_col='qseqid', delimiter="\t",
                   names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                   'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
                   'bitscore'])
print('Hits found: %d' % hits.shape[0])

# select according to evalue
good_hits = hits[(hits.evalue < 1E-7)]

print(good_hits.head())
# define a column with genbank id only and join with taxonomy
good_hits['sid'] = [int(s.split('|')[1]) for s in good_hits.sseqid]
good_hits = good_hits.join(taxonomy, on='sid')
matched_reads = good_hits.shape[0]
print('Hits passing coverage and identity filter: %d' % matched_reads)
org_count = good_hits.groupby('name').size()
org_count.order(ascending=False).to_csv('orgs_list.csv', header=True)