#!/opt/python2.7/bin/python2.7
'''Write viral and undetermined reads into different files.
Usage: sift_reads.py all_reads.fastq blast_results.tsv'''

import sys

try:
    args = sys.argv[1:]
    all_reads, blast_res = args
except ValueError:
    sys.exit('usage: %s all_reads.fastq blast_results.tsv' % sys.argv[0])

import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# selects reads with coverage and identity higher than 75
df = pd.read_csv('unique.tsv', sep='\t')
viral_ids = set(df[(df.qcovs > 75) & (df.pident > 75)].qseqid)

viral_c = 0
undet_c = 0

all_handle = open(all_reads)
undet_handle = open('undetermined_reads.fastq', 'w')
viral_handle = open('viral_reads.fastq', 'w')

# Using FastqGeneralIterator allows fast performance
for title, seq, qual in FastqGeneralIterator(all_handle):
    if title.split()[0] not in viral_ids:
        undet_c += 1
        undet_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        if undet_c % 100000 == 0:
            print >> sys.stderr, 'written %d undet reads' % undet_c
    else:
        viral_c += 1
        viral_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        if viral_c % 10000 == 0:
            print >> sys.stderr, 'written %d viral reads' % viral_c


undet_handle.close()
viral_handle.close()
print >> sys.stderr, 'written %d undet reads' % undet_c
print >> sys.stderr, 'written %d viral reads' % viral_c
