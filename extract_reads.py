#!/opt/python2.7/bin/python2.7
import sys
from Bio import SeqIO

ids_file = sys.argv[1]
fasta_file = sys.argv[2]

ids = set([l.strip() for l in open(ids_file)])

seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

matching = [seq_dict[k] for k in ids]

SeqIO.write(matching, 'matching.fastq', 'fastq')

