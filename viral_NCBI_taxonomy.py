#!/usr/bin/env python
# This was used to extract from gi_taxis (~350 million sequences) only those
# appearing in the rins database that we are using.

import sys
from Bio import SeqIO

# input files
vir_db = '/data/databases/rins_viral_database.fasta'
ncbi_tax_db = '/data/databases/gi_taxid_nucl.dmp'
names_tax = '/data/databases/names_scientific.dmp'
# output file
rins_tax = '/data/databases/rins_taxonomy.csv'

# dictionary gene_id: fasta_id of the rins database sequences
print >> sys.stderr, 'Reading %s' % vir_db
vir_gi = {s.id.split('_')[1]:s.id for s in SeqIO.parse(vir_db, 'fasta')}


# this puts gene_id, tax_id, sequence_id
print >> sys.stderr, 'Reading %s' % ncbi_tax_db
glob_dic = {}
with open(ncbi_tax_db) as infile:
    for line in infile:
        gi, tax = line.split()
        if gi in vir_gi:
            glob_dic[gi] = [gi, tax, vir_gi[gi].replace(',', '')]
            
# this takes the organism name
print >> sys.stderr, 'Reading %s' % names_tax
names_dic = {}
with open(names_tax) as infile:
    for line in infile:
        tid, name = line.split('|')[0].strip(), line.split('|')[1].strip().replace(',', '')
        names_dic[tid] = name

print 'GI,taxid,long_name,scientific_name'

for k, v in glob_dic.items():
    tid = v[1]
    v.append(names_dic[tid])
    print ','.join(v)
