#!/usr/bin/env python
# !/opt/python2.7/bin/python2.7
'''Search and download complete viral genomes from NCBI with edirect tools.
Check previously downloaded sequences and manually picked in file
picked_seqs.txt.
'''
import sys
import os
import re
from Bio import Entrez

Entrez.email = 'osvaldo.zagordi@uzh.ch'


def run_child(exe_name, arg_string):
    '''use subrocess to run an external program with arguments'''
    import subprocess
    print(exe_name + ' ' + arg_string)
    if not arg_string.startswith(' '):
        arg_string = ' ' + arg_string

    try:
        retcode = subprocess.call(exe_name + arg_string, shell=True)
        if retcode > 0:
            sys.exit("Child %s %s terminated by signal %d" %
                     (exe_name, arg_string, retcode))
    except OSError as ee:
        sys.exit("Execution of %s failed: %s" % (exe_name, ee))

    return retcode

c_code = run_child('econtact', '-email \"%s\"' % Entrez.email)

print('Parsing sequences already present')
p_code = run_child('grep',
                   ' \">\" viral_database.fasta | cut -f 2 -d \"|\" > old_ids')
old_ids = [l.strip() for l in open('old_ids')]
os.remove('old_ids')
assert len(old_ids) == len(set(old_ids)), 'Duplicate in old seqs: solve it'
print('%d sequences found' % len(old_ids))

print("Running NCBI search...")
# Run the search again
# Viruses, Taxonomy ID: 10239
# Human adenovirus A, Taxonomy ID: 129875 (only for testing, 7 hits)
# Mastadenovirus, Taxonomy ID: 10509 (only for testing, 341 hits)
# Cellular organisms, Taxonomy ID: 131567 (to avoid chimeras)
txid = '10509'  # change here for viruses or smaller taxa
search_text = "txid%s [orgn] " % txid + \
              "AND \\\"complete genome\\\" [Title] " + \
              "NOT txid131567 [orgn]"
s_code = run_child('esearch',
                   '-db nucleotide -query \"%s\" > vir_search' % search_text)

# get the number of hits
for l in open('vir_search'):
    if l.lstrip().startswith('<Count'):
        count = int(re.search('Count>(\d+)</Count', l).group(1))
print('NCBI search returned %d hits' % count)

print('Saving Gi <-> TaxId <-> Acc relationship of the search')
s_code = run_child('efetch',
                   '-format docsum < vir_search | ' +
                   'xtract -pattern DocumentSummary ' +
                   '-element Gi TaxId Caption > tmp.dmp')

s_code = run_child('cut', '-f 1,2 tmp.dmp > viral_gi_taxid.dmp')
# save ids for later
viral_ids = [l.split()[0] for l in open('viral_gi_taxid.dmp')]

print('Saving accession numbers')
s_code = run_child('cut', '-f 3 tmp.dmp > viral_search_NCBI_acc.txt')
os.remove('tmp.dmp')

print('Now the manually picked sequences')
id_file = 'picked_seqs.txt'
try:
    picked_ids = [line.strip() for line in open(id_file)]
except IOError:
    picked_ids = []
ids_count = len(picked_ids)
print ('%d sequences listed in %s file' % (ids_count, id_file))

to_add = set(picked_ids) | set(viral_ids)
to_add = to_add - set(old_ids)

print('%d sequences must be added' % len(to_add))

if to_add:
    print('adding the missing %d' % len(to_add))
    print('first to fasta')
    s_code = run_child('efetch',
                       '-db nuccore -id ' + ','.join(to_add) +
                       ' -format fasta >> viral_database.fasta')
    print('then to taxonomy')
    s_code = run_child('efetch',
                       '-db nuccore -id ' + ','.join(to_add) +
                       ' -format docsum | \
                        xtract -pattern DocumentSummary -element Gi TaxId \
                        >> viral_gi_taxid.dmp')
else:
    print('No new sequences to add')
os.remove('vir_search')
