#!/usr/bin/env python3.4
# !/opt/python2.7/bin/python2.7
'''Search and download complete viral genomes or refseq proteins from NCBI with
edirect tools. Check previously downloaded sequences and manually picked in
file picked_seqs.txt.
'''
import sys
import os
import re
import datetime

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

Entrez_email = 'osvaldo.zagordi@uzh.ch'
c_code = run_child('econtact', '-email \"%s\"' % Entrez_email)

try:
    db_type = sys.argv[1]
except KeyError as ee:
    sys.exit('usage: %s db_type [nuccore, protein]' % sys.argv[0])

if db_type not in ['nuccore', 'protein']:
    sys.exit('db_type: %s, must be nuccore or protein' % db_type)

os.chdir('/data/databases/viral_%s' % db_type)

print('Parsing sequences already present')
p_code = run_child('grep',
                   ' \"^>\" viral_database.fasta | cut -f 2 -d \"|\" > old_ids')
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
# txid = '10239'  # change here for viruses or smaller taxa
txids = {
    '39759': 'Deltaviruses',
    '35237': 'dsDNA viruses, no RNA stage',
    '35325': 'dsRNA viruses',
    '35268': 'retrotranscribingviruses',
    '12877': 'Satellites',
    '29258': 'ssDNA viruses',
    '439488': 'ssRNA viruses',
    '686617': 'unassigned viruses',
    '451344': 'unclassified archaeal viruses',
    '12333': 'unclassified phages',
    '552364': 'unclassified virophages',
    '12429': 'unclassified viruses'
    }
for txid, v in txids.items():
    if db_type == 'nuccore':
        search_text = "txid%s [orgn] " % txid + \
                      "AND \\\"complete genome\\\" [Title] " + \
                      "NOT txid131567 [orgn]"

    elif db_type == 'protein':
        search_text = "refseq [filter] AND txid%s [orgn] " % txid + \
                      "NOT txid131567 [orgn]"

    s_code = run_child('esearch',
                       '-db %s -query \"%s\" > vir_search_%s'
                       % (db_type, search_text, txid))

# get the number of hits per taxa
count = 0
for txid, v in txids.items():
    for l in open('vir_search_%s' % txid):
        if l.lstrip().startswith('<Count'):
            ch = int(re.search('Count>(\d+)</Count', l).group(1))
            print('Found %d hits from %s' % (ch, v))
            count += ch
            break
print('Summing these, NCBI search returned %d hits' % count)


# Now run the search on the whole virus branch
txid = '10239'
if db_type == 'nuccore':
    search_text = "txid%s [orgn] " % txid + \
                  "AND \\\"complete genome\\\" [Title] " + \
                  "NOT txid131567 [orgn]"

elif db_type == 'protein':
    search_text = "refseq [filter] AND txid%s [orgn] " % txid + \
                  "NOT txid131567 [orgn]"

s_code = run_child('esearch',
                   '-db %s -query \"%s\" > vir_search'
                   % (db_type, search_text))

for l in open('vir_search'):
    if l.lstrip().startswith('<Count'):
        ch = int(re.search('Count>(\d+)</Count', l).group(1))
        print('NCBI search on all viruses returned %d hits' % ch)

try:
    os.remove('tmp.dmp')
except OSError:
    pass

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

id_file = 'picked_seqs.txt'
print('Now the manually picked sequences from file:', id_file)

try:
    picked_ids = [line.strip().split()[0] for line in open(id_file)]
except IOError:
    picked_ids = []
ids_count = len(picked_ids)
print ('{} sequences listed in {} file'.format(ids_count, id_file))

# ids to add are the union of picked plus those in ncbi minus those present
to_add = set(picked_ids) | set(viral_ids)
to_add = to_add - set(old_ids)

print('{} sequences must be added'.format(len(to_add)))

if to_add:
    print('adding the missing {}'.format(len(to_add)))
    print('first to fasta')
    s_code = run_child('efetch',
                       '-db %s -id ' % db_type + ','.join(to_add) +
                       ' -format fasta >> viral_database.fasta')
    print('then to taxonomy')
    s_code = run_child('efetch',
                       '-db %s -id ' % db_type + ','.join(to_add) +
                       ' -format docsum | \
                        xtract -pattern DocumentSummary -element Gi TaxId \
                        >> viral_gi_taxid.dmp')
else:
    print('No new sequences to add')
os.remove('vir_search')

# update blast database
cml_str = 'makeblastdb -in viral_database.fasta -dbtype '
if db_type == 'nuccore':
    cml_str += 'nucl '
else:
    cml_str += 'prot '

tl = len(to_add) + len(set(old_ids))
dt = datetime.date.today().isoformat()
cml_str += "-hash_index \
-title \"Viral database {} sequences {}\" \
-out viral_db \
-logfile blast.log -parse_seqids -taxid_map viral_gi_taxid.dmp".format(tl, dt)

run_child(cml_str)
