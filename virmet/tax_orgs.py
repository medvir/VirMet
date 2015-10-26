#!/opt/python2.7/bin/python2.7
import sys
import pandas as pd

# taxonomy file with only sequences from rins_viral_database
#taxo_file = '/data/databases/rins_taxonomy.csv'
#taxonomy = pd.read_csv(taxo_file, index_col='GI')

# blastn hits to the viral database
try:
    hit_file = sys.argv[1]
except KeyError:
    sys.exit('usage: %s hit_file' % sys.argv[0])
hits = pd.read_csv(hit_file, index_col='qseqid',# delim_whitespace=True)
                   delimiter="\t")
print 'Hits found: %d' % hits.shape[0]
# select according to identity and coverage
good_hits = hits[(hits.pident > 75.) & (hits.qcovs > 75.)]
# define a column with genbank id only
#good_hits['sid'] = [int(s.split('_')[1]) for s in good_hits.sseqid]
#good_hits = good_hits.join(taxonomy, on='sid')
matched_reads = good_hits.shape[0]
print 'Hits passing coverage and identity filter: %d' % matched_reads
org_count = good_hits.groupby('sscinames').size()
org_count.order(ascending=False).to_csv('orgs_list.csv', header=True)
