#!/usr/bin/env python
import sys
import pandas as pd

# taxonomy file with only sequences from rins_viral_database
taxo_file = '/data/databases/rins_taxonomy.csv'
taxonomy = pd.read_csv(taxo_file, index_col='GI')

# blastn hits to the viral database
try:
    hit_file = sys.argv[1]
except KeyError:
    sys.exit('usage: %s hit_file' % sys.argv[0])
hits = pd.read_csv(hit_file, index_col='qseqid', delim_whitespace=True)
hits['pcov'] = 100 * (hits.qend + 1.0 - hits.qstart) / hits.length

# select according to identity and coverage
good_hits = hits[(hits.pident > 75.) & (hits.pcov > 75.)]
# define a column with genbank id only
good_hits['sid'] = [int(s.split('_')[1]) for s in good_hits.sseqid]
good_hits = good_hits.join(taxonomy, on='sid')
matched_reads = good_hits.shape[0]
org_count = good_hits.groupby('scientific_name').size()
org_count.to_csv('orgs_list.csv', header=True)
