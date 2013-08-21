#!/usr/bin/env python
import sys
import warnings
import pandas as pd

# taxonomy file with only sequences from rins_viral_database
taxo_file = '/data/databases/rins_taxonomy.csv'
taxids = pd.read_csv(taxo_file)

# Genbank taxonomy file names: delimiter also at the end of each line,
# so index_col=False or usecols
# squeezed to remove tabs with `tr -d "\t"` and grepped on scientific name
names_file = '/data/databases/scientific_names.dmp'
names_names = ['taxid', 'name', 'unique_name', 'name_class']
names = pd.read_csv(names_file, delimiter="|", header=0, names=names_names,
index_col=0, usecols=[0, 1, 2, 3])

# Genbank taxonomy file nodes
nodes_file = '/data/databases/squeezed_nodes.dmp'
node_names = ['tax_id',	'parent_tax_id', 'rank'] #  'embl_code', 'division_id',
# 'inherited_div_flag', 'genetic_code_id','inherited_GC_flag',
# 'mitochondrial_genetic_code_id', 'inherited_MGC_flag',
# 'GenBank_hidden_flag', 'hidden_subtree_root_flag', 'comments']
nodes = pd.read_csv(nodes_file, delimiter="|", index_col=0, names=node_names,
usecols=[0, 1, 2])

def get_upper_level(taxids, current_level):
    if len(taxids) != len(set(taxids)):
         warnings.warn('Not unique taxa')
    print 'Asking for %d taxa, %d different taxa' % \
        (len(taxids), len(set(taxids)))
    ids = pd.Series(taxids)
    x = nodes.loc[ids]
    x = x[x['rank'] == current_level]['parent_tax_id']
    print '%d were found at the requested %s level' % (x.shape[0], current_level)
    print '%d are unique' % x.nunique()
    y = names.loc[x]
    y_count = y.groupby('name').size()
    y_count.to_csv('above_%s_list.csv' % current_level, header=True)
    print '\n'
    return x
    

# Extracts only nodes in rins_taxonomy and corresponding to species
rins_ti = taxids['taxid']
rins_genera = get_upper_level(rins_ti, 'species')
rins_families = get_upper_level(rins_genera, 'genus')
rins_orders = get_upper_level(rins_families, 'family')
