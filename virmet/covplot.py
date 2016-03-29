#!/usr/bin/env python3.4

'''Runs on all samples of a MiSeq run or on a single fastq file'''
import os
import sys
import glob
import logging
import pandas as pd
from virmet.common import run_child, single_process, DB_DIR



def main(args):
    ''''''
    outdir = args.outdir
    organism = args.organism
    try:
        os.chdir(outdir)
    except FileNotFoundError:
        sys.exit('Check the path of output directory')
    for f in ['orgs_list.tsv', 'unique.tsv']:
        assert os.path.exists(f)
    orgs_list = pd.read_csv('orgs_list.tsv', sep='\t', header=0)
    criterion = orgs_list['organism'].map(lambda x: x.startswith(organism))
    matching_orgs = orgs_list[criterion]
    print(str(matching_orgs.iloc[0].organism))
    # blast results
    unique = pd.read_csv('unique.tsv', sep='\t', header=0)
    criterion = unique['sscinames'].map(lambda x: x.startswith(organism))
    matching_reads = unique[criterion]
