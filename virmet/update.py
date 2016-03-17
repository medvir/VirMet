#!/usr/bin/env python3.4
''' updates viral and bacterial database based on a new query to ncbi and a
manually added list of GIs
'''
import os
import sys
import time
import logging
import subprocess
import urllib.request
import pandas as pd
from virmet.common import run_child, viral_query, bacterial_query, get_gids, \
download_genomes, DB_DIR


def bactupdate(picked=None):
    bact_dir = os.path.join(DB_DIR, 'bacteria')
    os.chdir(bact_dir)
    os.rename('bacteria_refseq_info.tsv', 'old_bacteria_refseq_info.tsv')
    old_urls = bacterial_query(download=False,
                               info_file='old_bacteria_refseq_info.tsv')
    logging.info('%d bacterial assemblies present in refseq' % (len(old_urls)))
    new_urls = bacterial_query(download=True,
                               info_file='bacteria_refseq_info.tsv')
    logging.info('%d bacterial assemblies now in refseq' % (len(new_urls)))
    to_add = set(new_urls) - set(old_urls)
    for t in to_add:
        logging.info('genome from %s will be added' % t)
    download_bact_urls(to_add)
    if picked is None:
        return

    present_ids = itertools.chain.from_iterable([get_ids(f) \
        for f in glob.glob('fasta/bact*gz')])
    picked_ids = [l.strip() for l in open(picked)]
    to_add = set(present_ids) - set(picked_ids)

    for i, gid in enumerate(to_add):
        fileout = 'fasta/bact%d.fasta.gz' % ((i % 3) + 1)
        run_child('efetch',
                  '-db nuccore -id %s -format fasta | gzip - >> %s' % (gid, fileout))
    logging.info('added %d sequences from file %s' % (i, picked))


def virupdate(viral_type, picked=None):
    if viral_type == 'n':
        db_type = 'nuccore'
    elif viral_type == 'p':
        db_type = 'protein'
    viral_dir = os.path.join(DB_DIR, 'viral_%s' % db_type)

    # this query downloads a new viral_seqs_info.tsv and parses the GI
    logging.info('interrogating NCBI again')
    viral_query(viral_type)
    info_file = os.path.join(viral_dir, 'viral_seqs_info.tsv')
    info_seqs = pd.read_csv(info_file, sep='\t',
                            names=['Gi', 'TaxId', 'Caption', 'Slen', 'Organism', 'Title'])
    new_ids = [str(gi) for gi in info_seqs['Gi'].tolist()]
    logging.info('NCBI reports %d sequences' % len(new_ids))

    # read ids already present in fasta file
    fasta_db = os.path.join(viral_dir, 'viral_database.fasta')
    present_ids = get_gids(fasta_db)
    logging.info('fasta file has %d sequences' % len(present_ids))

    # sequences given manually by specifying file with GI
    if picked:
        manual_ids = [l.strip() for l in open(picked)]
        logging.info('%d sequences specified manually' % len(manual_ids))
    else:
        manual_ids = []

    # update fasta: ids to add are union of picked plus those in ncbi minus those present
    ids_to_add = set(manual_ids) | set(new_ids)
    ids_to_add = ids_to_add - set(present_ids)
    if len(ids_to_add) == 0:
        logging.info('no sequences to add to fasta file')
    elif len(ids_to_add) > 200:
        logging.error('cannot add %d sequences, exiting' % len(ids_to_add))
        sys.exit('too many sequences to add: run `virmet fetch` first')
    else:
        logging.info('adding %d sequences to fasta file' % len(ids_to_add))
        s_code = run_child('efetch',
                           '-db %s -id ' % db_type + ','.join(ids_to_add) +
                           ' -format fasta >> %s' % fasta_db)

    # update viral_seqs_info.tsv and taxonomy
    ids_to_add = set(present_ids) | set(manual_ids)
    ids_to_add = ids_to_add - set(new_ids)
    if len(ids_to_add) == 0:
        logging.info('no sequences to add to viral_seqs_info')
    else:
        logging.info('adding %d line(s) to viral_seqs_info.tsv' % len(ids_to_add))
        # loop needed as efetch with format docsum only takes one id at a time
        # (change introduced in edirect 3.30, December 2015)
        # slow, but other solutions seem complicated with edirect
        for ita in ids_to_add:
            cml = '-db %s -id %s' % (db_type, ita)
            cml = cml + ' -format docsum | xtract -pattern DocumentSummary \
            -element Gi TaxId Caption Slen Organism Title >> %s' % info_file
            run_child('efetch', cml)

    logging.info('updating taxonomy')
    s_code = run_child('cut', '-f 1,2 %s > %s' % (info_file, os.path.join(viral_dir, 'viral_gi_taxid.dmp')))

    # perform tests
    gids_1 = set(get_gids('viral_database.fasta'))
    gids_2 = set([l.split()[0] for l in open('viral_gi_taxid.dmp')])
    assert gids_1 == gids_2, 'taxonomy/viral_seqs_info not matching with fasta'

def main(args):
    logging.info('now in update_db')
    if args.viral and args.bact:
        logging.error('update either viral or bacterial in a single call')
        sys.exit('update either viral or bacterial in a single call')
    virupdate(args.viral, args.picked)
    if args.bact:
        bactupdate(args.picked)
