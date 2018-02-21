#!/usr/bin/env python3
"""Update viral and bacterial database based on a new query to ncbi and a manually added list of GIs


"""
import os
import sys
import logging
import pandas as pd
from virmet.common import run_child, viral_query, bact_fung_query, get_accs, download_genomes, DB_DIR


def bact_fung_update(query_type=None, picked=None):
    """
    """
    import glob
    import itertools

    cont_dir = os.path.join(DB_DIR, query_type)
    os.chdir(cont_dir)
    logging.info('updating %s, now in %s', query_type, cont_dir)
    # read old info
    os.rename('%s_refseq_info.tsv' % query_type, 'old_%s_refseq_info.tsv' % query_type)
    old_urls = bact_fung_query(query_type=query_type, download=False, info_file='old_%s_refseq_info.tsv' % query_type)

    logging.info('%d assemblies were present in refseq', len(old_urls))
    # download new info
    new_urls = bact_fung_query(query_type=query_type, download=True)
    logging.info('%d assemblies are now in refseq', len(new_urls))
    to_add = set(new_urls) - set(old_urls)
    to_add = list(to_add)

    if not to_add:
        logging.info('no new sequences in %s database', query_type)
        print('no new sequences in %s database' % query_type, file=sys.stderr)

    for t in to_add:
        logging.debug('genome from %s will be added', t)
        if query_type == 'bacteria':
            download_genomes(to_add, prefix='tmp', n_files=3)
            for i in [1, 2, 3]:
                run_child('bgzip -c fasta/tmp%d.fasta >> fasta/bact%d.fasta.gz' % (i, i))
                os.remove('fasta/bact%d.fasta.gz' % i)
        elif query_type == 'fungi':
            download_genomes(to_add, prefix='tmp', n_files=1)
            run_child('bgzip -c fasta/tmp1.fasta >> fasta/fungi1.fasta.gz')
            os.remove('fasta/fungi1.fasta.gz')

    if picked is None:
        return

    # present_ids = itertools.chain.from_iterable([get_gids(f) for f in glob.glob('fasta/*.fasta.gz')])
    present_ids = itertools.chain.from_iterable([get_accs(f) for f in glob.glob('fasta/*.fasta.gz')])
    picked_ids = [l.strip() for l in open(picked)]
    to_add = set(present_ids) - set(picked_ids)

    if not to_add:
        logging.info('no new sequence manually added')
        print('no new sequence manually added', file=sys.stderr)

    for i, gid in enumerate(to_add):
        if query_type == 'bacteria':
            fileout = 'fasta/bact%d.fasta.gz' % ((i % 3) + 1)
        elif query_type == 'fungi':
            fileout = 'fasta/fungi%d.fasta.gz' % ((i % 1) + 1)
        run_child('bgzip -c <(efetch -db nuccore -id %s -format fasta) >> %s' % (gid, fileout), exe='/bin/bash')
    logging.info('added %d sequences from file %s', i, picked)
    if query_type == 'bacteria':
        for i in [1, 2, 3]:
            run_child('bgzip -r fasta/bact%d.fasta.gz')
    elif query_type == 'fungi':
        run_child('bgzip -r fasta/fungi1.fasta.gz')


def virupdate(viral_type, picked=None):
    if viral_type == 'n':
        db_type = 'nuccore'
    elif viral_type == 'p':
        db_type = 'protein'
    viral_dir = os.path.join(DB_DIR, 'viral_%s' % db_type)

    # this query downloads a new viral_seqs_info.tsv and parses the GI
    logging.info('interrogating NCBI again')
    os.chdir(viral_dir)
    cml_search = viral_query(viral_type)
    run_child(cml_search)
    efetch_xtract = 'efetch -format docsum < ncbi_search | xtract'
    efetch_xtract += ' -pattern DocumentSummary -element Caption TaxId Slen Organism Title > viral_seqs_info.tsv'
    run_child(efetch_xtract)
    info_file = os.path.join(viral_dir, 'viral_seqs_info.tsv')
    info_seqs = pd.read_csv(info_file, sep='\t', names=['Caption', 'TaxId', 'Slen', 'Organism', 'Title'])
    new_ids = [str(acc) for acc in info_seqs['Caption'].tolist()]
    logging.info('NCBI reports %d sequences', len(new_ids))

    # read ids already present in fasta file
    fasta_db = os.path.join(viral_dir, 'viral_database.fasta')
    present_ids = get_accs(fasta_db)
    logging.info('fasta file has %d sequences', len(present_ids))

    # sequences given manually by specifying file with GI
    if picked:
        manual_ids = [l.strip() for l in open(picked)]
        logging.info('%d sequences specified manually', len(manual_ids))
    else:
        manual_ids = []

    # update fasta: ids to add are union of picked plus those in ncbi minus those present
    ids_to_add = set(manual_ids) | set(new_ids)
    ids_to_add = ids_to_add - set(present_ids)
    if not ids_to_add:
        logging.info('no sequences to add to fasta file')
        print('no sequences to add to fasta file', file=sys.stderr)
    elif len(ids_to_add) > 2000:
        logging.error('cannot add %d sequences, exiting', len(ids_to_add))
        sys.exit('too many sequences to add: run `virmet fetch` first')
    else:
        logging.info('adding %d sequences to fasta file', len(ids_to_add))
        s_code = run_child('efetch -db %s -id ' % db_type + ','.join(ids_to_add) + ' -format fasta >> %s' % fasta_db)
        logging.debug(s_code)

    # update viral_seqs_info.tsv and taxonomy
    ids_to_add = set(present_ids) | set(manual_ids)
    ids_to_add = ids_to_add - set(new_ids)
    if not ids_to_add:
        logging.info('no sequences to add to viral_seqs_info')
        print('no sequences to add to viral_seqs_info', file=sys.stderr)
    else:
        logging.info('adding %d line(s) to viral_seqs_info.tsv', len(ids_to_add))
        # loop needed as efetch with format docsum only takes one id at a time
        # (change introduced in edirect 3.30, December 2015)
        # slow, but other solutions seem complicated with edirect
        for ita in ids_to_add:
            cml = 'efetch -db %s -id %s' % (db_type, ita)
            cml = cml + ' -format docsum | xtract -pattern DocumentSummary \
            -element Caption TaxId Slen Organism Title >> %s' % info_file
            run_child(cml)

    logging.info('updating taxonomy')
    s_code = run_child('cut -f 1,2 %s > %s' % (info_file, os.path.join(viral_dir, 'viral_accn_taxid.dmp')))

    # perform tests
    gids_1 = set(get_accs('viral_database.fasta'))
    gids_2 = set([l.split()[0] for l in open('viral_accn_taxid.dmp')])
    assert gids_1 == gids_2, 'taxonomy/viral_seqs_info not matching with fasta'

def main(args):
    logging.info('now in update_db')
    if bool(args.viral) + args.bact + args.fungal > 1:
        logging.error('update either viral or bacterial or fungal in a single call')
        sys.exit('update either viral or bacterial or fungal in a single call')
    if args.viral:
        virupdate(args.viral, args.picked)
    elif args.bact:
        bact_fung_update(query_type='bacteria', picked=args.picked)
    elif args.fungal:
        bact_fung_update(query_type='fungi', picked=args.picked)
