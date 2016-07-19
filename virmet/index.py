#!/usr/bin/env python3
'''index sequences with blast or STAR
'''

import sys
import glob
import os.path
import logging
import datetime
from shutil import rmtree
import multiprocessing as mp
from virmet.common import run_child, DB_DIR


def single_bwa_index(index_params):
    '''run a single bwa indexing job'''
    in_fasta, index_prefix = index_params
    cml = 'index -p %s %s &> %s_bwa_index.log' % (index_prefix, in_fasta, index_prefix)
    run_child('bwa', cml, exe='/bin/bash')
    return 'index %s done' % index_prefix


def main(args):
    '''only function doing all the indexing'''
    logging.info('now in index')

    if args.viral == 'n':
        target_dir = os.path.join(DB_DIR, 'viral_nuccore')
        os.chdir(target_dir)
        dt = datetime.date.today().isoformat()
        cml = "-in viral_database.fasta -dbtype nucl -hash_index \
        -title \"Viral database indexed {}\" \
        -out viral_db \
        -logfile blast.log -parse_seqids -taxid_map viral_gi_taxid.dmp".format(dt)
        run_child('makeblastdb', cml)

    if args.viral == 'p':
        target_dir = os.path.join(DB_DIR, 'viral_protein')
        os.chdir(target_dir)
        dt = datetime.date.today().isoformat()
        cml = "-in viral_database.fasta -dbtype prot -hash_index \
        -title \"Viral database indexed {}\" \
        -out viral_db \
        -logfile blast.log -parse_seqids -taxid_map viral_gi_taxid.dmp".format(dt)
        run_child('makeblastdb', cml)

    index_pairs = []  # holds (fasta, index) tuples to run in parallel
    if args.bact:
        bwa_dir = os.path.join(DB_DIR, 'bacteria', 'bwa')
        try:
            os.mkdir(bwa_dir)
        except FileExistsError as err:
            logging.warning('FileExistsError: %s' % err)
        for i in [1, 2, 3]:
            fasta_file = os.path.join(DB_DIR, 'bacteria', 'fasta', 'bact%d.fasta.gz' % i)
            index_prefix = os.path.join(bwa_dir, 'bact%d' % i)
            index_pairs.append((fasta_file, index_prefix))

    if args.human:
        bwa_dir = os.path.join(DB_DIR, 'human', 'bwa')
        try:
            os.mkdir(bwa_dir)
        except FileExistsError as err:
            logging.warning('FileExistsError: %s' % err)
        fasta_file = os.path.join(DB_DIR, 'human', 'fasta', 'GRCh38.fasta.gz')
        index_prefix = os.path.join(bwa_dir, 'humanGRCh38')
        index_pairs.append((fasta_file, index_prefix))

    if args.fungal:
        bwa_dir = os.path.join(DB_DIR, 'fungi', 'bwa')
        try:
            os.mkdir(bwa_dir)
        except FileExistsError as err:
            logging.warning('FileExistsError: %s' % err)
        fasta_file = os.path.join(DB_DIR, 'fungi', 'fasta', 'fungi1.fasta.gz')
        index_prefix = os.path.join(bwa_dir, 'fungi1')
        index_pairs.append((fasta_file, index_prefix))

    if args.bovine:
        bwa_dir = os.path.join(DB_DIR, 'bovine', 'bwa')
        try:
            os.mkdir(bwa_dir)
        except FileExistsError as err:
            logging.warning('FileExistsError: %s' % err)
        fasta_file = os.path.join(DB_DIR, 'bovine', 'fasta', 'bt_ref_Bos_taurus_UMD_3.1.1.fasta.gz')
        index_prefix = os.path.join(bwa_dir, 'bt_ref')
        index_pairs.append((fasta_file, index_prefix))

    # run in parallel
    # TODO: use single_process
    pool = mp.Pool()
    results = pool.map(single_bwa_index, index_pairs)
    for r in results:
        logging.info(r)

    # TODO parallelize this too
    for fasta_file, prefix in index_pairs:
        run_child('samtools', 'faidx %s' % fasta_file)
