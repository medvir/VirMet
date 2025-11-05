#!/usr/bin/env python3
"""index sequences with blast or bwa"""

import datetime
import logging
import multiprocessing as mp
import os

from virmet.common import run_child, n_proc


def single_bwa_index(index_params):
    """run a single bwa indexing job"""
    in_fasta, index_prefix = index_params
    size_indx = int(os.path.getsize(in_fasta) / 8)
    cml = "bwa index -b %d -p %s %s &> %s_bwa_index.log" % (
        size_indx,
        index_prefix,
        in_fasta,
        index_prefix,
    )
    run_child(cml)
    return "index %s done" % index_prefix


def single_samtols(index_params):
    """run a single samtools faidx job"""
    in_fasta, index_prefix = index_params
    cml = "samtools faidx %s" % in_fasta
    run_child(cml)
    return "samtools %s done" % index_prefix


def main(args):
    """only function doing all the indexing"""
    DB_DIR = os.path.expandvars(args.dbdir)
    logging.info("now in index")
    logging.info("Database real path: %s" % os.path.realpath(DB_DIR))
    if args.viral == "n":
        target_dir = os.path.join(DB_DIR, "viral_nuccore")
        dt = datetime.date.today().isoformat()
        cml = 'makeblastdb -in {pth_dir}/viral_database.fasta \
        -dbtype nucl -hash_index \
        -title "Viral database indexed {dt_info}" \
        -out {pth_dir}/viral_db \
        -logfile {pth_dir}/blast.log \
        -parse_seqids -taxid_map {pth_dir}/viral_accn_taxid.dmp'.format(
            pth_dir=target_dir, dt_info=dt
        )
        run_child(cml)

    if args.viral == "p":
        target_dir = os.path.join(DB_DIR, "viral_protein")
        dt = datetime.date.today().isoformat()
        cml = 'makeblastdb -in {pth_dir}/viral_database.fasta \
        -dbtype prot -hash_index \
        -title "Viral database indexed {dt_info}" \
        -out {pth_dir}/viral_db \
        -logfile {pth_dir}/blast.log \
        -parse_seqids -taxid_map {pth_dir}/viral_accn_taxid.dmp'.format(
            pth_dir=target_dir, dt_info=dt
        )
        run_child(cml)

    if args.bact_fungal:
        idx_dir = os.path.join(DB_DIR, "bact_fungi")
        run_child(f"k2 build --threads {n_proc} --db {idx_dir}")

    index_pairs = []  # holds (fasta, index) tuples to run in parallel
    if args.human:
        bwa_dir = os.path.join(DB_DIR, "human", "bwa")
        try:
            os.mkdir(bwa_dir)
        except FileExistsError as err:
            logging.warning("FileExistsError: %s" % err)
        fasta_file = os.path.join(DB_DIR, "human", "fasta", "GRCh38.fasta.gz")
        index_prefix = os.path.join(bwa_dir, "humanGRCh38")
        index_pairs.append((fasta_file, index_prefix))

    if args.bovine:
        bwa_dir = os.path.join(DB_DIR, "bovine", "bwa")
        try:
            os.mkdir(bwa_dir)
        except FileExistsError as err:
            logging.warning("FileExistsError: %s" % err)
        fasta_file = os.path.join(
            DB_DIR,
            "bovine",
            "fasta",
            "ref_Bos_taurus.fasta.gz",
        )
        index_prefix = os.path.join(bwa_dir, "bt_ref")
        index_pairs.append((fasta_file, index_prefix))

    # Indexing
    pool_idx = mp.Pool(processes=n_proc)
    results_idx = pool_idx.map(single_bwa_index, index_pairs)
    for r in results_idx:
        logging.info(r)
    # Run in parallel
    pool = mp.Pool(processes=n_proc)
    results = pool.map(single_samtols, index_pairs)
    for r in results:
        logging.info(r)
