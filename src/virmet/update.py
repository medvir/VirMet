#!/usr/bin/env python3
"""Update viral and bacterial database based on a new query to ncbi and a manually added list of GIs"""

import glob
import itertools
import logging
import os
import sys
import warnings
from collections import Counter

import pandas as pd

from virmet.common import (
    n_proc,
    get_accs,
    run_child,
    viral_query,
)


def virupdate(DB_DIR, viral_type, picked=None, update_min_date=None):
    if viral_type == "n":
        db_type = "nuccore"
    elif viral_type == "p":
        db_type = "protein"
    else:
        raise ValueError(f'Invalid db_type value: "{db_type}".')
    viral_dir = os.path.join(DB_DIR, f"viral_{db_type}")

    # this query downloads a new viral_seqs_info.tsv and parses the GI
    logging.info("interrogating NCBI again")
    os.chdir(viral_dir)
    cml_search = viral_query(viral_type, update_min_date)
    run_child(cml_search)
    efetch_xtract = "efetch -format docsum < ncbi_search | xtract"
    efetch_xtract += " -pattern DocumentSummary -element Caption TaxId Slen Organism Title AccessionVersion > viral_seqs_info.tsv"
    run_child(efetch_xtract)
    info_file = os.path.join(viral_dir, "viral_seqs_info.tsv")
    info_seqs = pd.read_csv(
        info_file,
        sep="\t",
        names=[
            "Caption",
            "TaxId",
            "Slen",
            "Organism",
            "Title",
            "AccessionVersion",
        ],
    )
    new_ids = [str(acc) for acc in info_seqs["Caption"].tolist()]
    logging.info("NCBI reports %d sequences" % len(new_ids))

    # read ids already present in fasta file
    fasta_db = os.path.join(viral_dir, "viral_database.fasta")
    present_ids = get_accs(fasta_db)
    logging.info("fasta file has %d sequences" % len(present_ids))

    # sequences given manually by specifying file with GI
    if picked:
        manual_ids = [l.strip() for l in open(picked)]
        logging.info("%d sequences specified manually" % len(manual_ids))
    else:
        manual_ids = []

    # update fasta: ids to add are union of picked plus those in ncbi minus those present
    ids_to_add = set(manual_ids) | set(new_ids)
    ids_to_add = ids_to_add - set(present_ids)
    if not ids_to_add:
        logging.info("no sequences to add to fasta file")
        print("no sequences to add to fasta file", file=sys.stderr)
    elif len(ids_to_add) > 2000:
        logging.error("cannot add %d sequences, exiting" % len(ids_to_add))
        sys.exit("too many sequences to add: run `virmet fetch` first")
    else:
        logging.info("adding %d sequences to fasta file" % len(ids_to_add))
        s_code = run_child(
            "efetch -db %s -id " % db_type
            + ",".join(ids_to_add)
            + " -format fasta >> %s" % fasta_db
        )
        logging.debug(s_code)

    # update viral_seqs_info.tsv and taxonomy
    ids_to_add = set(present_ids) | set(manual_ids)
    ids_to_add = ids_to_add - set(new_ids)
    if not ids_to_add:
        logging.info("no sequences to add to viral_seqs_info")
        print("no sequences to add to viral_seqs_info", file=sys.stderr)
    else:
        logging.info(
            "adding %d line(s) to viral_seqs_info.tsv" % len(ids_to_add)
        )
        # loop needed as efetch with format docsum only takes one id at a time
        # (change introduced in edirect 3.30, December 2015)
        # slow, but other solutions seem complicated with edirect
        for ita in ids_to_add:
            cml = "efetch -db %s -id %s" % (db_type, ita)
            cml = (
                cml
                + " -format docsum | xtract -pattern DocumentSummary \
            -element Caption TaxId Slen Organism Title >> %s"
                % info_file
            )
            run_child(cml)

    logging.info("updating taxonomy")
    s_code = run_child(
        "cut -f 1,2 %s > %s"
        % (info_file, os.path.join(viral_dir, "viral_accn_taxid.dmp"))
    )

    # perform tests
    gids_1 = Counter(get_accs("viral_database.fasta"))
    gids_2 = Counter([l.split()[0] for l in open("viral_accn_taxid.dmp")])
    assert set(gids_1) == set(gids_2), (
        "taxonomy/viral_seqs_info not matching with fasta"
    )
    duplicates = [k for k, v in gids_1.items() if v > 1]
    if duplicates:
        warnings.warn(
            "Duplicate sequences in viral_database.fasta: %s"
            % " ".join(duplicates)
        )
        logging.warning(
            "Duplicate sequences in viral_database.fasta: %s",
            " ".join(duplicates),
        )
    for l in open("viral_database.fasta"):
        if ">" in l and not l.startswith(">") or l.count(">") > 1:
            warnings.warn("Invalid line in viral_database.fasta: %s" % l)
            logging.warning("Invalid line in viral_database.fasta: %s" % l)


def main(args):
    DB_DIR = os.path.expandvars(args.dbdir)
    logging.info("now in update_db")
    logging.info("Database real path: %s" % os.path.realpath(DB_DIR))
    if bool(args.viral) + args.bact + args.fungal > 1:
        logging.error(
            "update either viral or bacterial or fungal in a single call"
        )
        sys.exit("update either viral or bacterial or fungal in a single call")
    if args.viral:
        virupdate(DB_DIR, args.viral, args.picked, args.update_min_date)
