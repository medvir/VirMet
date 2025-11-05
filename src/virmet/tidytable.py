#!/usr/bin/env python3
"""
Read and print summary of reads category for the whole run
"""

import glob
import os
import sys

import pandas as pd


def run_tidytable(outdir):
    """Everything is defined here"""
    run = os.path.abspath(outdir).split("/")[-1].split("virmet_output_")[1]
    if not os.path.isdir(outdir):
        sys.exit("Where is the output dir? Check the path.")

    sample_dirs = glob.glob("%s/*_S*" % outdir)
    all_reads = pd.DataFrame()
    all_orgs = pd.DataFrame()
    for sd in sample_dirs:
        # parse and save stat files
        stat_file = os.path.join(sd, "stats.tsv")
        df = pd.read_csv(
            stat_file, sep="\t", header=None, names=["category", "reads"]
        )
        df["sample"] = os.path.basename(os.path.dirname(sd + "/"))
        df["run"] = run
        all_reads = pd.concat([all_reads, pd.DataFrame(df)])
        # parse and save orgs_list files
        orgs_file = os.path.join(sd, "orgs_list.tsv")
        if os.path.isfile(orgs_file):
            df = pd.read_csv(orgs_file, sep="\t", header=0)
            df["sample"] = os.path.basename(os.path.dirname(sd + "/"))
            df["run"] = run
            all_orgs = pd.concat([all_orgs, pd.DataFrame(df)])

    all_orgs.to_csv("%s/orgs_species_found.tsv" % outdir, sep="\t", index=False)
    all_reads.to_csv("%s/run_reads_summary.tsv" % outdir, sep="\t", index=False)
