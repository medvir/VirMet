#!/usr/bin/env python3
"""
Read and print summary of reads category for the whole run
"""

import os
import sys
import glob
import pandas as pd


def main(args):
    """Everything is defined here"""
    outdir = args.outdir
    run = os.path.abspath(outdir).split("/")[-1].split("virmet_output_")[1]
    try:
        os.chdir(outdir)
    except FileNotFoundError:
        sys.exit("Where is the output dir? Check the path.")

    sample_dirs = glob.glob("*_S*")
    all_reads = pd.DataFrame()
    all_orgs = pd.DataFrame()
    for sd in sample_dirs:
        # parse and save stat files
        stat_file = os.path.join(sd, "stats.tsv")
        df = pd.read_csv(
            stat_file, sep="\t", header=None, names=["category", "reads"]
        )
        df["sample"] = sd
        df["run"] = run
        all_reads = all_reads.append(df)
        # parse and save orgs_list files
        orgs_file = os.path.join(sd, "orgs_list.tsv")
        if os.path.isfile(orgs_file):
            df = pd.read_csv(orgs_file, sep="\t", header=0)
            df["sample"] = sd
            df["run"] = run
            all_orgs = all_orgs.append(df)
        else:
            continue

    all_orgs.to_csv("orgs_species_found.tsv", sep="\t", index=False)
    all_reads.to_csv("run_reads_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    args = {"outdir": sys.argv[1]}
    main(args)
