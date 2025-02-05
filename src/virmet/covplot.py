#!/usr/bin/env python3

"""Plot coverage for a specific organism by realigning all viral reads

:param outdir: directory with results of ``wolfpack``
:param organism: string identifying the desired organism, the program will identify the best
matching sequence among those starting with ``organism``
"""

import datetime
import logging
import os
import subprocess
import sys
from warnings import warn
from importlib import resources

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqRecord import SeqRecord

from virmet.common import DB_DIR, run_child


def best_species(orgs_list, org_name, best_spec_field_type="ssciname"):
    """Go through tsv file with organism | reads columns (sorted decreasingly by reads) and extract the one with
    most reads by those starting with ``org_name``.
    Example: let's assume that ``orgs_list`` is as follows::

        Human adenovirus 16     4024
        Human adenovirus 21     1615
        Simian adenovirus 33    1120
        Human adenovirus 7d2    1035
        Human poliovirus 2      663
        Human adenovirus 3+7    360
        Human poliovirus 3      262

    Then, if ``org_name='Human adenovirus'``, ``best_species`` will return ``Human adenovirus 16``.
    If ``org_name='Human adenovirus 2'``, ``best_species`` will return ``Human adenovirus 21` and so on.

    :param orgs_list: already opened tab separated file with [organism | reads] columns in decreasing order by reads
    :param org_name: string used to select the organism with the most reads among those starting with ``org_name``

    :returns: string of the best organism
    """
    logging.info("Reading an orgs_list file of size %s" % str(orgs_list.shape))
    # assert decreasing sorted
    diff = orgs_list["reads"] - orgs_list["reads"].shift(1)
    assert (diff > 0).sum() == 0, diff
    criterion = (
        orgs_list.loc[:, "ssciname"].str.startswith(org_name).fillna(False)
    )
    matching_orgs = orgs_list[criterion]
    logging.info("Found %d matchings" % orgs_list.shape[0])

    # organism matching that given on command line with most reads is the first
    # W.O. this assumes descending order of reads
    if best_spec_field_type == "ssciname":
        return str(matching_orgs.iloc[0].ssciname)
    else:
        return str(matching_orgs.iloc[0].stitle)


def infer_species(
    orgs_file, read_len=151, reads_cutoff=3, covered_score_cutoff=10
):
    """Go through orgs_file.tsv and infer the viral species
    that are worth considering.

        Criteria of exclusion:
        - sum(reads) of a specific specie < reads_cutoff
        - weighted.mean(covered_score, reads) < covered_score_cutoff
        - org_ssciname contains one of the blocked / excluded patterns (phage, bovine, etc.)


    :param orgs_file: path to tab separated file with [organism | reads | seq_len | etc] columns.
    :param read_len: lenght of the sequencing reads that are expeced.
    :param reads_cutoff: threshold (minimum number of mapped reads) to filter out the organism.
    :param covered_score_cutoff: threshold (minimum covered score) to filter out the organism.

    :returns: dataframe of organisms that require covplots
    """
    orgs_list = pd.read_csv(orgs_file, sep="\t", header=0)

    # Calculate number of reads for each organism
    orgs_readsum = orgs_list.groupby("ssciname")["reads"].sum()
    orgs_readsum.drop(
        orgs_readsum[orgs_readsum < reads_cutoff].index, inplace=True
    )
    # Filter out viruses not fulfilling the first criteria
    orgs_list.drop(
        orgs_list[~orgs_list["ssciname"].isin(orgs_readsum.index)].index,
        inplace=True,
    )

    # Calculate the covered score
    orgs_list["covered_percent"] = (
        100 * orgs_list["covered_region"] / orgs_list["seq_len"]
    )
    orgs_list["covered_percent_expected"] = 100 * (
        1 - np.exp(-orgs_list["reads"] * read_len / orgs_list["seq_len"])
    )
    orgs_list["covered_score"] = (
        100
        * orgs_list["covered_percent"]
        / orgs_list["covered_percent_expected"]
    )
    # Calculate the weighted mean for each species
    orgs_coverage = orgs_list.groupby("ssciname").apply(
        lambda x: np.average(x["covered_score"], weights=x["reads"]),
        include_groups=False,
    )
    orgs_coverage.drop(
        orgs_coverage[orgs_coverage < covered_score_cutoff].index, inplace=True
    )
    # Filter out viruses not fulfilling the second griteria
    orgs_list.drop(
        orgs_list[~orgs_list["ssciname"].isin(orgs_coverage.index)].index,
        inplace=True,
    )

    # Define patterns that should be excluded or blocked
    phages_patterns = r"emesvirus zinderi|tunavirus|phage|escherichia|streptococcus|staphylococcus|bacillus|actinomyces|ostreococcus|myoviridae|clostridium|shigella|haemophilus"
    uninteresting_patterns = r"endogenous retrovirus|baboon|bovine|ungulate|bosavirus|betabaculovirus|porcellio scaber|sheep"
    # Drop organisms with blocked or excluded patterns in the name
    orgs_list.drop(
        orgs_list[
            orgs_list["ssciname"].str.lower().str.contains(phages_patterns)
        ].index,
        inplace=True,
    )
    orgs_list.drop(
        orgs_list[
            orgs_list["ssciname"]
            .str.lower()
            .str.contains(uninteresting_patterns)
        ].index,
        inplace=True,
    )

    # Return remaining dataframe
    orgs_list.drop(
        labels=["covered_percent", "covered_percent_expected", "covered_score"],
        axis=1,
        inplace=True,
    )
    return orgs_list


def run_covplot(outdir, n_proc):
    """Extract the best species, realign reads, run ``covplot.R`` script to create the plot"""

    assert os.path.isdir(outdir), "Ensure that output directory exists"

    extended_org_file = os.path.join(outdir, "orgs_list.tsv")
    orgs_list = infer_species(extended_org_file)

    for organism in set(orgs_list["ssciname"]):
        best_spec = best_species(orgs_list, organism, "ssciname")

        # parse blast results
        blast_file = os.path.join(outdir, "unique.tsv.gz")
        unique = pd.read_csv(blast_file, sep="\t", header=0, compression="gzip")
        matching_reads = unique[unique["ssciname"] == best_spec]
        if matching_reads.empty:
            best_spec = best_species(orgs_list, organism, "sstitle")
            matching_reads = unique[unique["stitle"] == best_spec]

        best_seqids = (
            matching_reads.groupby("sseqid").size().sort_values(ascending=False)
        )
        try:
            dsc, acc = str(best_seqids.index.tolist()[0]).split("|")[:2]
        except ValueError:
            dsc = "None"
            acc = str(best_seqids.index.tolist()[0])
        logging.info("Best hit in blast results: %s accession:%s" % (dsc, acc))

        # copy single genome, index, align viral_reads
        organism = organism.replace(" ", "_").replace("/", "_")
        try:
            os.mkdir(os.path.join(outdir, organism))
        except FileExistsError:
            warn(
                "directory %s exists already: delete it to run covplot from scratch"
                % organism
            )

        single_fasta = os.path.join(outdir, organism, "single.fasta")
        if os.path.exists(single_fasta):
            warn("Reusing single.fasta")
            best_seq = SeqIO.parse(single_fasta, "fasta")
        else:
            viral_db = os.path.join(
                DB_DIR, "viral_nuccore/viral_database.fasta"
            )
            time1 = datetime.datetime.now()
            with open(viral_db) as handle:
                best_seq = [
                    s for name, s in SimpleFastaParser(handle) if acc in name
                ]
            print("best_seq found", file=sys.stderr)
            print(datetime.datetime.now() - time1, file=sys.stderr)
            sr = [SeqRecord(Seq(best_seq[0]), id=acc, description="")]
            SeqIO.write(sr, single_fasta, "fasta")

        seq_len = len(list(best_seq)[0])

        bam_file = os.path.join(outdir, organism, "single_sorted.bam")
        if os.path.exists(bam_file):
            warn("Reusing alignment")
            logging.info("Refusing to rerun alignment")
        else:
            run_child("bwa-mem2 index %s" % single_fasta)
            logging.info("Aligning viral reads")
            run_child(
                "bwa-mem2 mem -t %d %s %s/viral_reads.fastq.gz 2> /dev/null | samtools view -@ %d -u - | samtools sort -@ %d -O bam -T tmp -o %s -"
                % (n_proc, single_fasta, outdir, n_proc, n_proc, bam_file)
            )
            run_child("samtools index %s -@ %d" % (bam_file, n_proc))
        n_reads = int(
            subprocess.check_output(
                'samtools stats %s -@ %d | grep ^SN | grep "reads mapped:" | cut -f 3'
                % (bam_file, n_proc),
                shell=True,
            ).strip()
        )
        depth_file = os.path.join(outdir, organism, "depth.txt")
        if os.path.exists(depth_file):
            warn("Reusing depth file")
        else:
            run_child(
                "samtools depth -a -q 0 -Q 0 %s > %s" % (bam_file, depth_file)
            )
        image_name = os.path.join(outdir, organism, organism + "_coverage.pdf")
        logging.info("Plotting coverage")
        perc_obs = subprocess.check_output(
            [
                "Rscript",
                resources.files("virmet").joinpath("covplot.R"),
                str(depth_file),
                str(acc),
                str(seq_len),
                str(image_name),
                str(n_reads),
            ]
        )
        try:
            perc_obs_string = perc_obs.decode("ascii").split()[1]
        except IndexError:
            perc_obs_string = "NA"
        print(
            "acc:%s seq_len:%s n_reads:%d perc_obs:%s"
            % (acc, seq_len, n_reads, perc_obs_string)
        )
