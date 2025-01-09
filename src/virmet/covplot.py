#!/usr/bin/env python3

"""Plot coverage for a specific organism by realigning all viral reads

:param outdir: directory with results of ``wolfpack``
:param organism: string identifying the desired organism, the program will identify the best
matching sequence among those starting with ``organism``
"""
import os
import logging
import subprocess
import sys
import shlex
from warnings import warn
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser
from virmet.common import run_child, DB_DIR

covpl_exe = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "scripts",
    "covplot.R",
)


def best_species(orgs_file, org_name, best_spec_field_type="ssciname"):
    """Go through tsv file with organism | reads columns (sorted decreasingly by reads) and extract the one with
    most reads by those starting with ``org_name``.
    Example: let's assume that ``org_file`` is as follows::

        Human adenovirus 16     4024
        Human adenovirus 21     1615
        Simian adenovirus 33    1120
        Human adenovirus 7d2    1035
        Human poliovirus 2      663
        Human adenovirus 3+7    360
        Human poliovirus 3      262

    Then, if ``org_name='Human adenovirus'``, ``best_species`` will return ``Human adenovirus 16``.
    If ``org_name='Human adenovirus 2'``, ``best_species`` will return ``Human adenovirus 21` and so on.

    :param orgs_file: path to tab separated file with [organism | reads] columns in decreasing order by reads
    :param org_name: string used to select the organism with the most reads among those starting with ``org_name``

    :returns: string of the best organism
    """
    orgs_list = pd.read_csv(orgs_file, sep="\t", header=0)
    logging.info("Reading an orgs_list file of size %s", str(orgs_list.shape))
    # assert decreasing sorted
    diff = orgs_list["reads"] - orgs_list["reads"].shift(1)
    assert (diff > 0).sum() == 0, diff
    # criterion is "startswith"
    # criterion = orgs_list['sscinames'].map(lambda x: x.startswith(organism))
    criterion = (
        orgs_list.loc[:, "ssciname"].str.startswith(org_name).fillna(False)
    )
    matching_orgs = orgs_list[criterion]
    logging.info("Found %d matchings", orgs_list.shape[0])

    # organism matching that given on command line with most reads is the first
    # W.O. this assumes descending order of reads
    # return str(matching_orgs.iloc[0].ssciname)
    if best_spec_field_type == "ssciname":
        return str(matching_orgs.iloc[0].ssciname)
    else:
        return str(matching_orgs.iloc[0].stitle)


def main(args):
    """Extract the best species, realign reads, run ``covplot.R`` script to create the plot"""
    import datetime

    outdir = args.outdir
    organism = args.organism

    assert os.path.isdir(outdir), "Where is the output dir? Check the path."

    org_file = os.path.join(outdir, "orgs_list.tsv")
    best_spec = best_species(org_file, organism, "ssciname")

    # parse blast results
    blast_file = os.path.join(outdir, "unique.tsv.gz")
    unique = pd.read_csv(blast_file, sep="\t", header=0, compression="gzip")
    matching_reads = unique[unique["ssciname"] == best_spec]
    if matching_reads.empty:
        best_spec = best_species(org_file, organism, "sstitle")
        matching_reads = unique[unique["stitle"] == best_spec]

    best_seqids = (
        matching_reads.groupby("sseqid").size().sort_values(ascending=False)
    )
    try:
        dsc, acc = str(best_seqids.index.tolist()[0]).split("|")[:2]
    except ValueError:
        dsc = "None"
        acc = str(best_seqids.index.tolist()[0])
    logging.info("Best hit in blast results: %s accession:%s", dsc, acc)

    # copy single genome, index, align viral_reads
    os.chdir(outdir)
    organism = organism.replace(" ", "_").replace("/", "_")
    try:
        os.mkdir(organism)
    except FileExistsError:
        warn(
            "directory %s exists already: delete it to run covplot from scratch"
            % organism
        )
    os.chdir(organism)

    if os.path.exists("single.fasta"):
        warn("Reusing single.fasta")
        best_seq = SeqIO.parse("single.fasta", "fasta")
    else:
        viral_db = os.path.join(DB_DIR, "viral_nuccore/viral_database.fasta")
        time1 = datetime.datetime.now()
        # best_seq = [s for s in SeqIO.parse(viral_db, 'fasta') if acc in s.id]
        with open(viral_db) as handle:
            best_seq = [
                s for name, s in SimpleFastaParser(handle) if acc in name
            ]
        print("best_seq found", file=sys.stderr)
        print(datetime.datetime.now() - time1, file=sys.stderr)
        sr = [SeqRecord(Seq(best_seq[0]), id=acc, description="")]
        SeqIO.write(sr, "single.fasta", "fasta")

    seq_len = len(list(best_seq)[0])

    bam_file = "single_sorted.bam"
    if os.path.exists(bam_file):
        warn("Reusing alignment")
        logging.info("Refusing to rerun alignment")
    else:
        run_child("bwa index single.fasta")
        logging.info("Aligning viral reads")
        run_child(
            "bwa mem -t 8 single.fasta ../viral_reads.fastq.gz 2> /dev/null | samtools view -@ 2 -u - | samtools sort -@ 2 -O bam -T tmp -o %s -"
            % bam_file
        )
        run_child("samtools index %s" % bam_file)
    n_reads = int(
        subprocess.check_output(
            'samtools stats %s | grep ^SN | grep "reads mapped:" | cut -f 3'
            % bam_file,
            shell=True,
        ).strip()
    )
    depth_file = "depth.txt"
    if os.path.exists(depth_file):
        warn("Reusing depth file")
    else:
        run_child(
            "samtools depth -a -q 0 -Q 0 %s > %s" % (bam_file, depth_file)
        )
    image_name = organism + "_coverage.pdf"
    logging.info("Plotting coverage")
    perc_obs = subprocess.check_output(
        shlex.split(
            "Rscript %s %s %s %s %s %d"
            % (covpl_exe, depth_file, acc, seq_len, image_name, n_reads)
        )
    )
    try:
        perc_obs_string = perc_obs.decode("ascii").split()[1]
    except IndexError:
        perc_obs_string = "NA"
    print(
        "acc:%s seq_len:%s n_reads:%d perc_obs:%s"
        % (acc, seq_len, n_reads, perc_obs_string)
    )

    return best_species


if __name__ == "__main__":
    import sys

    args_main = {"outdir": sys.argv[1], "organism": sys.argv[2]}
    main(args_main)
