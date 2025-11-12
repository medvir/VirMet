#!/usr/bin/env python3

"""Runs on all samples of a MiSeq run or on a single fastq file"""

import glob
import gzip
import logging
import multiprocessing as mp
import os
import re
import shlex
import shutil
import subprocess
import tempfile
import warnings

import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from virmet.__init__ import __version__
from virmet.common import run_child, n_proc
from virmet.covplot import run_covplot
from virmet.tidytable import run_tidytable

blast_cov_threshold = 75.0
blast_ident_threshold = 75.0


def strip(str_):
    """Make the strip method a function"""
    return str_.strip()


def span_coverage(row):
    """Return a set of covered positions given a blast alignment entry."""
    start, end = min(row.sstart, row.send), max(row.sstart, row.send)
    return set(range(start, end + 1))


def merge_coverage(series):
    """Join the covered ranges and return the length."""
    covered = set()
    for s in series.tolist():
        covered |= s
    return len(covered)


def run_blast(fasta_file, Database, n_proc, out_file, delete_fasta=False):
    """Run blastn with task megablast"""
    cml = (
        "blastn -task megablast \
            -query %s -db %s \
            -num_threads %d \
            -out %s \
            -outfmt '6 qseqid sseqid ssciname stitle pident qcovs score length mismatch gapopen qstart qend sstart send staxid'"
        % (fasta_file, Database, n_proc, out_file)
    )
    run_child(cml)
    if delete_fasta:
        os.remove(fasta_file)


def get_nodes_names(dir_name):
    """Look for dmp taxonomy files in dir_name, read and return them

    Some lines from https://github.com/zyxue/ncbitax2lin are used.
    """
    nodes_file = os.path.join(dir_name, "nodes.dmp.gz")
    logging.info("reading nodes file %s" % nodes_file)
    colnames = [
        "tax_id",
        "parent_tax_id",
        "rank",
        "embl_code",
        "division_id",
        "inherited_div_flag",
        "genetic_code_id",
        "inherited_GC__flag",
        "mitochondrial_genetic_code_id",
        "inherited_MGC_flag",
        "GenBank_hidden_flag",
        "hidden_subtree_root_flag",
        "comments",
    ]
    nodes = pd.read_csv(
        nodes_file,
        sep="|",
        header=None,
        index_col=False,
        names=colnames,
        compression="gzip",
    )
    # To get rid of flanking tab characters
    nodes["rank"] = nodes["rank"].apply(strip)
    nodes["embl_code"] = nodes["embl_code"].apply(strip)
    nodes["comments"] = nodes["comments"].apply(strip)
    nodes.set_index("tax_id", inplace=True)

    names_file = os.path.join(dir_name, "names.dmp.gz")
    logging.info("reading names file %s" % names_file)
    colnames = ["tax_id", "taxon_name", "alt_name", "class_name"]
    names = pd.read_csv(
        names_file,
        sep="|",
        header=None,
        index_col=False,
        names=colnames,
        compression="gzip",
    )

    names["taxon_name"] = names["taxon_name"].apply(strip)
    names["alt_name"] = names["alt_name"].apply(strip)
    names["class_name"] = names["class_name"].apply(strip)
    sci_names = names[names["class_name"] == "scientific name"]
    sci_names.set_index("tax_id", inplace=True)

    return nodes, sci_names


def get_parent_species(inrow, nodes, names):
    """Travel the whole taxonomy above query until species, return the species name

    :param inrow:
    :param nodes:
    :param names:

    """
    query = inrow["tax_id"]
    if query == 0:
        return "NA"
    i = 100  # failsafe method to avoid infinite loops (shame on me)
    while i:
        try:
            node_row = nodes.loc[query]
        except KeyError:  # databases are sometimes not aligned
            warnings.warn("Could not find parent for node %s" % query)
            return "NA"
        rank = node_row[
            "rank"
        ]  # rank also a method of DataFrame, can't use node_row.rank
        name_row = names.loc[query]
        org_name = name_row.taxon_name
        if rank == "species":
            break
        query = node_row.parent_tax_id
        i -= 1
    return org_name


def hunter(fq_file, out_dir, n_proc):
    """runs quality filter on a fastq file with fastp,
    returns output directory
    """

    logging.debug("hunter will run on %s processors" % n_proc)
    if "L001" in fq_file:
        s_dir = "_".join(os.path.split(fq_file)[1].split("_")[:2])
        s_dir = os.path.join(out_dir, s_dir)
        try:
            os.mkdir(s_dir)
        except FileExistsError:
            logging.debug("entering %s already existing" % s_dir)
    else:
        s_dir = out_dir

    # skip if this is a hot run
    if os.path.exists(os.path.join(s_dir, "fastp.html")):
        logging.info("hunter was already run in %s, skipping" % s_dir)
        return s_dir

    # trim, discard short reads and filter with fastp
    good_file = os.path.join(s_dir, "good.fastq")
    cml = (
        "fastp -i %s -o %s -A \
            --cut_front --cut_tail --cut_mean_quality 20 --length_required 75 \
            --average_qual 20 --low_complexity_filter --complexity_threshold 30 \
            --thread  %d --json %s/fastp.json --html %s/fastp.html > %s/fastp.err 2>&1"
        % (fq_file, good_file, n_proc, s_dir, s_dir, s_dir)
    )
    run_child(cml)

    # Parsing statistics to stats.tsv file
    stats_file = os.path.join(s_dir, "stats.tsv")
    oh = open(stats_file, "w+")
    with open("%s/fastp.json" % s_dir) as f:
        useful = [next(f) for _ in range(35)]
        raw_reads = int(
            re.search(r'"total_reads": *(\d+)', str(useful)).group(1)
        )
        too_short = int(
            re.search(r'"too_short_reads": *(\d+)', str(useful)).group(1)
        )
        low_entropy = int(
            re.search(r'"low_complexity_reads": *(\d+)', str(useful)).group(1)
        )
        low_quality = int(
            re.search(r'"low_quality_reads": *(\d+)', str(useful)).group(1)
        )
        passed_filter = int(
            re.search(r'"passed_filter_reads": *(\d+)', str(useful)).group(1)
        )
    oh.write("raw_reads\t%d\n" % raw_reads)
    oh.write("trimmed_too_short\t%d\n" % too_short)
    oh.write("low_entropy\t%d\n" % low_entropy)
    oh.write("low_quality\t%d\n" % low_quality)
    oh.write("passing_filter\t%d\n" % passed_filter)
    oh.close()

    sample_info = os.path.join(s_dir, "sample_info.txt")
    with open(sample_info, "a") as oh:
        oh.write(f"VirMet version: {__version__}\n")

    return s_dir


def victor(input_reads, contaminant, n_proc):
    """decontaminate reads by aligning against contaminants with bwa and removing
    reads with alignments
    """

    rf_head = input_reads.split(".")[0]
    cont_name = os.path.split(contaminant)[1]
    sam_name = "%s_%s.sam" % (rf_head, cont_name)
    err_name = "%s_%s.err" % (rf_head, cont_name)
    clean_name = os.path.splitext(sam_name)[0] + ".fastq"

    # skipping if hot run
    if os.path.exists(err_name):
        logging.info("decontamination already performed, skipping")
        return os.path.split(clean_name)[1]

    # alignment with bwa
    cont_real_link = os.path.realpath(contaminant)
    logging.info("Database real path: %s" % cont_real_link)
    cml = (
        "bwa mem -t %d -R '@RG\\tID:foo\\tSM:bar\\tLB:library1' -T 75 -M %s %s 2> \
    %s | samtools view -h -F 4 - > %s"
        % (n_proc, contaminant, input_reads, err_name, sam_name)
    )
    logging.debug(
        "running bwa %s %s on %d cores" % (cont_name, rf_head, n_proc)
    )
    run_child(cml)

    # reading sam file to remove reads with hits
    # test if an object is in set is way faster than in list
    mapped_reads = set(
        run_child('grep -v "^@" %s | cut -f 1' % sam_name).strip().split("\n")
    )
    try:  # if no matches, empty string is present
        mapped_reads.remove("")
    except KeyError:
        pass

    stats_file = os.path.join(os.path.split(input_reads)[0], "stats.tsv")
    oh = open(stats_file, "a")
    oh.write("matching_%s\t%d\n" % (cont_name, len(mapped_reads)))
    oh.close()

    output_handle = open(clean_name, "w")
    logging.debug(
        "Cleaning reads in %s with alignments in %s" % (input_reads, sam_name)
    )
    logging.debug("Writing to %s" % clean_name)
    if input_reads.endswith(".gz"):
        cont_handle = gzip.open(input_reads)
    else:
        cont_handle = open(input_reads)
    c = 0
    # Using FastqGeneralIterator allows fast performance
    for title, seq, qual in FastqGeneralIterator(cont_handle):
        if title.split()[0] not in mapped_reads:
            c += 1
            output_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
            if c % 100000 == 0:
                logging.debug("written %d clean reads" % c)
    logging.info("written %d clean reads" % c)
    output_handle.close()
    cont_handle.close()

    if os.path.split(input_reads)[1] != "good.fastq":
        os.remove(input_reads)

    return os.path.split(clean_name)[1]


def victor_bact_fungal(input_reads, decont_db, n_proc):
    """decontaminate reads by aligning against contaminants with kraken2"""

    rf_head = input_reads.split(".")[0]
    err_name = "%s_bact_fungi.err" % rf_head
    clean_name = "%s_bact_fungi.fastq" % rf_head
    report_name = "%s_bact_fungi.txt" % rf_head

    # skipping if hot run
    if os.path.exists(err_name):
        logging.info("decontamination already performed, skipping")
        return os.path.split(clean_name)[1]

    # alignment with kraken2
    cont_real_link = os.path.realpath(decont_db)
    logging.info("Database real path: %s" % cont_real_link)
    cml = (
        "kraken2 --db %s --threads %d --minimum-hit-groups 4 --confidence 0.5 --report %s --unclassified-out %s %s > %s 2>&1"
        % (decont_db, n_proc, report_name, clean_name, input_reads, err_name)
    )
    logging.debug("running kraken2 with %s on %d cores" % (rf_head, n_proc))
    run_child(cml)

    # Calculate statistics
    stats_kraken = pd.read_csv(
        report_name,
        sep="\t",
        names=["Percentage", "Total", "Assigned", "Code", "Taxid", "Sciname"],
    )
    try:
        mapped_bact = int(
            stats_kraken.loc[stats_kraken["Taxid"] == 2, "Total"].iloc[0]
        )
    except (IndexError, KeyError, ValueError):
        mapped_bact = 0
    try:
        mapped_fungi = int(
            stats_kraken.loc[stats_kraken["Taxid"] == 4751, "Total"].iloc[0]
        )
    except (IndexError, KeyError, ValueError):
        mapped_fungi = 0
    try:
        mapped_cellorg = (
            int(
                stats_kraken.loc[stats_kraken["Taxid"] == 131567, "Total"].iloc[
                    0
                ]
            )
            - mapped_bact
            - mapped_fungi
        )
    except (IndexError, KeyError, ValueError):
        mapped_cellorg = 0

    # Write statistics to output
    stats_file = os.path.join(os.path.split(input_reads)[0], "stats.tsv")
    oh = open(stats_file, "a")
    oh.write("matching_bacteria\t%d\n" % mapped_bact)
    oh.write("matching_fungi\t%d\n" % mapped_fungi)
    oh.write("matching_other_cellular_organisms\t%d\n" % mapped_cellorg)
    oh.close()

    if os.path.split(input_reads)[1] != "good.fastq":
        os.remove(input_reads)

    return os.path.split(clean_name)[1]


def viral_blast(file_in, n_proc, nodes, names, out_dir, DB_DIR):
    """runs blast against viral database"""

    viral_reads = os.path.join(out_dir, "viral_reads.fastq.gz")
    undet_reads = os.path.join(out_dir, "undetermined_reads.fastq.gz")
    # on hot start, blast again all decontaminated reads
    if os.path.exists(viral_reads) and os.path.exists(undet_reads):
        run_child("gunzip -c %s %s > %s" % (viral_reads, undet_reads, file_in))
        os.remove(viral_reads)
        os.remove(undet_reads)

    # streams will be used during the execution
    child_dir = os.path.split(file_in)[0]

    stats_file = os.path.join(child_dir, "stats.tsv")
    unique_file = os.path.join(child_dir, "unique.tsv")
    oh = open(stats_file, "a")

    decont_reads = os.path.join(child_dir, "hq_decont_reads.fastq")
    if not os.path.exists(decont_reads):
        os.rename(file_in, decont_reads)
    fasta_file = os.path.join(child_dir, "hq_decont_reads.fasta")
    run_child("seqtk seq -A %s > %s" % (decont_reads, fasta_file))
    try:
        tot_seqs = int(run_child('grep -c "^>" %s' % fasta_file).strip())
    except AttributeError:  # deals with empty file
        tot_seqs = 0
        logging.info("No reads left after decontamination")

    oh.write("reads_to_blast\t%d\n" % tot_seqs)

    if tot_seqs == 0:
        oh.write("viral_reads\t0\n")
        oh.write("undetermined_reads\t0\n")
        oh.close()
        return

    # blast needs access to taxdb files to retrieve organism name
    os.environ["BLASTDB"] = DB_DIR
    logging.info("runnning on %d cores" % n_proc)
    DB_real_path = os.path.realpath(
        os.path.join(DB_DIR, "viral_nuccore/viral_db")
    )
    logging.info("Database real path: %s" % DB_real_path)

    if (
        tot_seqs > 50000
    ):  # If there are less than 50000 sequences, it's not worth it
        max_n = (
            int((tot_seqs / n_proc / 2) + 1) * 2
        )  # Multiply by 2 because it's a fasta
        # We want to split in n_proc processors
        cml = "split -l %d %s %s/splitblastn --additional-suffix=.fasta" % (
            max_n,
            fasta_file,
            child_dir,
        )
        logging.info("Splitting fasta into smaller files for faster blastn")
        run_child(cml)
        # Find all splitted files and run blast on them
        splitted_files = glob.glob(os.path.join(child_dir, "splitblastn*"))
        myblast = []
        n_out = 1
        for split_file in splitted_files:
            tmp_output = os.path.join(child_dir, "tmp_blastn_%d.tsv" % n_out)
            myblast.append((split_file, DB_real_path, 1, tmp_output, True))
            n_out += 1
        # Run blastn commands in parallel for faster performance
        with mp.Pool(n_proc) as pool_blast:
            _ = pool_blast.starmap_async(run_blast, myblast, chunksize=2)
            # Wait to finish all the blastn before going to the next steps
            pool_blast.close()
            pool_blast.join()
        tmp_outputs = glob.glob(os.path.join(child_dir, "tmp_blastn_*"))
        with open(unique_file, "w") as outfile:
            outfile.write(
                "qseqid\tsseqid\tssciname\tstitle\tpident\tqcovs\tscore\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tstaxid\n"
            )
            for fname in tmp_outputs:
                with open(fname) as infile:
                    outfile.write(infile.read())
                os.remove(fname)
    else:
        run_blast(fasta_file, DB_real_path, n_proc, unique_file)
        with open(unique_file, "r") as original:
            data = original.read()
        with open(unique_file, "w") as modified:
            modified.write(
                "qseqid\tsseqid\tssciname\tstitle\tpident\tqcovs\tscore\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tstaxid\n"
                + data
            )

    # Take only the best hit, aka hsps
    run_child(
        "awk '!a[$1]++' %s >> %s/mtmp; mv %s/mtmp %s"
        % (unique_file, child_dir, child_dir, unique_file)
    )
    logging.debug("saving blast database info")
    cml = shlex.split("blastdbcmd -db %s -info" % DB_real_path)
    with open("%s/blast_info.txt" % child_dir, "a") as boh:
        subprocess.call(cml, stdout=boh)

    logging.debug("filtering and grouping by hit sequence")

    with open(unique_file, "r") as f:
        num_lines = sum(1 for _ in f)

    if num_lines == 1:
        warnings.warn("No hits, unique.tsv file is empty")
        return

    hits = pd.read_csv(
        unique_file, index_col="qseqid", delimiter="\t", na_filter=False
    )
    logging.debug("found %d hits" % hits.shape[0])
    # select according to identity and coverage, count occurrences
    good_hits = hits[
        (hits.pident > blast_ident_threshold)
        & (hits.qcovs > blast_cov_threshold)
    ]
    matched_reads = good_hits.shape[0]
    logging.debug(
        "%d hits passing coverage and identity filter" % matched_reads
    )
    oh.write("viral_reads\t%s\n" % matched_reads)
    unknown_reads = tot_seqs - matched_reads
    oh.write("undetermined_reads\t%d\n" % unknown_reads)
    oh.close()

    if matched_reads == 0:  # deals with no good_hits
        warnings.warn("No hits")
        return

    # create a column for accession number
    good_hits = good_hits.assign(
        accn=good_hits["sseqid"].apply(
            lambda x: re.search(r"([A-Z]+_?\d*\.?\d*)", x).group(1)
        )
    )
    good_hits = good_hits.rename(columns={"staxid": "tax_id"})

    viral_info_file = os.path.join(DB_DIR, "viral_nuccore/viral_seqs_info.tsv")
    viral_info = pd.read_table(
        viral_info_file, names=["accn", "TaxId", "seq_len", "Organism"]
    )
    good_hits = pd.merge(good_hits, viral_info, on="accn")
    # if blastn gives no taxid and scientific name, fill these col from viral_seqs_info.tsv file
    good_hits.loc[:, "ssciname"] = (
        good_hits.loc[:, "ssciname"].fillna(good_hits["Organism"]).astype(str)
    )
    good_hits.loc[:, "tax_id"] = (
        good_hits.loc[:, "tax_id"].fillna(good_hits["TaxId"]).astype(int)
    )
    # fill the species and the covered range on subject sequence
    good_hits["species"] = good_hits.apply(
        lambda row: get_parent_species(row, nodes, names), axis=1
    )
    good_hits["covered_region"] = good_hits.apply(
        lambda row: span_coverage(row), axis=1
    )
    if good_hits.isnull().any().any():
        logging.error(
            "There is 'nan' in the result of the blastn after selecting good hits."
        )

    # now summarise and write the covered region length
    ds = good_hits.groupby(
        ["accn", "stitle", "ssciname", "species", "tax_id"]
    ).agg({"covered_region": merge_coverage})
    ds["reads"] = good_hits.groupby(
        ["accn", "stitle", "ssciname", "species", "tax_id"]
    ).size()
    ds = ds.reset_index()

    viral_info = viral_info.drop(columns=["TaxId", "Organism"])
    ds = pd.merge(ds, viral_info)
    ds = ds.loc[
        :,
        [
            "species",
            "accn",
            "reads",
            "stitle",
            "ssciname",
            "covered_region",
            "seq_len",
        ],
    ]
    ds = ds.sort_values(
        by=["reads", "covered_region"], ascending=[False, False]
    )
    orgs_list = os.path.join(child_dir, "orgs_list.tsv")
    ds.to_csv(orgs_list, header=True, sep="\t", index=False)


def cleaning_up(cleaned_dir):
    """sift reads into viral/unknown, compresses and removes files"""

    # selects reads with coverage and identity higher than 75
    unique_file = os.path.join(cleaned_dir, "unique.tsv")
    all_reads = os.path.join(cleaned_dir, "hq_decont_reads.fastq")
    undet_reads = os.path.join(cleaned_dir, "undetermined_reads.fastq")

    with open(unique_file, "r") as f:
        num_lines = sum(1 for _ in f)

    if num_lines != 1:
        df = pd.read_csv(unique_file, sep="\t")
        viral_ids = set(
            df[
                (df.qcovs > blast_cov_threshold)
                & (df.pident > blast_ident_threshold)
            ].qseqid
        )
        viral_c = 0
        undet_c = 0
        all_handle = open(all_reads)
        viral_reads = os.path.join(cleaned_dir, "viral_reads.fastq")
        undet_handle = open(undet_reads, "w")
        viral_handle = open(viral_reads, "w")
        # Using FastqGeneralIterator allows fast performance
        for title, seq, qual in FastqGeneralIterator(all_handle):
            if title.split()[0] not in viral_ids:
                undet_c += 1
                undet_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                if undet_c % 100000 == 0:
                    logging.debug("written %d undet reads" % undet_c)
            else:
                viral_c += 1
                viral_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                if viral_c % 10000 == 0:
                    logging.debug("written %d viral reads" % viral_c)
        undet_handle.close()
        viral_handle.close()
        logging.info("written %d undet reads" % undet_c)
        logging.info("written %d viral reads" % viral_c)

        run_child("gzip -f %s" % viral_reads)
        run_child("gzip -f %s" % undet_reads)
        os.remove(all_reads)
    else:
        os.rename(all_reads, undet_reads)

    # removing and zipping
    for samfile in glob.glob(os.path.join(cleaned_dir, "*.sam")):
        os.remove(samfile)
    for rf in [
        "%s/good.fastq" % cleaned_dir,
        "%s/good_humanGRCh38_bt_ref_bact_fungi.txt" % cleaned_dir,
        "%s/good_humanGRCh38_bt_ref_bact_fungi.fastq" % cleaned_dir,
        "%s/hq_decont_reads.fasta" % cleaned_dir,
        "%s/fastp.json" % cleaned_dir,
    ]:
        try:
            os.remove(rf)
        except FileNotFoundError:
            pass

    run_child("gzip -f %s/unique.tsv" % cleaned_dir)


def main(args):
    """"""

    DB_DIR = os.path.expandvars(args.dbdir)
    # Specify all databases for decontamination steps
    contaminant_db = [
        os.path.join(DB_DIR, "human/bwa/humanGRCh38"),
        os.path.join(DB_DIR, "bovine/bwa/bt_ref"),
    ]

    if args.run:
        miseq_dir = args.run.rstrip("/")
        run_name = os.path.split(miseq_dir)[1]
        if (
            run_name.startswith(("1", "2"))
            and len(run_name.split("-")[-1]) == 5
            and run_name.split("_")[1].startswith("M")
        ):
            try:
                _, machine_name = run_name.split("_")[:2]
                logging.info(
                    "running on run %s from machine %s"
                    % (run_name, machine_name)
                )
            except ValueError:
                logging.info("running on directory %s" % miseq_dir)

            # Check where the base calls actually are
            possible_bc_dir = os.path.join(
                miseq_dir, "Data", "Intensities", "BaseCalls"
            )
            if os.path.isdir(possible_bc_dir):
                bc_dir = possible_bc_dir
            else:
                bc_dir = miseq_dir
        else:
            bc_dir = miseq_dir

        rel_fastq_files = glob.glob("%s/*_S*.fastq*" % bc_dir)
        samples_to_run = [
            os.path.split(fq)[1].split("_")[1] for fq in rel_fastq_files
        ]
        logging.info("samples to run: %s", " ".join(samples_to_run))
        if not args.noctrls:
            all_fastq_files = [os.path.abspath(f) for f in rel_fastq_files]
        else:
            all_fastq_files = [
                os.path.abspath(f)
                for f in rel_fastq_files
                if not re.search(r"ntc-|Undetermin", f)
            ]
    elif args.file:
        all_fastq_files = [os.path.abspath(args.file)]
        logging.info("running on a single file %s" % all_fastq_files[0])
        run_name = os.path.split(args.file)[1].split(".")[0]

    out_dir_final = "virmet_output_%s" % run_name
    out_dir_final = os.path.abspath(out_dir_final)

    if os.path.isdir(out_dir_final):
        logging.info("directory %s already exists" % out_dir_final)
        logging.info("performing hot run, results will be overwritten")
        out_dir = out_dir_final
    else:
        out_dir = os.path.join(
            tempfile.gettempdir(), "virmet_output_%s" % run_name
        )
        out_dir = os.path.abspath(out_dir)
        try:
            os.mkdir(out_dir)
        except OSError:
            shutil.rmtree(out_dir)
            os.mkdir(out_dir)

    # run hunter on all fastq files
    s_dirs = []
    for fq in all_fastq_files:
        logging.info("running hunter on %s" % fq)
        sd = hunter(fq, out_dir, n_proc)
        s_dirs.append(sd)

    # run mapping against contaminants to remove
    cont_reads = "good.fastq"  # first run on good.fastq
    for cont in contaminant_db:
        logging.info("decontamination against %s" % cont)
        for sample_dir in s_dirs:
            logging.info("--- now for sample %s" % sample_dir)
            input_vict = os.path.join(sample_dir, cont_reads)
            decont_reads = victor(
                input_reads=input_vict, contaminant=cont, n_proc=n_proc
            )
        cont_reads = decont_reads  # decontaminated reads are input for next round (equal across samples)

    # Decontamination against bacterial and fungal db
    logging.info("decontamination against bacterial and fungal database")
    bact_fung_db = f"{DB_DIR}/bact_fungi"
    for sample_dir in s_dirs:
        logging.info("--- now for sample %s" % sample_dir)
        input_vict = os.path.join(sample_dir, cont_reads)
        file_to_blast = victor_bact_fungal(
            input_reads=input_vict, decont_db=bact_fung_db, n_proc=n_proc
        )

    logging.info("blasting against viral database")
    logging.info("%d cores will be used" % n_proc)

    logging.info("reading taxonomy files")
    nodes, names = get_nodes_names(DB_DIR)

    for sample_dir in s_dirs:
        logging.info("now sample %s" % sample_dir)
        viral_blast(
            os.path.join(sample_dir, file_to_blast),
            n_proc,
            nodes,
            names,
            sample_dir,
            DB_DIR,
        )
        logging.info("sample %s blasted" % sample_dir)

    logging.info("summarising and cleaning up")
    for sample_dir in s_dirs:
        logging.info("now in %s" % sample_dir)
        cleaning_up(sample_dir)

    if not args.nocovplot:
        for sample_dir in s_dirs:
            run_covplot(sample_dir, n_proc, DB_DIR)

    # Organise results into a tidy table
    if len(s_dirs) > 1:
        run_tidytable(out_dir)
    else:
        orgs_file_path = os.path.join(s_dirs[0], "orgs_list.tsv")
        stats_file_path = os.path.join(s_dirs[0], "stats.tsv")
        if os.path.exists(orgs_file_path):
            shutil.copyfile(orgs_file_path, f"{out_dir}/orgs_species_found.tsv")
        shutil.copyfile(stats_file_path, f"{out_dir}/run_reads_summary.tsv")

    # Move results to final directory
    if os.path.abspath(out_dir) != os.path.abspath(out_dir_final):
        shutil.move(out_dir, out_dir_final)
    return out_dir_final
