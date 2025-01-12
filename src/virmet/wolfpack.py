#!/usr/bin/env python3

"""Runs on all samples of a MiSeq run or on a single fastq file"""

import glob
import logging
import os
import re
import shlex
import subprocess
import warnings

import pandas as pd
from pkg_resources import DistributionNotFound, get_distribution

from virmet.common import DB_DIR, run_child  # , single_process

try:
    __version__ = get_distribution("virmet").version
except DistributionNotFound:
    # package is not installed
    pass

contaminant_db = [
    "/data/virmet_databases/human/bwa/humanGRCh38",
    "/data/virmet_databases/bacteria_concise/bwa/bact1",
    "/data/virmet_databases/bacteria_concise/bwa/bact2",
    "/data/virmet_databases/bacteria_concise/bwa/bact3",
    "/data/virmet_databases/bacteria_concise/bwa/bact4",
    "/data/virmet_databases/bacteria_concise/bwa/bact5",
    "/data/virmet_databases/fungi/bwa/fungi1",
    "/data/virmet_databases/bovine/bwa/bt_ref",
]
ref_map = {
    "humanGRCh38": "/data/virmet_databases/human/fasta/GRCh38.fasta.gz",
    "bact1": "/data/virmet_databases/bacteria_concise/fasta/bact1.fasta.gz",
    "bact2": "/data/virmet_databases/bacteria_concise/fasta/bact2.fasta.gz",
    "bact3": "/data/virmet_databases/bacteria_concise/fasta/bact3.fasta.gz",
    "bact4": "/data/virmet_databases/bacteria_concise/fasta/bact4.fasta.gz",
    "bact5": "/data/virmet_databases/bacteria_concise/fasta/bact5.fasta.gz",
    "fungi1": "/data/virmet_databases/fungi/fasta/fungi1.fasta.gz",
    "bt_ref": "/data/virmet_databases/bovine/fasta/bt_ref_Bos_taurus_UMD_3.1.1.fasta.gz",
}

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


def get_nodes_names(dir_name):
    """Look for dmp taxonomy files in dir_name, read and return them

    Some lines from https://github.com/zyxue/ncbitax2lin are used.
    """
    nodes_file = os.path.join(dir_name, "nodes.dmp.gz")
    logging.info("reading nodes file %s", nodes_file)
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
    logging.info("reading names file %s", names_file)
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


def hunter(fq_file):
    """runs quality filter on a fastq file with seqtk and prinseq,
    simple parallelisation with xargs, returns output directory
    """
    # from virmet.common import prinseq_exe
    prinseq_exe = "prinseq-lite.pl"
    prinseq_exe = "prinseq"

    try:
        n_proc = min(os.cpu_count(), 16)
        if n_proc == 1:
            n_proc = 2
    except NotImplementedError:
        n_proc = 2

    logging.debug("hunter will run on %s processors", n_proc)
    if "L001" in fq_file:
        s_dir = "_".join(os.path.split(fq_file)[1].split("_")[:2])
        try:
            os.mkdir(s_dir)
        except FileExistsError:
            logging.debug("entering %s already existing", s_dir)
        os.chdir(s_dir)
        s_dir = os.getcwd()
    else:
        s_dir = os.getcwd()

    # skip if this is a hot run
    if os.path.exists("prinseq.err") and os.path.exists("prinseq.log"):
        logging.info("hunter was already run in %s, skipping", s_dir)
        os.chdir(os.pardir)
        return s_dir

    # first occurrence of stats.tsv
    oh = open("stats.tsv", "w+")
    # count raw reads
    if fq_file.endswith("gz"):
        out1 = run_child("gunzip -c %s | wc -l" % fq_file)
    else:
        out1 = run_child("wc -l %s" % fq_file)
    out1 = out1.strip().split()[0]
    n_reads = int(int(out1.strip()) / 4)
    oh.write("raw_reads\t%d\n" % n_reads)

    # trim and discard short reads, count
    logging.debug("trimming with seqtk")
    cml = "trimfq %s | seqtk seq -L 75 - > intermediate.fastq" % fq_file
    out1 = run_child("seqtk " + cml)
    out1 = run_child("wc -l intermediate.fastq")
    out1 = out1.strip().split()[0]

    long_reads = int(int(out1.strip()) / 4)
    short = n_reads - long_reads
    oh.write("trimmed_too_short\t%d\n" % short)

    # We want to split in n_proc processors, so each file has at most
    # (n_reads / n_proc) + 1 reads and 4 times as many lines
    # this fails if there are more cpus than reads!
    max_reads_per_file = int(n_reads / n_proc) + 1
    max_l = max_reads_per_file * 4
    # split and rename
    run_child("split -l %d intermediate.fastq splitted" % max_l)
    os.remove("intermediate.fastq")
    splitted = glob.glob("splitted*")
    n_splitted = len(splitted)
    for i, spf in enumerate(sorted(splitted)):
        os.rename(spf, "splitted%03d.fastq" % i)  # W.O. max 1000 files/cpus

    # filter with prinseq, parallelize with xargs
    logging.debug("filtering with prinseq")
    cml = (
        "-f %%03g 0 %d | xargs -P %d -I {} %s \
            -fastq splitted{}.fastq -lc_method entropy -lc_threshold 70 \
            -log prinseq{}.log -min_qual_mean 20 \
            -out_good ./good{} -out_bad ./bad{} > ./prinseq.err 2>&1"
        % (n_splitted - 1, n_splitted, prinseq_exe)
    )
    run_child("/usr/bin/seq " + cml)

    logging.debug("cleaning up")
    if glob.glob("good???.fastq"):
        run_child("cat good???.fastq > good.fastq")
        run_child("rm good???.fastq")

    if glob.glob("bad???.fastq"):
        run_child("cat bad???.fastq > bad.fastq")
        run_child("rm bad???.fastq")

    if glob.glob("prinseq???.log"):
        run_child("cat prinseq???.log > prinseq.log")
        run_child("rm prinseq???.log")

    run_child("rm splitted*fastq")

    # parsing number of reads deleted because of low entropy
    low_ent = 0
    min_qual = 0
    with open("prinseq.log") as f:
        for l in f:
            match_lc = re.search("lc_method\:\s(\d*)$", l)
            match_mq = re.search("min_qual_mean\:\s(\d*)$", l)
            if match_lc:
                low_ent += int(match_lc.group(1))
            elif match_mq:
                min_qual += int(match_mq.group(1))
    oh.write("low_entropy\t%d\n" % low_ent)
    oh.write("low_quality\t%d\n" % min_qual)

    out1 = run_child("wc -l good.fastq")
    out1 = out1.strip().split()[0]
    n_reads = int(int(out1) / 4)
    lost_reads = n_reads + low_ent + min_qual - long_reads
    if lost_reads > 0:
        logging.error("%d reads were lost", lost_reads)
        warnings.warn("%d reads were lost" % lost_reads, RuntimeWarning)
    oh.write("passing_filter\t%d\n" % n_reads)
    oh.close()

    with open("sample_info.txt", "a") as oh:
        oh.write("VirMet version: %s\n" % __version__)

    os.chdir(os.pardir)
    return s_dir


def victor(input_reads, contaminant):
    """decontaminate reads by aligning against contaminants with bwa and removing
    reads with alignments
    """
    import gzip

    from Bio.SeqIO.QualityIO import FastqGeneralIterator

    try:
        n_proc = min(os.cpu_count(), 16)
    except NotImplementedError:
        n_proc = 2

    rf_head = input_reads.split(".")[0]
    cont_name = os.path.split(contaminant)[1]
    sam_name = "%s_%s.sam" % (rf_head, cont_name)
    err_name = "%s_%s.err" % (rf_head, cont_name)
    clean_name = os.path.splitext(sam_name)[0] + ".fastq"

    # skipping if hot run
    if os.path.exists(err_name):
        logging.info("decontamination already performed, skipping")
        return clean_name

    # alignment with bwa
    cont_real_link = os.path.realpath(contaminant)
    logging.info("Database real path: %s" % cont_real_link)
    cml = (
        "bwa mem -t %d -R '@RG\\tID:foo\\tSM:bar\\tLB:library1' -T 75 -M %s %s 2> \
    %s | samtools view -h -F 4 - > %s"
        % (n_proc, contaminant, input_reads, err_name, sam_name)
    )
    logging.debug("running bwa %s %s on %d cores", cont_name, rf_head, n_proc)
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

    oh = open("stats.tsv", "a")
    oh.write("matching_%s\t%d\n" % (cont_name, len(mapped_reads)))
    oh.close()

    output_handle = open(clean_name, "w")
    logging.debug(
        "Cleaning reads in %s with alignments in %s", input_reads, sam_name
    )
    logging.debug("Writing to %s", clean_name)
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
                logging.debug("written %d clean reads", c)
    logging.info("written %d clean reads", c)
    output_handle.close()

    if input_reads != "good.fastq":
        os.remove(input_reads)

    return clean_name


def viral_blast(file_in, n_proc, nodes, names):
    """runs blast against viral database, parallelise with xargs"""
    import re
    import sys
    import warnings

    # on hot start, blast again all decontaminated reads
    if os.path.exists("viral_reads.fastq.gz") and os.path.exists(
        "undetermined_reads.fastq.gz"
    ):
        run_child(
            "zcat viral_reads.fastq.gz undetermined_reads.fastq.gz > %s"
            % file_in
        )
        os.remove("viral_reads.fastq.gz")
        os.remove("undetermined_reads.fastq.gz")

    # streams will be used during the execution
    oh = open("stats.tsv", "a")
    bh = open("unique.tsv", "w")
    bh.write(
        "qseqid\tsseqid\tssciname\tstitle\tpident\tqcovs\tscore\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tstaxid\n"
    )

    if not os.path.exists("hq_decont_reads.fastq"):
        os.rename(file_in, "hq_decont_reads.fastq")
    fasta_file = "hq_decont_reads.fasta"
    run_child("seqtk seq -A hq_decont_reads.fastq > %s" % fasta_file)
    try:
        tot_seqs = int(run_child('grep -c "^>" %s' % fasta_file).strip())
    except AttributeError:  # deals with empty file
        tot_seqs = 0
        logging.info("No reads left after decontamination")

    oh.write("reads_to_blast\t%d\n" % tot_seqs)

    if tot_seqs == 0:
        bh.close()
        oh.write("viral_reads\t0\n")
        oh.write("undetermined_reads\t0\n")
        oh.close()
        return

    max_n = (tot_seqs / n_proc) + 1

    # We want to split in n_proc processors, so each file has at most
    # (tot_seqs / n_proc) + 1 reads
    cml = (
        'awk -v "MAX_N=%d" \'BEGIN {n_seq=0;} /^>/ \
    {if(n_seq %% %d == 0){file=sprintf("splitted_clean_%%d.fasta", n_seq/%d);} \
    print >> file; n_seq++; next;} { print >> file; }\' %s'
        % (max_n, max_n, max_n, fasta_file)
    )
    run_child(cml)

    # blast needs access to taxdb files to retrieve organism name
    os.environ["BLASTDB"] = DB_DIR
    if sys.platform.startswith("linux"):
        xargs_thread = 0  # means on all available cores, caution
    elif sys.platform.startswith("darwin"):
        xargs_thread = n_proc  # darwin xargs does not accept -P 0
    else:
        logging.info(
            "could not detect system platform: runnning on %d cores", n_proc
        )
        xargs_thread = n_proc
    # if Darwin then xargs_thread must be n_proc
    DB_real_path = os.path.realpath(
        os.path.join(DB_DIR, "viral_nuccore/viral_db")
    )
    logging.info("Database real path: %s" % DB_real_path)
    cml = (
        "seq 0 %s | xargs -P %d -I {} blastn -task megablast \
           -query splitted_clean_{}.fasta -db %s \
           -out tmp_{}.tsv \
           -outfmt '6 qseqid sseqid ssciname stitle pident qcovs score length mismatch gapopen qstart qend sstart send staxid'"
        % (
            n_proc - 1,
            xargs_thread,
            os.path.join(DB_DIR, "viral_nuccore/viral_db"),
        )
    )
    logging.debug("running blast now")
    run_child(cml)

    logging.debug("saving blast database info")
    cml = shlex.split(
        "blastdbcmd -db /data/virmet_databases/viral_nuccore/viral_db -info"
    )
    with open("blast_info.txt", "a") as boh:
        subprocess.call(cml, stdout=boh)

    logging.debug("parsing best HSP for each query sequence")
    qseqid = ""
    # write to unique.tsv
    for tmpf in glob.glob("tmp_*.tsv"):
        i = tmpf.split("_")[1].split(".")[0]
        with open(tmpf) as f:
            for line in f:
                if line.split("\t")[0] != qseqid:
                    bh.write(line)
                    qseqid = line.split("\t")[0]
        os.remove(tmpf)
        os.remove("splitted_clean_%s.fasta" % i)
    bh.close()

    logging.debug("filtering and grouping by hit sequence")
    hits = pd.read_csv("unique.tsv", index_col="qseqid", delimiter="\t")
    logging.debug("found %d hits", hits.shape[0])
    # select according to identity and coverage, count occurrences
    good_hits = hits[
        (hits.pident > blast_ident_threshold)
        & (hits.qcovs > blast_cov_threshold)
    ]
    matched_reads = good_hits.shape[0]
    logging.debug("%d hits passing coverage and identity filter", matched_reads)
    oh.write("viral_reads\t%s\n" % matched_reads)
    unknown_reads = tot_seqs - matched_reads
    oh.write("undetermined_reads\t%d\n" % unknown_reads)
    oh.close()

    if matched_reads == 0:  # deals with no good_hits
        warnings.warn("No hits")
        return

    # create a column for accession number
    good_hits["accn"] = good_hits.apply(
        lambda row: re.search(r"([A-Z]+_?\d*)\.?\d*", row["sseqid"]).group(1),
        axis=1,
    )
    good_hits = good_hits.rename(columns={"staxid": "tax_id"})

    viral_info_file = os.path.join(DB_DIR, "viral_nuccore/viral_seqs_info.tsv")
    viral_info = pd.read_table(
        viral_info_file,
        names=["accn", "TaxId", "seq_len", "Organism", "Title", "accn_version"],
    )
    viral_info.drop(columns=["accn_version"])
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

    viral_info = viral_info.drop(columns=["TaxId", "Organism", "Title"])
    ds = pd.merge(ds, viral_info)
    # ds['covered_fraction'] = round(ds['covered_region'] / ds['seq_len'], 4)
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
    ds.to_csv("orgs_list.tsv", header=True, sep="\t", index=False)


def cleaning_up():
    """sift reads into viral/unknown, compresses and removes files"""
    import multiprocessing as mp

    from Bio.SeqIO.QualityIO import FastqGeneralIterator

    # selects reads with coverage and identity higher than 75
    df = pd.read_csv("unique.tsv", sep="\t")
    viral_ids = set(
        df[
            (df.qcovs > blast_cov_threshold)
            & (df.pident > blast_ident_threshold)
        ].qseqid
    )
    viral_c = 0
    undet_c = 0
    all_reads = "hq_decont_reads.fastq"
    all_handle = open(all_reads)
    undet_handle = open("undetermined_reads.fastq", "w")
    viral_handle = open("viral_reads.fastq", "w")
    # Using FastqGeneralIterator allows fast performance
    for title, seq, qual in FastqGeneralIterator(all_handle):
        if title.split()[0] not in viral_ids:
            undet_c += 1
            undet_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
            if undet_c % 100000 == 0:
                logging.debug("written %d undet reads", undet_c)
        else:
            viral_c += 1
            viral_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
            if viral_c % 10000 == 0:
                logging.debug("written %d viral reads", viral_c)
    undet_handle.close()
    viral_handle.close()
    logging.info("written %d undet reads", undet_c)
    logging.info("written %d viral reads", viral_c)

    run_child("gzip -f viral_reads.fastq")
    run_child("gzip -f undetermined_reads.fastq")
    os.remove(all_reads)

    cmls = []
    for samfile in glob.glob("*.sam"):
        stem = os.path.splitext(samfile)[0]
        cont = stem.split("_")[-1]
        if cont == "ref":  # hack because _ in bovine file name
            cont = "bt_ref"
        cml = (
            "samtools sort -O bam -l 0 -T /tmp -@ 4 %s | \
        samtools view -T %s -C -o %s.cram -@ 4 -"
            % (samfile, ref_map[cont], stem)
        )
        cmls.append(cml)

    # run in parallel
    pool = mp.Pool()
    results = pool.map(run_child, cmls)
    for r in results:
        logging.debug(r)

    # removing and zipping
    for samfile in glob.glob("*.sam"):
        os.remove(samfile)
    for rf in ["good.fastq", "bad.fastq", "hq_decont_reads.fasta"]:
        try:
            os.remove(rf)
        except FileNotFoundError:
            pass

    for gf in glob.glob("good_*fastq"):
        os.remove(gf)
    run_child("gzip -f unique.tsv")


def main(args):
    """"""

    if args.run:
        miseq_dir = args.run.rstrip("/")
        run_name = os.path.split(miseq_dir)[1]
        if (
            run_name.startswith(("1", "2"))
            and len(run_name.split("-")[-1]) == 5
            and run_name.split("_")[1].startswith("M")
        ):
            try:
                run_date, machine_name = run_name.split("_")[:2]
                logging.info(
                    "running on run %s from machine %s", run_name, machine_name
                )
            except ValueError:
                logging.info("running on directory %s", miseq_dir)
            bc_dir = os.path.join(miseq_dir, "Data/Intensities/BaseCalls/")
        else:
            bc_dir = miseq_dir

        rel_fastq_files = glob.glob("%s/*_S*.fastq*" % bc_dir)
        samples_to_run = [
            os.path.split(fq)[1].split("_")[1] for fq in rel_fastq_files
        ]
        logging.info("samples to run: %s", " ".join(samples_to_run))
        all_fastq_files = [os.path.abspath(f) for f in rel_fastq_files]
    elif args.file:
        all_fastq_files = [os.path.abspath(args.file)]
        logging.info("running on a single file %s", all_fastq_files[0])
        run_name = os.path.split(args.file)[1].split(".")[0]

    out_dir = "virmet_output_%s" % run_name

    try:
        os.mkdir(out_dir)
    except OSError:
        logging.error("directory %s exists", out_dir)
    os.chdir(out_dir)

    # run hunter on all fastq files
    s_dirs = []
    for fq in all_fastq_files:
        logging.info("running hunter on %s", fq)
        sd = hunter(fq)
        s_dirs.append(sd)

    # run mapping against contaminants to remove
    cont_reads = "good.fastq"  # first run on good.fastq
    for cont in contaminant_db:
        logging.info("decontamination against %s", cont)
        for sample_dir in s_dirs:
            logging.info("--- now for sample %s", sample_dir)
            os.chdir(sample_dir)
            decont_reads = victor(input_reads=cont_reads, contaminant=cont)
            os.chdir(os.pardir)
        cont_reads = decont_reads  # decontaminated reads are input for next round (equal across samples)

    logging.info("blasting against viral database")
    file_to_blast = cont_reads  # last output of victor is input for blast
    try:
        n_proc = min(os.cpu_count(), 12)
    except NotImplementedError:
        n_proc = 2
    logging.info("%d cores that will be used", n_proc)

    logging.info("reading taxonomy files")
    nodes, names = get_nodes_names(DB_DIR)

    for sample_dir in s_dirs:
        os.chdir(sample_dir)
        logging.info("now sample %s", sample_dir)
        viral_blast(file_to_blast, n_proc, nodes, names)
        logging.info("sample %s blasted", sample_dir)
        os.chdir(os.pardir)

    logging.info("summarising and cleaning up")
    for sample_dir in s_dirs:
        os.chdir(sample_dir)
        logging.info("now in %s", sample_dir)
        cleaning_up()
        os.chdir(os.pardir)

    os.chdir(os.pardir)
    return out_dir


if __name__ == "__main__":
    import sys

    assert os.path.exists(sys.argv[1])
    all_nodes, all_names = get_nodes_names(DB_DIR)
    viral_blast(sys.argv[1], 2, all_nodes, all_names)
