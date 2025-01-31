#!/usr/bin/env python3

import gzip
import logging
import multiprocessing as mp
import os
import random
import subprocess
import urllib.request
from time import sleep
from urllib.request import Request, urlopen

import pandas as pd

# TODO: This should be updated to a more global location rather than user based.
DB_DIR = os.path.expandvars("/data/virmet_databases")
DB_DIR_UPDATE = os.path.expandvars("/data/virmet_databases_update")
N_FILES_BACT = 5
MAX_TAXID = 10000  # max number of sequences belonging to the same taxid in compressed viral database
n_proc = min(os.cpu_count() or 8, 16)


# decorator for taken from RepoPhlan
# https://bitbucket.org/nsegata/repophlan/src/5804db9d341165f72c4f7e448691ce90a1465764/repophlan_get_microbes.py?at=ncbi2017&fileviewer=file-view-default
def retry(tries, delay=3, backoff=2):
    """Decorator taken from RepoPhlan, used to retry failed download from NCBI, for example."""
    if backoff <= 1:
        raise ValueError("backoff must be greater than 1")

    # tries = math.floor(tries)
    if tries < 0:
        raise ValueError("tries must be 0 or greater")

    if delay <= 0:
        raise ValueError("delay must be greater than 0")

    def deco_retry(f):
        def f_retry(*args, **kwargs):
            mtries, mdelay = tries, delay  # make mutable
            # rv = f(*args, **kwargs) # first attempt
            while mtries > 1:
                try:
                    return f(*args, **kwargs)
                except Exception as e:
                    if "No such file or directory" in str(e) and "550" in str(
                        e
                    ):
                        logging.warning(
                            "No remote file found (some ffn and faa are know to be missing remotely). Aborting the download of %s. %s",
                            str(args[0]),
                            str(e),
                        )
                        raise e
                    else:
                        logging.warning(
                            "%s: %s, Retrying in %d seconds...",
                            str(args[0]),
                            str(e),
                            mdelay,
                        )
                sleep(mdelay)
                mtries -= 1  # consume an attempt
                mdelay *= backoff  # make future wait longer
            return f(*args, **kwargs)  # Ran out of tries :-(

        return f_retry  # true decorator -> decorated function

    return deco_retry  # @retry(arg[, ...]) -> true decorator


def run_child(cmd):
    """Use subrocess.check_output to run an external program with arguments."""
    logging.info("Running instance of %s", cmd.split()[0])
    try:
        output = subprocess.check_output(
            cmd,
            universal_newlines=True,
            shell=True,
            stderr=subprocess.STDOUT,
            executable="/bin/bash",
        )
        logging.debug("Completed")
    except subprocess.CalledProcessError as ee:
        logging.error(
            "Execution of %s failed with returncode %d: %s",
            cmd,
            ee.returncode,
            ee.output,
        )
        logging.error(cmd)
        output = None
    return output


@retry(tries=8, delay=5, backoff=1.5)
def ftp_down(remote_url, local_url=None):
    """Download files, correctly handling both gzipped and uncompressed files."""

    # from io import BytesIO

    if local_url:
        outname = local_url
    else:
        outname = remote_url.split("/")[-1]
    logging.debug("Downloading %s", remote_url)
    # compressing
    if not remote_url.endswith(".gz") and outname.endswith(".gz"):
        raise NotImplementedError(
            "compressing on the fly not implemented (yet?)"
        )

    # decompressing
    elif remote_url.endswith(".gz") and not outname.endswith(".gz"):
        # print(remote_url)
        # print(outname)
        if os.path.exists(outname):
            outhandle = open(outname, "a")
        else:
            outhandle = open(outname, "w")
        logging.debug("Downloading %s", remote_url)
        with (
            urlopen(
                Request(remote_url, headers={"Accept-Encoding": "gzip"}),
                timeout=30,
            ) as response,
            gzip.GzipFile(fileobj=response) as f,
        ):
            outhandle.write(f.read().decode("utf-8"))

    # keeping the compression status
    elif remote_url.endswith(".gz") and outname.endswith(".gz"):
        if os.path.exists(outname):
            outhandle = open(outname, "ab")
        else:
            outhandle = open(outname, "wb")
        logging.debug("Downloading %s", remote_url)
        with urllib.request.urlopen(remote_url, timeout=300) as f:
            outhandle.write(f.read())

    # uncompressed to uncompressed
    else:
        if os.path.exists(outname):
            outhandle = open(outname, "a")
        else:
            outhandle = open(outname, "w")
        logging.debug("Downloading %s", remote_url)
        with urllib.request.urlopen(remote_url, timeout=30) as f:
            # print(f.read().decode('utf-8', 'replace').encode('utf-8', 'replace'), file=outhandle)
            outhandle.write(
                f.read().decode("utf-8", "replace")
            )  # .encode('utf-8', 'replace'), file=outhandle)

    # outhandle.close()
    return outhandle


def random_reduction(viral_mode):
    # This code, identify sequences from the same species using their taxid.
    # Based on the given thresholds (10,000) subsample sequences of the corresponding taxid.
    # Currently this code has the absolute threshold and any taxid with above
    # 10,000 reference sequences are filtered out.
    # This is mainly done because of SARS-CoV-2 sequences covering above 90% of
    # the viral database as of begining of 2023.
    # The highest frequent virus is SARS-CoV-2 with 90% sequences of the database.

    if viral_mode == "n":
        logging.info("compressing viral nuccore sequences")
        target_dir = os.path.join(DB_DIR_UPDATE, "viral_nuccore")
    elif viral_mode == "p":
        logging.info("compressing viral protein sequences")
        target_dir = os.path.join(DB_DIR_UPDATE, "viral_protein")
    else:
        raise ValueError('Invalid viral_mode: "{viral_mode}".')
    logging.info(
        "Database real path for compression: %s", os.path.realpath(target_dir)
    )

    viral_info_file = os.path.join(target_dir, "viral_seqs_info.tsv")
    viral_fasta_file = os.path.join(target_dir, "viral_database.fasta")

    viral_info = pd.read_table(
        viral_info_file,
        names=["accn", "TaxId", "seq_len", "Organism", "Title", "accn_version"],
    )
    viral_info.drop_duplicates(subset=["accn_version"], inplace=True)
    # Do compression at the texid ID level,
    # check if any taxid has more than 10,000 reference sequences move them to a new pandas table
    # remove them from the main pandas table
    TaxId_to_counter_df = viral_info["TaxId"].value_counts().reset_index()
    TaxId_to_counter_df.rename(
        columns={"TaxId": "TaxId_count", "index": "TaxId_num"}, inplace=True
    )
    TaxId_to_percentage = (
        viral_info["TaxId"].value_counts(normalize=True).reset_index()
    )
    TaxId_to_counter_df["percentage"] = TaxId_to_percentage["TaxId"]
    TaxId_to_counter_filterred_df = TaxId_to_counter_df[
        TaxId_to_counter_df["TaxId_count"] > MAX_TAXID
    ]

    viral_info_subsampled = viral_info.copy()
    random.seed(100)
    removed_set = set()
    for _, row in TaxId_to_counter_filterred_df.iterrows():
        taxid_to_subsample = int(row["TaxId_num"])
        viral_info_to_subsample_df = viral_info_subsampled[
            viral_info_subsampled["TaxId"] == taxid_to_subsample
        ]
        accn_set = set(viral_info_to_subsample_df["accn_version"])
        # subsample to make it about 1% of the database
        selected_accn_ls = random.sample(accn_set, MAX_TAXID)

        # filter out the unselected IDs from viral_info
        # (viral_info_subsampled['accn'] in selected_accn_ls)
        # Keep refseq sequence: all accn of refseq sequences contain char  '_'
        # viral_info_subsampled['accn'].contain('_')  viral_info_subsampled['accn'].str.contains('ON8')
        viral_info_subsampled = viral_info_subsampled.loc[
            (viral_info_subsampled["TaxId"] != taxid_to_subsample)
            | viral_info_subsampled["accn_version"].isin(selected_accn_ls)
            | viral_info_subsampled["accn_version"].str.contains("_")
        ]
        # create a text file of the non_selected ACC ID as extra information
        removed_set = removed_set | (accn_set - set(selected_accn_ls))

    viral_info_subsampled.drop_duplicates(subset=["accn_version"], inplace=True)
    viral_info_subsampled.accn_version.to_csv(
        os.path.join(target_dir, "outfile.csv"),
        sep="\n",
        index=False,
        header=False,
    )
    # extract the selected accession number from the fasta file using seqtk
    subsample_fasta_command = (
        "seqtk subseq %s/viral_database.fasta %s/outfile.csv >  %s/viral_database_subsampled.fasta"
        % (target_dir, target_dir, target_dir)
    )
    run_child(subsample_fasta_command)
    os.rename(
        viral_fasta_file,
        os.path.join(target_dir, "viral_database_original_rmdup.fasta"),
    )
    os.rename(
        os.path.join(target_dir, "viral_database_subsampled.fasta"),
        viral_fasta_file,
    )
    TaxId_to_counter_filterred_df.to_csv(
        os.path.join(target_dir, "filtered_taxids.csv"), sep=",", index=False
    )


def viral_query(viral_db, update_min_date=None):
    # Viruses, Taxonomy ID: 10239
    # Human adenovirus A, Taxonomy ID: 129875 (only for testing, 7 hits)
    # Mastadenovirus, Taxonomy ID: 10509 (only for testing, 440 hits)
    # Alphatorquevirus Taxonomy ID: 687331
    # Cellular organisms, Taxonomy ID: 131567 (to avoid chimeras)
    txid = "10239"  # change here for viruses or smaller taxa
    query_text = (
        f'-query "txid{txid} [orgn] AND ('
        '\\"complete genome\\" [Title] OR '
        '\\"complete segment\\" [Title] OR '
        "srcdb_refseq[prop])"
    )
    query_text += ' NOT \\"cellular organisms\\"[Organism] NOT AC_000001[PACC] : AC_999999[PACC]"'

    if update_min_date:
        logging.info(
            "Viral Database Update is performed with sequences added to NCBI after %s .\n",
            update_min_date,
        )
        query_text += "-datetype PDAT -mindate %s" % str(
            update_min_date
        )  # -datetype MDAT -mindate %s'

    query_text += " > ncbi_search"

    if viral_db == "n":
        target_dir = os.path.join(DB_DIR_UPDATE, "viral_nuccore")
        search_text = "-db nuccore " + query_text
    elif viral_db == "p":
        target_dir = os.path.join(DB_DIR_UPDATE, "viral_protein")
        search_text = "-db protein " + query_text
    else:
        raise ValueError(f"Invalid viral_db value: '{viral_db}'.")

    os.makedirs(target_dir, exist_ok=True)
    logging.info("Database real path: %s", os.path.realpath(target_dir))
    return "esearch " + search_text


def bact_fung_query(
    query_type=None, download=True, info_file=None, target_folder="./"
):
    """Download/read bacterial and fungal genomes in refseq as explained in
    FAQ 12 here http://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#asmsumfiles.

    If info_file is not given, it will be inferred from query_type;
    if download is true, it will be used to save the downloaded info;
    if download is false, an already present file will be read.
    """
    if query_type not in ["bacteria", "fungi"]:
        raise SystemExit("Choose bacteria or fungi")
    if not info_file:
        info_file = "%s_refseq_info.tsv" % query_type
    info_file = os.path.join(target_folder, info_file)
    if download:
        url = (
            "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/%s/assembly_summary.txt"
            % query_type
        )
        bh = open(info_file, "w")
        with urllib.request.urlopen(url) as f:
            print(f.read().decode("utf-8"), file=bh)
        bh.close()
    querinfo = pd.read_csv(
        info_file,
        sep="\t",
        header=0,
        skiprows=1,
        na_values = "na",
        dtype={"excluded_from_refseq": str, "relation_to_type_material": str, "pubmed_id": str, "refseq_category": str} 
    )
    querinfo.rename(
        columns={"#assembly_accession": "assembly_accession"}, inplace=True
    )
    if query_type == "bacteria":
        gb = querinfo[
            (querinfo.assembly_level == "Complete Genome")
            & (querinfo.version_status == "latest")
        ]
    elif query_type == "fungi":
        gb = querinfo[
            (
                (querinfo.assembly_level == "Complete Genome")
                | (querinfo.assembly_level == "Chromosome")
            )
            & (pd.notna(querinfo.refseq_category))
            & (querinfo.version_status == "latest")
            & (querinfo.genome_rep == "Full")
            & (querinfo.release_type == "Major")
        ]
    else:
        raise ValueError(f"Invalid query_type value: '{query_type}'.")

    gb.set_index("assembly_accession", inplace=True)
    gb = gb[gb['ftp_path'].notna()]
    x = gb["ftp_path"].apply(
        lambda col: col + "/" + col.split("/")[-1] + "_genomic.fna.gz"
    )
    gb = gb.assign(ftp_genome_path=x)
    all_urls = list(gb["ftp_genome_path"])
    assert len(all_urls) == len(gb)
    return all_urls


def download_genomes(all_urls, prefix, n_files=1):
    """Download genomes given a list of urls, randomly assigning them to one of several (n_files) fasta files.

    It assigns sequences randomly to (three) sets using answer here
    http://stackoverflow.com/questions/2659900/python-slicing-a-list-into-n-nearly-equal-length-partitions
    """
    logging.info("writing %d genome assemblies to fasta files", len(all_urls))
    random.shuffle(all_urls)
    q, r = divmod(len(all_urls), n_files)  # quotient, remainder
    indices = [q * i + min(i, r) for i in range(n_files + 1)]
    seqs_urls = [all_urls[indices[i] : indices[i + 1]] for i in range(n_files)]
    os.makedirs("fasta", exist_ok=True)

    dl_pairs = []
    for i, seqs in enumerate(seqs_urls):
        fasta_out = f"fasta/{prefix}{i + 1}.fasta"
        # if os.path.exists(fasta_out):
        #    os.remove(fasta_out)
        dl_pairs.append((fasta_out, seqs))

    # run download in parallel
    with mp.Pool() as pool:
        results = pool.map(multiple_download, dl_pairs)
    return results


def multiple_download(dl_pair):
    fasta_out, urls = dl_pair
    logging.debug("About to download %d files", len(urls))
    for url in urls:
        downloaded_handle = ftp_down(url, fasta_out)
        downloaded_handle.close()
    return len(urls)


def get_gids(fasta_file):
    if fasta_file.endswith(".gz"):
        cml = f'zcat {fasta_file} | grep "^>" | cut -f 2 -d "|"'
        gids = run_child(cml).strip().split("\n")
    else:
        cml = f'grep "^>" {fasta_file} | cut -f 2 -d "|"'
        gids = run_child(cml).strip().split("\n")
    return gids


def get_accs(fasta_file):
    if fasta_file.endswith(".gz"):
        cml = f'zcat {fasta_file} | grep "^>" | tr -d ">" | cut -f 1 -d "."'
        accs = run_child(cml).strip().split("\n")
    else:
        cml = f'grep "^>" {fasta_file} | tr -d ">" | cut -f 1 -d "."'
        accs = run_child(cml).strip().split("\n")
    return accs
