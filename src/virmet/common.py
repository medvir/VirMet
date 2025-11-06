#!/usr/bin/env python3

import gzip
import logging
import os
import random
import re
import requests
import subprocess
from Bio import Entrez
from datetime import datetime, timedelta
from time import sleep
from urllib.request import Request, urlopen

import pandas as pd
import numpy as np

MAX_TAXID = 10000  # max number of sequences belonging to the same taxid in compressed viral database
n_proc = int((os.cpu_count() or 8) * 0.75)


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


def safe_entrez_read(**kwargs):
    retries = 5
    for attempt in range(retries):
        try:
            with Entrez.esearch(**kwargs) as handle:
                return Entrez.read(handle)
        except RuntimeError as e:
            logging.warning(f"Entrez failed: {e}")
            if attempt < retries - 1:
                sleep(3)
            else:
                raise


def run_child(cmd):
    """Use subrocess.check_output to run an external program with arguments."""
    logging.info("Running instance of %s" % cmd.split()[0])
    try:
        output = subprocess.check_output(
            cmd,
            universal_newlines=True,
            shell=True,  # nosec B602: Required for shell piping.
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
    # compressing
    if not remote_url.endswith(".gz") and outname.endswith(".gz"):
        raise NotImplementedError(
            "compressing on the fly not implemented (yet?)"
        )

    # decompressing
    elif remote_url.endswith(".gz") and not outname.endswith(".gz"):
        if os.path.exists(outname):
            outhandle = open(outname, "a")
        else:
            outhandle = open(outname, "w")
        logging.debug("Downloading %s", remote_url)
        with (
            urlopen(  # nosec B310: opening url
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
        with urlopen(remote_url, timeout=300) as f:  # nosec B310: opening url
            outhandle.write(f.read())

    # uncompressed to uncompressed
    else:
        if os.path.exists(outname):
            outhandle = open(outname, "a")
        else:
            outhandle = open(outname, "w")
        logging.debug("Downloading %s", remote_url)
        with urlopen(remote_url, timeout=30) as f:  # nosec B310: opening url
            outhandle.write(f.read().decode("utf-8", "replace"))

    outhandle.close()


def find_file(base_url, file_regex):
    """
    Given a base directory URL, find the file matching the regex pattern.

    Parameters:
    - base_url: str, the URL to the directory to scan.
    - file_regex: regex pattern to match file names
    """
    logging.info(f"Scanning {base_url}")
    response = requests.get(base_url, timeout=10000)  # nosec B310: opening url
    response.raise_for_status()

    # Look for file matches
    match = list(set(re.findall(file_regex, response.text)))[0]
    if not match:
        raise ValueError(
            f"No matching files found at {base_url} with pattern {file_regex}"
        )

    logging.info(f"Found file: {match}.")

    return match


def random_reduction(viral_mode, DB_DIR_UPDATE):
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
        "Database real path for compression: %s" % os.path.realpath(target_dir)
    )

    viral_info_file = os.path.join(target_dir, "viral_seqs_info.tsv")
    viral_fasta_file = os.path.join(target_dir, "viral_database.fasta")

    viral_info = pd.read_table(
        viral_info_file, names=["accn_version", "TaxId", "seq_len", "Organism"]
    )
    viral_info.drop_duplicates(subset=["accn_version"], inplace=True)
    # Do compression at the texid ID level,
    # check if any taxid has more than 10,000 reference sequences move them to a new pandas table
    # remove them from the main pandas table
    TaxId_to_counter_df = viral_info["TaxId"].value_counts().reset_index()
    TaxId_to_counter_df.rename(
        columns={"TaxId": "TaxId_num", "count": "TaxId_count"}, inplace=True
    )
    TaxId_to_percentage = (
        viral_info["TaxId"].value_counts(normalize=True).reset_index()
    )
    TaxId_to_counter_df["percentage"] = TaxId_to_percentage["proportion"]
    TaxId_to_counter_filterred_df = TaxId_to_counter_df[
        TaxId_to_counter_df["TaxId_count"] > MAX_TAXID
    ]

    viral_info_subsampled = viral_info.copy()
    random.seed(100)
    for _, row in TaxId_to_counter_filterred_df.iterrows():
        taxid_to_subsample = int(row["TaxId_num"])
        viral_info_to_subsample_df = viral_info_subsampled[
            viral_info_subsampled["TaxId"] == taxid_to_subsample
        ]
        accn_set = set(viral_info_to_subsample_df["accn_version"])
        # subsample to make it about 1% of the database
        selected_accn_ls = random.sample(list(accn_set), MAX_TAXID)  # nosec B311: escape bandit

        # filter out the unselected IDs from viral_info
        # (viral_info_subsampled['accn'] in selected_accn_ls)
        # Keep refseq sequence: all accn of refseq sequences contain char  '_'
        viral_info_subsampled = viral_info_subsampled.loc[
            (viral_info_subsampled["TaxId"] != taxid_to_subsample)
            | viral_info_subsampled["accn_version"].isin(selected_accn_ls)
            | viral_info_subsampled["accn_version"].str.contains("_")
        ]

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


def viral_query(DB_DIR_UPDATE, viral_db, update_min_date=None):
    # Viruses, Taxonomy ID: 10239
    # Human adenovirus A, Taxonomy ID: 129875 (only for testing, 7 hits)
    # Mastadenovirus, Taxonomy ID: 10509 (only for testing, 440 hits)
    # Alphatorquevirus Taxonomy ID: 687331
    # Cellular organisms, Taxonomy ID: 131567 (to avoid chimeras)
    txid = "10239"  # change here for viruses or smaller taxa
    query_text = (
        'txid%s [orgn] AND \
        ("complete genome" [Title] OR \
        "complete segment" [Title] OR \
        srcdb_refseq[prop])'
        % txid
    )
    query_text += ' NOT "cellular organisms"[Organism] NOT \
        AC_000001[PACC] : AC_999999[PACC]'

    if update_min_date:
        logging.info(
            "Viral Database Update is performed with sequences added to NCBI after %s .\n"
            % update_min_date,
        )
        min_date_chosen = update_min_date
    else:
        min_date_chosen = "1980/01/01"

    if viral_db == "n":
        target_dir = os.path.join(DB_DIR_UPDATE, "viral_nuccore")
        db_text = "nuccore"
    elif viral_db == "p":
        target_dir = os.path.join(DB_DIR_UPDATE, "viral_protein")
        db_text = "protein"
    else:
        raise ValueError(f"Invalid viral_db value: '{viral_db}'.")

    # Create output folder if it doesn't exist
    os.makedirs(target_dir, exist_ok=True)
    logging.info("Database real path: %s" % os.path.realpath(target_dir))

    # Obtain Accessions until 2 days ago so all databases should be updated & match
    Entrez.email = "virmet@virmet.ch"
    try:
        Entrez.api_key = os.environ["NCBI_API_KEY"]
    except KeyError:
        logging.info("Please, export a NCBI_API_KEY in your .bashrc")
    handle = Entrez.esearch(
        db=db_text,
        idtype="acc",
        term=query_text,
        datetype="pdat",
        mindate=min_date_chosen,
        maxdate=str((datetime.today() - timedelta(2)).strftime("%Y/%m/%d")),
        retstart=0,
        retmax=1,
    )
    record = Entrez.read(handle)
    tot_accs = int(record["Count"])
    batch_size = 6000

    # Retreive Accession numbers
    ncbi_accs = np.empty(tot_accs, dtype="U16")
    for i in range(0, tot_accs, batch_size):
        sleep(2)
        record = safe_entrez_read(
            db=db_text,
            idtype="acc",
            term=query_text,
            datetype="pdat",
            mindate=min_date_chosen,
            maxdate=str((datetime.today() - timedelta(2)).strftime("%Y/%m/%d")),
            retstart=i,
            retmax=batch_size,
        )
        id_list = np.array(record["IdList"])
        ncbi_accs[i : i + len(id_list)] = id_list
    return ncbi_accs


def get_accs(fasta_file):
    cml = f"seqkit seq -i -n {fasta_file} --threads {n_proc}"
    accs = run_child(cml).strip().split("\n")
    return accs
