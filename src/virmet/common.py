#!/usr/bin/env python3

import os
import random
import logging
from time import sleep
import urllib.request
from urllib.request import urlopen, Request
import subprocess
import multiprocessing as mp
import pandas as pd

DB_DIR = '/data/virmet_databases/'


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
                    if "No such file or directory" in e.reason and "550" in e.reason:
                        msg = "No remote file found (some ffn and faa are know to be missing remotely). Aborting the download of %s. %s" % (
                            str(args[0]), str(e))
                        logging.warning(msg)
                        raise e
                    else:
                        msg = "%s: %s, Retrying in %d seconds..." % (
                            str(args[0]), str(e), mdelay)
                        logging.warning(msg)
                sleep(mdelay)
                mtries -= 1  # consume an attempt
                mdelay *= backoff  # make future wait longer
            return f(*args, **kwargs)  # Ran out of tries :-(
        return f_retry  # true decorator -> decorated function
    return deco_retry  # @retry(arg[, ...]) -> true decorator


def run_child(cmd, exe='/bin/bash'):
    """Use subrocess.check_output to run an external program with arguments."""
    logging.info('Running instance of %s' % cmd.split()[0])
    try:
        output = subprocess.check_output(cmd, universal_newlines=True, shell=True, stderr=subprocess.STDOUT)
        logging.debug('Completed')
    except subprocess.CalledProcessError as ee:
        logging.error("Execution of %s failed with returncode %d: %s", cmd, ee.returncode, ee.output)
        logging.error(cmd)
        output = None
    return output


@retry(tries=8, delay=5, backoff=1.5)
def ftp_down(remote_url, local_url=None):
    """Download files, correctly handling both gzipped and uncompressed files."""
    import gzip
    # from io import BytesIO

    if local_url:
        outname = local_url
    else:
        outname = remote_url.split('/')[-1]
    logging.debug('Downloading %s', remote_url)
    # compressing
    if not remote_url.endswith('.gz') and outname.endswith('.gz'):
        raise NotImplementedError('compressing on the fly not implemented (yet?)')

    # decompressing
    elif remote_url.endswith('.gz') and not outname.endswith('.gz'):
        if os.path.exists(outname):
            outhandle = open(outname, 'a')
        else:
            outhandle = open(outname, 'w')
        logging.debug('Downloading %s', remote_url)
        with urlopen(Request(remote_url, headers={"Accept-Encoding": "gzip"}), timeout=30) as response, \
                gzip.GzipFile(fileobj=response) as f:
            outhandle.write(f.read().decode('utf-8'))

    # keeping the compression status
    elif remote_url.endswith('.gz') and outname.endswith('.gz'):
        if os.path.exists(outname):
            outhandle = open(outname, 'ab')
        else:
            outhandle = open(outname, 'wb')
        logging.debug('Downloading %s', remote_url)
        with urllib.request.urlopen(remote_url, timeout=30) as f:
            outhandle.write(f.read())

    # uncompressed to uncompressed
    else:
        if os.path.exists(outname):
            outhandle = open(outname, 'a')
        else:
            outhandle = open(outname, 'w')
        logging.debug('Downloading %s', remote_url)
        with urllib.request.urlopen(remote_url, timeout=30) as f:
            # print(f.read().decode('utf-8', 'replace').encode('utf-8', 'replace'), file=outhandle)
            outhandle.write(f.read().decode('utf-8', 'replace'))  # .encode('utf-8', 'replace'), file=outhandle)

    # outhandle.close()
    return outhandle


def viral_query(viral_db):
    # Viruses, Taxonomy ID: 10239
    # Human adenovirus A, Taxonomy ID: 129875 (only for testing, 7 hits)
    # Mastadenovirus, Taxonomy ID: 10509 (only for testing, 440 hits)
    # Alphatorquevirus Taxonomy ID: 687331
    # Cellular organisms, Taxonomy ID: 131567 (to avoid chimeras)
    txid = '10239'  # change here for viruses or smaller taxa
    if viral_db == 'n':
        target_dir = os.path.join(DB_DIR, 'viral_nuccore')
        search_text = '-db nuccore -query \"txid%s [orgn] AND \\"complete genome\\" [Title] NOT txid131567 [orgn]\" > ncbi_search' % txid
    elif viral_db == 'p':
        target_dir = os.path.join(DB_DIR, 'viral_protein')
        search_text = '-db protein -query \"refseq [filter] AND txid%s [orgn] NOT txid131567 [orgn]\" > ncbi_search' % txid
    try:
        os.mkdir(target_dir)
    except FileExistsError:
        pass
    os.chdir(target_dir)
    return 'esearch ' + search_text


def bact_fung_query(query_type=None, download=True, info_file=None):
    """Download/read bacterial and fungal genomes in refseq as explained in
    FAQ 12 here http://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#asmsumfiles.

    If info_file is not given, it will be inferred from query_type;
    if download is true, it will be used to save the downloaded info;
    if download is false, an already present file will be read.
    """
    if query_type not in ['bacteria', 'fungi']:
        raise SystemExit('Choose bacteria or fungi')
    if not info_file:
        info_file = '%s_refseq_info.tsv' % query_type
    if download:
        url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/%s/assembly_summary.txt' % query_type
        bh = open(info_file, 'w')
        with urllib.request.urlopen(url) as f:
            print(f.read().decode('utf-8'), file=bh)
        bh.close()
    querinfo = pd.read_csv(info_file, sep='\t', header=0, skiprows=1)
    querinfo.rename(columns={'# assembly_accession': 'assembly_accession'}, inplace=True)
    if query_type == 'bacteria':
        gb = querinfo[(querinfo.assembly_level == 'Complete Genome') &
                      (querinfo.version_status == 'latest')]
    elif query_type == 'fungi':
        gb = querinfo[((querinfo.assembly_level == 'Complete Genome') | (querinfo.assembly_level == 'Chromosome')) &
                      (querinfo.refseq_category != 'na') &
                      (querinfo.version_status == 'latest') &
                      (querinfo.genome_rep == 'Full') &
                      (querinfo.release_type == 'Major')]
    gb.set_index('assembly_accession')
    x = gb['ftp_path'].apply(lambda col: col + '/' + col.split('/')[-1] + '_genomic.fna.gz')
    gb = gb.assign(ftp_genome_path=x)
    all_urls = list(gb['ftp_genome_path'])
    assert len(all_urls) == len(gb)
    return all_urls


def download_genomes(all_urls, prefix, n_files=1):
    """Download genomes given a list of urls, randomly assigning them to one of several (n_files) fasta files.

    It assigns sequences randomly to (three) sets using answer here
    http://stackoverflow.com/questions/2659900/python-slicing-a-list-into-n-nearly-equal-length-partitions
    """
    logging.info('writing %d genome assemblies to fasta files', len(all_urls))
    random.shuffle(all_urls)
    q, r = divmod(len(all_urls), n_files)  # quotient, remainder
    indices = [q * i + min(i, r) for i in range(n_files + 1)]
    seqs_urls = [all_urls[indices[i]:indices[i + 1]] for i in range(n_files)]
    try:
        os.mkdir('fasta')
    except FileExistsError:
        pass

    dl_pairs = []
    for i, seqs in enumerate(seqs_urls):
        fasta_out = 'fasta/%s%d.fasta' % (prefix, i + 1)
        if os.path.exists(fasta_out):
            os.remove(fasta_out)
        dl_pairs.append((fasta_out, seqs))

    # run download in parallel
    pool = mp.Pool()
    results = pool.map(multiple_download, dl_pairs)
    return results


def multiple_download(dl_pair):
    fasta_out, urls = dl_pair
    logging.debug('About to download %d files', len(urls))
    for url in urls:
        downloaded_handle = ftp_down(url, fasta_out)
        downloaded_handle.close()
    return len(urls)


def get_gids(fasta_file):
    if fasta_file.endswith('.gz'):
        cml = 'zcat %s | grep \"^>\" | cut -f 2 -d \"|\"' % fasta_file
        gids = run_child(cml).strip().split('\n')
    else:
        cml = 'grep \"^>\" %s | cut -f 2 -d \"|\"' % fasta_file
        gids = run_child(cml).strip().split('\n')
    return gids


def get_accs(fasta_file):
    if fasta_file.endswith('.gz'):
        cml = 'zcat %s | grep \"^>\" | tr -d \">\" | cut -f 1 -d \".\"' % fasta_file
        accs = run_child(cml).strip().split('\n')
    else:
        cml = 'grep \"^>\" %s | tr -d \">\" | cut -f 1 -d \".\"' % fasta_file
        accs = run_child(cml).strip().split('\n')
    return accs
