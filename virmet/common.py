#!/usr/bin/env python3.4

import os
import random
import logging
import urllib.request
from urllib.request import urlopen, Request
import subprocess
import multiprocessing as mp
import pandas as pd

DB_DIR = '/data/virmet_databases/'
prinseq_exe = '/usr/local/bin/prinseq-lite.pl'


def run_child(exe_name, arg_string, exe='/bin/sh'):
    '''use subrocess.check_output to run an external program with arguments'''
    try:
        output = subprocess.check_output(exe_name + ' ' + arg_string, universal_newlines=True, shell=True, executable=exe)
    except subprocess.CalledProcessError as ee:
        logging.error("Execution of %s failed with returncode %d: %s" % (exe_name, ee.returncode, ee.output))
        logging.error(exe_name + ' ' + arg_string)
        output = None
    return output


def single_process(ec_pair):
    '''single process via run_child, used to parallelize
    '''
    exe, cml = ec_pair
    out = run_child(exe, cml)
    return out


def ftp_down(remote_url, local_url=None):
    '''Handles correctly gzipped and uncompressed files'''
    import gzip
    from io import BytesIO

    if local_url:
        outname = local_url
    else:
        outname = remote_url.split('/')[-1]

    # compressing
    if not remote_url.endswith('.gz') and outname.endswith('.gz'):
        raise NotImplementedError('compressing on the fly not implemented (yet?)')

    # decompressing
    elif remote_url.endswith('.gz') and not outname.endswith('.gz'):
        if os.path.exists(outname):
            outhandle = open(outname, 'a')
        else:
            outhandle = open(outname, 'w')
        with urlopen(Request(remote_url, headers={"Accept-Encoding": "gzip"}), timeout=30) as response, \
                gzip.GzipFile(fileobj=response) as f:
            outhandle.write(f.read().decode('utf-8'))

    # keeping the compression status
    elif remote_url.endswith('.gz') and outname.endswith('.gz'):
        if os.path.exists(outname):
            outhandle = open(outname, 'ab')
        else:
            outhandle = open(outname, 'wb')
        with urllib.request.urlopen(remote_url, timeout=30) as f:
            outhandle.write(f.read())

    # uncompressed to uncompressed
    else:
        if os.path.exists(outname):
            outhandle = open(outname, 'a')
        else:
            outhandle = open(outname, 'w')
        with urllib.request.urlopen(remote_url, timeout=30) as f:
            print(f.read().decode('utf-8'), file=outhandle)

    outhandle.close()


def viral_query(viral_db):
    # Viruses, Taxonomy ID: 10239
    # Human adenovirus A, Taxonomy ID: 129875 (only for testing, 7 hits)
    # Mastadenovirus, Taxonomy ID: 10509 (only for testing, 440 hits)
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
    run_child('esearch', search_text)
    efetch_xtract = '-format docsum < ncbi_search | xtract'
    efetch_xtract += ' -pattern DocumentSummary -element Gi TaxId Caption Slen Organism Title > viral_seqs_info.tsv'
    run_child('efetch', efetch_xtract)
    logging.info('downloaded viral seqs info in %s' % target_dir)


def bacterial_query(download=True, info_file='bacteria_refseq_info.tsv'):
    ''' download bacterial genomes in refseq as explained in FAQ 12 here
    http://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#asmsumfiles
    '''
    if download:
        url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt'
        bh = open(info_file, 'w')
        with urllib.request.urlopen(url) as f:
            print(f.read().decode('utf-8'), file=bh)
        bh.close()
    bactinfo = pd.read_csv(info_file, sep='\t', header=0)
    bactinfo.rename(columns={'# assembly_accession': 'assembly_accession'}, inplace=True)
    gb = bactinfo[(bactinfo.assembly_level == 'Complete Genome') & (bactinfo.version_status == 'latest')]
    gb.set_index('assembly_accession')
    x = gb['ftp_path'].apply(lambda col: col + '/' + col.split('/')[5] + '_genomic.fna.gz')
    gb.loc[:, 'ftp_genome_path'] = pd.Series(x, index=gb.index)
    all_urls = list(gb['ftp_genome_path'])
    assert len(all_urls) == len(gb)
    return all_urls


def fungal_query(download=True, info_file='fungi_refseq_info.tsv'):
    ''' download fungal genomes in refseq in a similar way to bacterial
    '''
    if download:
        url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt'
        bh = open(info_file, 'w')
        with urllib.request.urlopen(url) as f:
            print(f.read().decode('utf-8'), file=bh)
        bh.close()
    funginfo = pd.read_csv(info_file, sep='\t', header=0)
    funginfo.rename(columns={'# assembly_accession': 'assembly_accession'}, inplace=True)
    gb = funginfo[(funginfo.refseq_category != 'na') &
                  (funginfo.version_status == 'latest') &
                  (funginfo.genome_rep == 'Full') &
                  (funginfo.release_type == 'Major')]
    gb.set_index('assembly_accession')
    x = gb['ftp_path'].apply(lambda col: col + '/' + col.split('/')[5] + '_genomic.fna.gz')
    gb.loc[:, 'ftp_genome_path'] = pd.Series(x, index=gb.index)
    all_urls = list(gb['ftp_genome_path'])
    assert len(all_urls) == len(gb)
    return all_urls


def download_genomes(all_urls, prefix, n_files=1):
    ''' download genomes given a list of urls, randomly assigning them
    to one of several (n_files) fasta files
    '''
    # assign sequences randomly to (three) sets using answer here
    # http://stackoverflow.com/questions/2659900/python-slicing-a-list-into-n-nearly-equal-length-partitions
    logging.info('writing %d genome assemblies to fasta files' % len(all_urls))
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


def multiple_download(dl_pair):
    fasta_out, urls = dl_pair
    for j, url in enumerate(urls):
        ftp_down(url, fasta_out)
        if (j + 1) % 500 == 0:
            logging.debug('file %s, %d genomes downloaded' % (fasta_out, j + 1))
    logging.debug('file %s, %d genomes downloaded' % (fasta_out, j + 1))
    return


def get_gids(fasta_file):
    if fasta_file.endswith('.gz'):
        cml = '%s | grep \"^>\" | cut -f 2 -d \"|\"' % fasta_file
        gids = run_child('zcat', cml).strip().split('\n')
    else:
        cml = '\"^>\" %s | cut -f 2 -d \"|\"' % fasta_file
        gids = run_child('grep', cml).strip().split('\n')
    return gids
