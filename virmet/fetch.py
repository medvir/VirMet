#!/usr/bin/env python3.4

import os
import sys
import logging
import subprocess
import pandas as pd

from virmet.common import viral_query, bacterial_query, fungal_query, \
ftp_down, run_child, download_genomes, get_gids, DB_DIR


def main(args):
    import random

    logging.info('now in fetch_data')

    if args.viral == 'n':
        logging.info('downloading viral nuccore sequences')
        viral_query('n')
        target_dir = os.path.join(DB_DIR, 'viral_nuccore')
        os.chdir(target_dir)
        cml = '-format fasta < ncbi_search > viral_database.fasta'
        run_child('efetch', cml)
        logging.info('saving viral nuccore taxonomy')
        # viral_seqs_info.tsv contains Gi TaxId
        run_child('cut', '-f 1,2 viral_seqs_info.tsv > viral_gi_taxid.dmp')
        gids_1 = set(get_gids('viral_database.fasta'))
        gids_2 = set([l.split()[0] for l in open('viral_gi_taxid.dmp')])
        assert gids_1 == gids_2

    elif args.viral == 'p':
        logging.info('downloaded viral protein sequences')
        viral_query('p')
        target_dir = os.path.join(DB_DIR, 'viral_protein')
        os.chdir(target_dir)
        cml = ' -format fasta < ncbi_search > viral_database.fasta'
        run_child('efetch', cml)

        logging.info('saving viral protein taxonomy')
        cml = '-format docsum < ncbi_search | xtract -pattern DocumentSummary \
        -element Gi TaxId Caption > tmp.dmp'
        run_child('efetch', cml)
        logging.info('saving viral nuccore taxonomy')
        # viral_seqs_info.tsv contains Gi TaxId
        run_child('cut', '-f 1,2 viral_seqs_info.tsv > viral_gi_taxid.dmp')
        gids_1 = set(get_gids('viral_database.fasta'))
        gids_2 = set([l.split()[0] for l in open('viral_gi_taxid.dmp')])
        assert gids_1 == gids_2

    if args.viral:
        os.chdir(DB_DIR)
        ftp_down('ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz')
        run_child('tar', 'xvfz taxdb.tar.gz')
        os.remove('taxdb.tar.gz')

    if args.bact:
        target_dir = os.path.join(DB_DIR, 'bacteria')
        try:
            os.mkdir(target_dir)
        except FileExistsError:
            pass
        os.chdir(target_dir)

        # first download summary file with all ftp paths and return urls
        all_urls = bacterial_query()
        logging.info('%d bacterial genomes were found' % len(all_urls))
        # then download genomic_fna.gz files
        download_genomes(all_urls, prefix='bact', n_files=3)
        for j in [1, 2, 3]:
            run_child('bgzip', 'fasta/bact%d.fasta' % j)

    elif args.human:
        target_dir = os.path.join(DB_DIR, 'human')
        try:
            os.mkdir(target_dir)
        except FileExistsError:
            pass
        os.chdir(target_dir)
        try:
            os.mkdir('fasta')
        except FileExistsError:
            pass
        os.chdir('fasta')
        fasta_url = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/GRCh38.primary_assembly.genome.fa.gz'
        gtf_url = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.primary_assembly.annotation.gtf.gz'
        logging.info('Downloading human annotation')
        ftp_down(gtf_url)
        logging.info('Downloading human genome and bgzip compressing')
        if os.path.exists('GRCh38.fasta'):
            os.remove('GRCh38.fasta')
        ftp_down(fasta_url, 'GRCh38.fasta')
        run_child('bgzip', 'GRCh38.fasta')


    elif args.fungal:
        target_dir = os.path.join(DB_DIR, 'fungi')
        try:
            os.mkdir(target_dir)
        except FileExistsError:
            pass
        os.chdir(target_dir)

        # first download summary file with all ftp paths and return urls
        all_urls = fungal_query()
        logging.info('%d fungal genomes were found' % len(all_urls))
        # then download genomic_fna.gz files
        download_genomes(all_urls, prefix='fungi', n_files=1)
        run_child('bgzip', 'fasta/fungi1.fasta')

    elif args.bovine:
        target_dir = os.path.join(DB_DIR, 'bovine')
        try:
            os.mkdir(target_dir)
        except FileExistsError:
            pass
        os.chdir(target_dir)
        try:
            os.mkdir('fasta')
        except FileExistsError:
            pass
        os.chdir('fasta')
        chromosomes = ['chr%d' % chrom for chrom in range(1, 30)]
        chromosomes.extend(['chrMT', 'chrX', 'unplaced'])  # Y IS MISSING
        logging.info('Downloading bovine genome')
        local_file_name = os.path.join(target_dir, 'fasta', 'bt_ref_Bos_taurus_UMD_3.1.1.fasta')
        if os.path.exists(local_file_name):
            os.remove(local_file_name)
        for chrom in chromosomes:
            logging.debug('Downloading bovine chromosome %s' % chrom)
            fasta_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/Bos_taurus/Assembled_chromosomes/seq/bt_ref_Bos_taurus_UMD_3.1.1_%s.fa.gz' % chrom
            ftp_down(fasta_url, local_file_name)
        run_child('bgzip', local_file_name)
        logging.info('Downloading gff annotation file')
        gff3_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/Bos_taurus/GFF/ref_Bos_taurus_UMD_3.1.1_top_level.gff3.gz'
        ftp_down(gff3_url)
