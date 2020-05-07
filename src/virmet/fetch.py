#!/usr/bin/env python3
"""Define all the functions that download reference sequences and other files."""
import os
import logging

from virmet.common import viral_query, bact_fung_query, ftp_down, run_child, \
download_genomes, get_accs, DB_DIR, N_FILES_BACT


def fetch_viral(viral_mode):
    """Download nucleotide or protein database."""
    # define the search nuccore/protein
    if viral_mode == 'n':
        logging.info('downloading viral nuccore sequences')
        target_dir = os.path.join(DB_DIR, 'viral_nuccore')
        cml_search = viral_query('n')
    elif viral_mode == 'p':
        logging.info('downloaded viral protein sequences')
        target_dir = os.path.join(DB_DIR, 'viral_protein')
        cml_search = viral_query('p')
    # run the search and download
    os.chdir(target_dir)
    run_child(cml_search)
    cml_fetch_fasta = 'efetch -format fasta < ncbi_search > viral_database.fasta'
    run_child(cml_fetch_fasta)
    cml_efetch_xtract = 'efetch -format docsum < ncbi_search | xtract'
    cml_efetch_xtract += ' -pattern DocumentSummary -element Caption TaxId Slen Organism Title > viral_seqs_info.tsv'
    run_child(cml_efetch_xtract)
    logging.info('downloaded viral seqs info in %s', target_dir)
    logging.info('saving viral taxonomy')
    # viral_seqs_info.tsv contains Accn TaxId
    cml = 'cut -f 1,2 viral_seqs_info.tsv > viral_accn_taxid.dmp'
    run_child(cml)
    accs_1 = set(get_accs('viral_database.fasta'))
    accs_2 = set([l.split()[0] for l in open('viral_accn_taxid.dmp')])
    assert accs_1 == accs_2, accs_1 ^ accs_2
    logging.info('taxonomy and fasta sequences match')

    os.chdir(DB_DIR)
    logging.info('downloading taxonomy databases')
    download_handle = ftp_down('ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz')
    download_handle.close()
    run_child('tar xvfz taxdb.tar.gz')
    os.remove('taxdb.tar.gz')
    download_handle = ftp_down('ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz')
    download_handle.close()
    run_child('tar xvfz taxdump.tar.gz')
    for ftd in ['taxdump.tar.gz', 'merged.dmp', 'gencode.dmp', 'division.dmp', 'delnodes.dmp', 'citations.dmp']:
        try:
            os.remove(ftd)
        except OSError:
            logging.warning('Could not find file %s', ftd)


def fetch_bacterial():
    """Download the three bacterial sequence databases."""
    target_dir = os.path.join(DB_DIR, 'bacteria')
    try:
        os.mkdir(target_dir)
    except FileExistsError:
        pass
    os.chdir(target_dir)

    # first download summary file with all ftp paths and return urls
    all_urls = bact_fung_query(query_type='bacteria')
    logging.info('%d bacterial genomes were found', len(all_urls))
    # then download genomic_fna.gz files
    download_genomes(all_urls, prefix='bact', n_files=N_FILES_BACT)
    for j in range(1, N_FILES_BACT+1):
        run_child('bgzip fasta/bact%d.fasta' % j)


def fetch_human():
    """Download human genome and annotations."""
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
    fasta_url = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz'
    gtf_url = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.primary_assembly.annotation.gtf.gz'    
    logging.info('Downloading human annotation')
    download_handle = ftp_down(gtf_url)
    download_handle.close()
    logging.info('Downloading human genome and bgzip compressing')
    if os.path.exists('GRCh38.fasta'):
        os.remove('GRCh38.fasta')
    download_handle = ftp_down(fasta_url, 'GRCh38.fasta')
    download_handle.close()
    run_child('bgzip GRCh38.fasta')


def fetch_fungal():
    """Download fungal sequences."""
    target_dir = os.path.join(DB_DIR, 'fungi')
    try:
        os.mkdir(target_dir)
    except FileExistsError:
        pass
    os.chdir(target_dir)

    # first download summary file with all ftp paths and return urls
    all_urls = bact_fung_query(query_type='fungi')
    logging.info('%d fungal genomes were found', len(all_urls))
    # then download genomic_fna.gz files
    download_genomes(all_urls, prefix='fungi', n_files=1)
    run_child('bgzip fasta/fungi1.fasta')

def fetch_bovine():
    """Download cow genome and annotations."""
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
    chromosomes.extend(['chrX']) # chrY is missing 
    logging.info('Downloading bovine genome')
    local_file_name = os.path.join(target_dir, 'fasta', 'ref_Bos_taurus_GCF_002263795.1_ARS-UCD1.2.fasta')
    if os.path.exists(local_file_name):
        os.remove(local_file_name)
    for chrom in chromosomes:
        logging.debug('Downloading bovine chromosome %s', chrom)
        fasta_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.1_ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/%s.fna.gz' % chrom
        download_handle = ftp_down(fasta_url, local_file_name)
        download_handle.close()
        logging.debug('Downloaded bovine chromosome %s', chrom)
    fasta_url ='ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.1_ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_assembly_structure/non-nuclear/assembled_chromosomes/FASTA/chrMT.fna.gz'
    download_handle = ftp_down(fasta_url, local_file_name)
    download_handle.close()
    logging.debug('Downloaded bovine chromosome MT')
    fasta_url ='ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.1_ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_assembly_structure/Primary_Assembly/unplaced_scaffolds/FASTA/unplaced.scaf.fna.gz'
    download_handle = ftp_down(fasta_url, local_file_name)
    download_handle.close()
    logging.debug('Downloaded bovine chromosome unplaced')
    
    run_child('bgzip %s' % local_file_name)
    logging.info('Downloading gff annotation file')
    gff_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.1_ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_genomic.gff.gz'
    download_handle = ftp_down(gff_url)
    download_handle.close()


def main(args):
    """What the main does."""
    logging.info('now in fetch_data')
    if args.viral:
        fetch_viral(args.viral)
    if args.bact:
        fetch_bacterial()
    elif args.human:
        fetch_human()
    elif args.fungal:
        fetch_fungal()
    elif args.bovine:
        fetch_bovine()
