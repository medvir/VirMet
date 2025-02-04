#!/usr/bin/env python3
"""Define all the functions that download reference sequences and other files."""

import logging
import os
import re

from virmet.common import (
    DB_DIR_UPDATE,
    N_FILES_BACT,
    n_proc,
    bact_fung_query,
    download_genomes,
    ftp_down,
    get_accs,
    random_reduction,
    run_child,
    viral_query,
)

DB_DIR = DB_DIR_UPDATE


def fetch_viral(viral_mode, compression=True):
    """Download nucleotide or protein database."""
    # define the search nuccore/protein
    if viral_mode == "n":
        logging.info("downloading viral nuccore sequences")
        target_dir = os.path.join(DB_DIR, "viral_nuccore")
        ncbi_acc = viral_query("n")
    elif viral_mode == "p":
        logging.info("downloaded viral protein sequences")
        target_dir = os.path.join(DB_DIR, "viral_protein")
        ncbi_acc = viral_query("p")
    else:
        raise ValueError(f"Invalid viral mode: {viral_mode}")
    
    # run the download
    logging.info("Database real path: ", os.path.realpath(target_dir))
    os.chdir(target_dir)

    ncbi_file = os.path.join(target_dir, "ncbi_batch.txt")
    viral_database = os.path.join(target_dir, "viral_database.fasta")
    viral_seqs_info = os.path.join(target_dir, "viral_seqs_info.tsv")
    out_file = os.path.join(target_dir, "ncbi_dataset")
    for i in range(0,len(ncbi_acc), 10000):
        newfile = open(ncbi_file,'w')
        for acc_numb in ncbi_acc[i:i+10000]:
            newfile.write(acc_numb+"\n")
        newfile.close()
        cml_download = "datasets download virus genome accession --inputfile %s --api-key $NCBI_API_KEY" % ncbi_file
        run_child(cml_download)
        cml_extract = (
            "unzip %s.zip; cat %s/data/genomic.fna >> %s ; \
            dataformat tsv virus-genome --inputfile %s/data/data_report.jsonl --fields accession,virus-tax-id,virus-name >> %s" 
            % (out_file, out_file, viral_database, out_file, viral_seqs_info)
            )
        run_child(cml_extract)
    os.remove(ncbi_file)
    os.remove(out_file)

    # Report information about the download
    with open(viral_seqs_info) as f:
        total_acc_output = sum(1 for _ in f)
    logging.info("downloaded viral database info in %s", target_dir)
    logging.info("downloaded %s from a total of %s accession numbers", total_acc_output, len(ncbi_acc))
   
    # viral_seqs_info.tsv contains Accn TaxId
    cml = "cut -f 1,2 viral_seqs_info.tsv > viral_accn_taxid.dmp"
    run_child(cml)
    accs_1 = set(get_accs("viral_database.fasta"))
    accs_2 = set([line.split()[0] for line in open("viral_accn_taxid.dmp")])
    assert accs_1 == accs_2, accs_1 ^ accs_2
    logging.info("taxonomy and fasta sequences match")

    rmdup_cmd = "cat viral_database.fasta | seqkit rmdup -i -o viral_database_rmdup.fasta -D duplicated_names.txt"
    run_child(rmdup_cmd)
    os.rename("viral_database.fasta", "viral_database_original.fasta")
    os.rename("viral_database_rmdup.fasta", "viral_database.fasta")

    if compression:
        logging.info("Compress the database\n")
        random_reduction(viral_mode)

    os.chdir(DB_DIR)
    logging.info("downloading taxonomy databases")
    download_handle = ftp_down(
        "https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz"
    )
    download_handle.close()
    run_child("tar xvfz taxdb.tar.gz")
    os.remove("taxdb.tar.gz")
    download_handle = ftp_down(
        "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
    )
    download_handle.close()
    run_child("tar xvfz taxdump.tar.gz")
    for ftd in [
        "taxdump.tar.gz",
        "merged.dmp",
        "gencode.dmp",
        "division.dmp",
        "delnodes.dmp",
        "citations.dmp",
    ]:
        try:
            os.remove(ftd)
        except OSError:
            logging.warning("Could not find file ", ftd)
    try:
        run_child(f"bgzip -@ {n_proc} names.dmp")
        run_child(f"bgzip -@ {n_proc} nodes.dmp")
    except Exception:
        logging.debug("Could not find files names.dmp, nodes.dmp.")


def fetch_bacterial():
    """Download the three bacterial sequence databases."""
    target_dir = os.path.join(DB_DIR, "bacteria")
    logging.info("Database real path: %s" % os.path.realpath(target_dir))
    os.makedirs(target_dir, exist_ok=True)
    os.chdir(target_dir)

    # first download summary file with all ftp paths and return urls
    all_urls = bact_fung_query(query_type="bacteria", target_folder=target_dir)
    # print(all_urls)
    logging.info("%d bacterial genomes were found", len(all_urls))
    # then download genomic_fna.gz files
    mid = len(all_urls) // 2
    download_genomes(all_urls[:mid], prefix="bact", n_files=N_FILES_BACT)
    print("second half starts")
    download_genomes(all_urls[mid:], prefix="bact", n_files=N_FILES_BACT)
    for j in range(1, N_FILES_BACT + 1):
        run_child(f"bgzip -@ {n_proc} fasta/bact{j}.fasta")


def fetch_human():
    """Download human genome and annotations."""
    target_dir = os.path.join(DB_DIR, "human")
    logging.info("Database real path: %s", os.path.realpath(target_dir))
    os.makedirs(target_dir, exist_ok=True)
    os.chdir(target_dir)
    os.makedirs("fasta", exist_ok=True)
    os.chdir("fasta")
    fasta_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/GRCh38.primary_assembly.genome.fa.gz"
    gtf_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v47.primary_assembly.annotation.gtf.gz"
    logging.info("Downloading human annotation")
    download_handle = ftp_down(gtf_url)
    download_handle.close()
    logging.info("Downloading human genome and bgzip compressing")
    if os.path.exists("GRCh38.fasta"):
        os.remove("GRCh38.fasta")
    download_handle = ftp_down(fasta_url, "GRCh38.fasta")
    download_handle.close()
    run_child(f"bgzip -@ {n_proc} GRCh38.fasta")


def fetch_fungal():
    """Download fungal sequences."""
    target_dir = os.path.join(DB_DIR, "fungi")
    logging.info("Database real path: %s", os.path.realpath(target_dir))
    os.makedirs(target_dir, exist_ok=True)
    os.chdir(target_dir)

    # first download summary file with all ftp paths and return urls
    all_urls = bact_fung_query(query_type="fungi", target_folder=target_dir)
    logging.info("%d fungal genomes were found", len(all_urls))
    # then download genomic_fna.gz files
    download_genomes(all_urls, prefix="fungi", n_files=1)
    run_child(f"bgzip -@ {n_proc} fasta/fungi1.fasta")


def fetch_bovine():
    """Download cow genome and annotations."""
    target_dir = os.path.join(DB_DIR, "bovine")
    logging.info("Database real path: %s", os.path.realpath(target_dir))
    os.makedirs(target_dir, exist_ok=True)
    os.chdir(target_dir)
    os.makedirs("fasta", exist_ok=True)
    os.chdir("fasta")
    chromosomes = [f"chr{chrom}" for chrom in range(1, 30)]
    chromosomes.extend(["chrX"])  # chrY is missing
    logging.info("Downloading bovine genome")
    local_file_name = os.path.join(
        target_dir, "fasta", "ref_Bos_taurus_GCF_002263795.3_ARS-UCD2.0.fasta"
    )
    if os.path.exists(local_file_name):
        os.remove(local_file_name)
    for chrom in chromosomes:
        logging.debug("Downloading bovine chromosome %s", chrom)
        fasta_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.3_ARS-UCD2.0/GCF_002263795.3_ARS-UCD2.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/{chrom}.fna.gz"
        download_handle = ftp_down(fasta_url, local_file_name)
        download_handle.close()
        logging.debug("Downloaded bovine chromosome ", chrom)
    fasta_url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.3_ARS-UCD2.0/GCF_002263795.3_ARS-UCD2.0_assembly_structure/non-nuclear/assembled_chromosomes/FASTA/chrMT.fna.gz"
    download_handle = ftp_down(fasta_url, local_file_name)
    download_handle.close()
    logging.debug("Downloaded bovine chromosome MT")
    fasta_url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.3_ARS-UCD2.0/GCF_002263795.3_ARS-UCD2.0_assembly_structure/Primary_Assembly/unplaced_scaffolds/FASTA/unplaced.scaf.fna.gz"
    download_handle = ftp_down(fasta_url, local_file_name)
    download_handle.close()
    logging.debug("Downloaded bovine chromosome unplaced")

    run_child(f"bgzip -@ {n_proc} {local_file_name}")
    logging.info("Downloading gff annotation file")
    gff_url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.3_ARS-UCD2.0/GCF_002263795.3_ARS-UCD2.0_genomic.gff.gz"
    download_handle = ftp_down(gff_url)
    download_handle.close()


def main(args):
    """What the main does."""
    logging.info("now in fetch_data")
    if args.viral:
        # print(args.no_db_compression)
        fetch_viral(args.viral, compression=not args.no_db_compression)
    if args.bact:
        fetch_bacterial()
    elif args.human:
        fetch_human()
    elif args.fungal:
        fetch_fungal()
    elif args.bovine:
        fetch_bovine()
