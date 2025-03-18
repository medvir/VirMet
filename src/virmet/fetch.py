#!/usr/bin/env python3
"""Define all the functions that download reference sequences and other files."""

import logging
import os
import re

from virmet.common import (
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

def fetch_viral(DB_DIR, viral_mode, compression=True):
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
    logging.info("Database real path: %s" % os.path.realpath(target_dir))

    ncbi_file = os.path.join(target_dir, "ncbi_batch.txt")
    viral_database = os.path.join(target_dir, "viral_database.fasta")
    viral_seqs_info = os.path.join(target_dir, "viral_seqs_info.tsv")
    out_file = os.path.join(target_dir, "ncbi_dataset")
    for i in range(0,len(ncbi_acc), 10000):
        n_missing = min(len(ncbi_acc)-i, 10000)
        newfile = open(ncbi_file,'w')
        for acc_numb in ncbi_acc[i:i+n_missing]:
            newfile.write(re.sub(r"[']", "", str(acc_numb))+"\n")
        newfile.close()
        cml_download = "datasets download virus genome accession --inputfile %s --api-key $NCBI_API_KEY --filename %s.zip " % (ncbi_file, out_file)
        run_child(cml_download)
        cml_extract = (
            "unzip -o %s.zip -d %s/; rm %s.zip; cat %s/data/genomic.fna >> %s ; \
            dataformat tsv virus-genome --inputfile %s/data/data_report.jsonl --fields accession,virus-tax-id,length,virus-name --elide-header >> %s; \
            rm -r %s" 
            % (out_file, target_dir, out_file, out_file, viral_database, out_file, viral_seqs_info, out_file)
            )
        run_child(cml_extract)
    os.remove(ncbi_file)

    # Report information about the download
    with open(viral_seqs_info) as f:
        total_acc_output = sum(1 for _ in f)
    logging.info("downloaded viral database info in %s" % target_dir)
    logging.info("downloaded %d from a total of %d accession numbers" % (total_acc_output, len(ncbi_acc)))
   
    # viral_seqs_info.tsv contains Accn TaxId
    cml = "cut -f 1,2 %s > %s/viral_accn_taxid.dmp" % (viral_seqs_info, target_dir)
    run_child(cml)
    accs_1 = set(get_accs(viral_database))
    vir_acc = open("%s/viral_accn_taxid.dmp" % target_dir)
    accs_2 = set([line.split()[0] for line in vir_acc])
    vir_acc.close()
    assert accs_1 == accs_2, accs_1 ^ accs_2
    logging.info("taxonomy and fasta sequences match")

    rmdup_cmd = f'seqkit rmdup {viral_database} --threads {n_proc} -i -o {target_dir}/viral_database_rmdup.fasta -D {target_dir}/duplicated_names.txt'
    run_child(rmdup_cmd)
    os.rename(viral_database, "%s/viral_database_original.fasta" % target_dir)
    os.rename("%s/viral_database_rmdup.fasta" % target_dir, viral_database)

    if compression:
        logging.info("Compress the database\n")
        random_reduction(viral_mode, DB_DIR)

    logging.info("downloading taxonomy databases")
    ftp_down(
        "https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz",
        os.path.join(DB_DIR, "taxdb.tar.gz")
    )
    run_child("tar xvfz %s/taxdb.tar.gz -C %s --overwrite" % (DB_DIR, DB_DIR))
    os.remove("%s/taxdb.tar.gz" % DB_DIR)
    ftp_down(
        "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
        os.path.join(DB_DIR, "taxdump.tar.gz")
    )
    run_child("tar xvfz %s/taxdump.tar.gz -C %s --overwrite" % (DB_DIR, DB_DIR))
    for ftd in [
        "%s/taxdump.tar.gz" % DB_DIR,
        "%s/merged.dmp" % DB_DIR,
        "%s/gencode.dmp" % DB_DIR,
        "%s/division.dmp" % DB_DIR,
        "%s/delnodes.dmp" % DB_DIR,
        "%s/citations.dmp" % DB_DIR,
        "%s/images.dmp" % DB_DIR,
        "%s/taxonomy4blast.sqlite3" % DB_DIR,
    ]:
        try:
            os.remove(ftd)
        except OSError:
            logging.warning("Could not find file %s" % ftd)
    try:
        run_child(f"bgzip -@ {n_proc} -f {DB_DIR}/names.dmp")
        run_child(f"bgzip -@ {n_proc} -f {DB_DIR}/nodes.dmp")
    except Exception:
        logging.debug(f"Could not find files {DB_DIR}/names.dmp, {DB_DIR}/nodes.dmp.")


def fetch_bacterial(DB_DIR):
    """Download the three bacterial sequence databases."""
    target_dir = os.path.join(DB_DIR, "bacteria")
    logging.info("Database real path: %s" % os.path.realpath(target_dir))
    os.makedirs(target_dir, exist_ok=True)

    # first download summary file with all ftp paths and return urls
    all_urls = bact_fung_query(DB_DIR, query_type="bacteria", target_folder=target_dir)
    logging.info("%d bacterial genomes were found" % len(all_urls))
    # then download genomic_fna.gz files
    mid = len(all_urls) // 2
    download_genomes(all_urls[:mid], prefix="bact", target_dir=target_dir, n_files=N_FILES_BACT)
    print("second half starts")
    download_genomes(all_urls[mid:], prefix="bact", target_dir=target_dir, n_files=N_FILES_BACT)
    for j in range(1, N_FILES_BACT + 1):
        run_child(f"bgzip -@ {n_proc} -f {target_dir}/fasta/bact{j}.fasta")


def fetch_human(DB_DIR):
    """Download human genome and annotations."""
    target_dir = os.path.join(DB_DIR, "human")
    logging.info("Database real path: %s" % os.path.realpath(target_dir))
    os.makedirs(target_dir, exist_ok=True)
    out_dir = os.path.join(target_dir, "fasta")
    os.makedirs(out_dir, exist_ok=True)
    fasta_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/GRCh38.primary_assembly.genome.fa.gz"
    gtf_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v47.primary_assembly.annotation.gtf.gz"
    logging.info("Downloading human annotation")
    ftp_down(gtf_url, os.path.join(out_dir, "gencode.v47.primary_assembly.annotation.gtf.gz"))
    logging.info("Downloading human genome and bgzip compressing")
    fasta_path = os.path.join(out_dir, "GRCh38.fasta")
    if os.path.exists(fasta_path):
        os.remove(fasta_path)
    ftp_down(fasta_url, fasta_path)
    run_child(f"bgzip -@ {n_proc} -f {fasta_path}")


def fetch_fungal(DB_DIR):
    """Download fungal sequences."""
    target_dir = os.path.join(DB_DIR, "fungi")
    logging.info("Database real path: %s" % os.path.realpath(target_dir))
    os.makedirs(target_dir, exist_ok=True)

    # first download summary file with all ftp paths and return urls
    all_urls = bact_fung_query(DB_DIR, query_type="fungi", target_folder=target_dir)
    logging.info("%d fungal genomes were found" % len(all_urls))
    # then download genomic_fna.gz files
    download_genomes(all_urls, prefix="fungi", target_dir=target_dir, n_files=1)
    run_child(f"bgzip -@ {n_proc} -f {target_dir}/fasta/fungi1.fasta")


def fetch_bovine(DB_DIR):
    """Download cow genome and annotations."""
    target_dir = os.path.join(DB_DIR, "bovine")
    logging.info("Database real path: %s" % os.path.realpath(target_dir))
    os.makedirs(target_dir, exist_ok=True)
    os.makedirs(os.path.join(target_dir, "fasta"), exist_ok=True)
    chromosomes = [f"chr{chrom}" for chrom in range(1, 30)]
    chromosomes.extend(["chrX"])  # chrY is missing
    logging.info("Downloading bovine genome")
    local_file_name = os.path.join(
        target_dir, "fasta", "ref_Bos_taurus_GCF_002263795.3_ARS-UCD2.0.fasta"
    )
    if os.path.exists(local_file_name):
        os.remove(local_file_name)
    for chrom in chromosomes:
        logging.debug("Downloading bovine chromosome %s" % chrom)
        fasta_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.3_ARS-UCD2.0/GCF_002263795.3_ARS-UCD2.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/{chrom}.fna.gz"
        ftp_down(fasta_url, local_file_name)
        logging.debug("Downloaded bovine chromosome %s" % chrom)
    fasta_url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.3_ARS-UCD2.0/GCF_002263795.3_ARS-UCD2.0_assembly_structure/non-nuclear/assembled_chromosomes/FASTA/chrMT.fna.gz"
    ftp_down(fasta_url, local_file_name)
    logging.debug("Downloaded bovine chromosome MT")
    fasta_url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.3_ARS-UCD2.0/GCF_002263795.3_ARS-UCD2.0_assembly_structure/Primary_Assembly/unplaced_scaffolds/FASTA/unplaced.scaf.fna.gz"
    ftp_down(fasta_url, local_file_name)
    logging.debug("Downloaded bovine chromosome unplaced")

    run_child(f"bgzip -@ {n_proc} -f {local_file_name}")
    logging.info("Downloading gff annotation file")
    gff_url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.3_ARS-UCD2.0/GCF_002263795.3_ARS-UCD2.0_genomic.gff.gz"
    ftp_down(gff_url, os.path.join(target_dir, "fasta", "GCF_002263795.3_ARS-UCD2.0_genomic.gff.gz"))


def main(args):
    """What the main does."""
    DB_DIR = os.path.expandvars(args.dbdir)
    logging.info("now in fetch_data")
    if args.viral:
        # print(args.no_db_compression)
        fetch_viral(DB_DIR, args.viral, compression=not args.no_db_compression)
    if args.bact:
        fetch_bacterial(DB_DIR)
    elif args.human:
        fetch_human(DB_DIR)
    elif args.fungal:
        fetch_fungal(DB_DIR)
    elif args.bovine:
        fetch_bovine(DB_DIR)
