#!/usr/bin/env python3
"""Define all the functions that download reference sequences and other files."""

import logging
import os
import pathlib
import shutil
import subprocess

from virmet.common import (
    DB_DIR_UPDATE,
    N_FILES_BACT,
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
        target_dir = pathlib.Path(DB_DIR) / "viral_nuccore"
        cml_search = viral_query("n")
    elif viral_mode == "p":
        logging.info("downloaded viral protein sequences")
        target_dir = pathlib.Path(DB_DIR) / "viral_protein"
        cml_search = viral_query("p")
    else:
        raise ValueError(f"Invalid viral mode: {viral_mode}")
    # run the search and download
    logging.info("Database real path: %s", os.path.realpath(target_dir))
    # os.chdir(target_dir)
    search_result = target_dir / "ncbi_search"

    with search_result.open("w") as stdin:
        run_child(
            [shutil.which("esearch")] + cml_search,
            stdout=stdin,
        )
    # cml_fetch_fasta = (
    #     "efetch -format fasta < ncbi_search > viral_database.fasta"
    # )
    viral_database_path = target_dir / "viral_databse.fasta"
    with (
        search_result.open("r") as stdin,
        viral_database_path.open("w") as stdout,
    ):
        run_child(
            [
                shutil.which("efetch"),
                "-format",
                "fasta",
            ],
            stdin=stdin,
            stdout=stdout,
        )

    with search_result.open("rb") as stdin:
        docsum = run_child(
            [
                shutil.which("efetch"),
                "-format",
                "docsum",
            ],
            stdin=stdin,
            stdout=subprocess.PIPE,
        )

    viral_seqs_info_path = target_dir / "viral_seqs_info.tsv"
    with viral_seqs_info_path.open("wb") as stdout:
        run_child(
            [
                shutil.which("xtract"),
                "-pattern",
                "DocumentSummary",
                "-element",
                "Caption",
                "TaxId",
                "Slen",
                "Organism",
                "Title",
                "AccessionVersion",
            ],
            stdin=docsum.stdout,
            stdout=stdout,
        )

    # cml_efetch_xtract = (
    #     "efetch -format docsum < ncbi_search | xtract "
    #     " -pattern DocumentSummary "
    #     "-element Caption TaxId Slen Organism Title AccessionVersion "
    #     "> viral_seqs_info.tsv"
    # )
    # run_child(cml_efetch_xtract)
    logging.info("downloaded viral seqs info in %s", target_dir)
    logging.info("saving viral taxonomy")
    # viral_seqs_info.tsv contains Accn TaxId
    viral_accn_taxid_file = target_dir / "viral_accn_taxid.dmp"
    with (
        viral_seqs_info_path.open("rb") as stdin,
        viral_accn_taxid_file.open("wb") as stdout,
    ):
        run_child(
            [
                shutil.which("cut"),
                "-f",
                "1,2",
            ],
            stdin=stdin,
            stdout=stdout,
        )
    # cml = "cut -f 1,2 viral_seqs_info.tsv > viral_accn_taxid.dmp"
    # run_child(cml)
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
        "ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz"
    )
    download_handle.close()
    run_child("tar xvfz taxdb.tar.gz")
    os.remove("taxdb.tar.gz")
    download_handle = ftp_down(
        "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
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
            logging.warning("Could not find file %s", ftd)
    try:
        run_child("bgzip names.dmp")
        run_child("bgzip nodes.dmp")
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
        run_child(f"bgzip fasta/bact{j}.fasta")


def fetch_human():
    """Download human genome and annotations."""
    target_dir = os.path.join(DB_DIR, "human")
    logging.info("Database real path: %s", os.path.realpath(target_dir))
    os.makedirs(target_dir, exist_ok=True)
    os.chdir(target_dir)
    os.makedirs("fasta", exist_ok=True)
    os.chdir("fasta")
    fasta_url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz"
    gtf_url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.primary_assembly.annotation.gtf.gz"
    logging.info("Downloading human annotation")
    download_handle = ftp_down(gtf_url)
    download_handle.close()
    logging.info("Downloading human genome and bgzip compressing")
    if os.path.exists("GRCh38.fasta"):
        os.remove("GRCh38.fasta")
    download_handle = ftp_down(fasta_url, "GRCh38.fasta")
    download_handle.close()
    run_child("bgzip GRCh38.fasta")


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
    run_child("bgzip fasta/fungi1.fasta")


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
        target_dir, "fasta", "ref_Bos_taurus_GCF_002263795.2_ARS-UCD1.3.fasta"
    )
    if os.path.exists(local_file_name):
        os.remove(local_file_name)
    for chrom in chromosomes:
        logging.debug("Downloading bovine chromosome %s", chrom)
        fasta_url = f"ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.2_ARS-UCD1.3/GCF_002263795.2_ARS-UCD1.3_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/{chrom}.fna.gz"
        download_handle = ftp_down(fasta_url, local_file_name)
        download_handle.close()
        logging.debug("Downloaded bovine chromosome %s", chrom)
    fasta_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.2_ARS-UCD1.3/GCF_002263795.2_ARS-UCD1.3_assembly_structure/non-nuclear/assembled_chromosomes/FASTA/chrMT.fna.gz"
    download_handle = ftp_down(fasta_url, local_file_name)
    download_handle.close()
    logging.debug("Downloaded bovine chromosome MT")
    fasta_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.2_ARS-UCD1.3/GCF_002263795.2_ARS-UCD1.3_assembly_structure/Primary_Assembly/unplaced_scaffolds/FASTA/unplaced.scaf.fna.gz"
    download_handle = ftp_down(fasta_url, local_file_name)
    download_handle.close()
    logging.debug("Downloaded bovine chromosome unplaced")

    run_child(f"bgzip {local_file_name}")
    logging.info("Downloading gff annotation file")
    gff_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.2_ARS-UCD1.3/GCF_002263795.2_ARS-UCD1.3_genomic.gff.gz"
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
