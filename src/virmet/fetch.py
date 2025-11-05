#!/usr/bin/env python3
"""Define all the functions that download reference sequences and other files."""

import logging
import os
import re

from virmet.common import (
    n_proc,
    ftp_down,
    get_accs,
    random_reduction,
    run_child,
    viral_query,
    find_file,
)


def fetch_viral(DB_DIR, viral_mode, n_proc, compression=True):
    """Download nucleotide or protein database."""
    # define the search nuccore/protein
    if viral_mode == "n":
        logging.info("downloading viral nuccore sequences")
        target_dir = os.path.join(DB_DIR, "viral_nuccore")
        ncbi_acc = viral_query(DB_DIR, "n")
    elif viral_mode == "p":
        logging.info("downloading viral protein sequences")
        target_dir = os.path.join(DB_DIR, "viral_protein")
        ncbi_acc = viral_query(DB_DIR, "p")
    else:
        raise ValueError(f"Invalid viral mode: {viral_mode}")

    # run the download
    logging.info("Database real path: %s" % os.path.realpath(target_dir))

    ncbi_file = os.path.join(target_dir, "ncbi_batch.txt")
    viral_database = os.path.join(target_dir, "viral_database.fasta")
    viral_seqs_info = os.path.join(target_dir, "viral_seqs_info.tsv")
    out_file = os.path.join(target_dir, "ncbi_dataset")
    for i in range(0, len(ncbi_acc), 10000):
        n_missing = min(len(ncbi_acc) - i, 10000)
        newfile = open(ncbi_file, "w")
        for acc_numb in ncbi_acc[i : i + n_missing]:
            newfile.write(re.sub(r"[']", "", str(acc_numb)) + "\n")
        newfile.close()
        cml_download = (
            "datasets download virus genome accession --inputfile %s --api-key $NCBI_API_KEY --filename %s.zip "
            % (ncbi_file, out_file)
        )
        run_child(cml_download)
        cml_extract = (
            "unzip -o %s.zip -d %s/; rm %s.zip; cat %s/data/genomic.fna >> %s ; \
            dataformat tsv virus-genome --inputfile %s/data/data_report.jsonl --fields accession,virus-tax-id,length,virus-name --elide-header >> %s; \
            rm -r %s"
            % (
                out_file,
                target_dir,
                out_file,
                out_file,
                viral_database,
                out_file,
                viral_seqs_info,
                out_file,
            )
        )
        run_child(cml_extract)
    os.remove(ncbi_file)

    # Report information about the download
    with open(viral_seqs_info) as f:
        total_acc_output = sum(1 for _ in f)
    logging.info("downloaded viral database info in %s" % target_dir)
    logging.info(
        "downloaded %d from a total of %d accession numbers"
        % (total_acc_output, len(ncbi_acc))
    )

    # viral_seqs_info.tsv contains Accn TaxId
    cml = "cut -f 1,2 %s > %s/viral_accn_taxid.dmp" % (
        viral_seqs_info,
        target_dir,
    )
    run_child(cml)
    accs_1 = set(get_accs(viral_database))
    vir_acc = open("%s/viral_accn_taxid.dmp" % target_dir)
    accs_2 = set([line.split()[0] for line in vir_acc])
    vir_acc.close()
    if accs_1 == accs_2:
        logging.info("taxonomy and fasta sequences match")
    else:
        logging.info("taxonomy and fasta sequences do not match")
        logging.info("fasta sequences: %s" % len(accs_1))
        logging.info("taxonomy information: %s" % len(accs_2))

    rmdup_cmd = f"seqkit rmdup {viral_database} --threads {n_proc} -i -o {target_dir}/viral_database_rmdup.fasta -D {target_dir}/duplicated_names.txt"
    run_child(rmdup_cmd)
    os.rename(viral_database, "%s/viral_database_original.fasta" % target_dir)
    os.rename("%s/viral_database_rmdup.fasta" % target_dir, viral_database)

    if compression:
        logging.info("Compress the database\n")
        random_reduction(viral_mode, DB_DIR)

    logging.info("downloading taxonomy databases")
    ftp_down(
        "https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz",
        os.path.join(DB_DIR, "taxdb.tar.gz"),
    )
    run_child("tar xvfz %s/taxdb.tar.gz -C %s --overwrite" % (DB_DIR, DB_DIR))
    os.remove("%s/taxdb.tar.gz" % DB_DIR)
    ftp_down(
        "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
        os.path.join(DB_DIR, "taxdump.tar.gz"),
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
        "%s/readme.txt" % DB_DIR,
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
        logging.debug(
            f"Could not find files {DB_DIR}/names.dmp, {DB_DIR}/nodes.dmp."
        )


def fetch_bact_fungal(DB_DIR, n_proc):
    """Download the bacterial and fungal databases from Kraken2."""
    target_dir = os.path.join(DB_DIR, "bact_fungi")
    logging.info("Database real path: %s" % os.path.realpath(target_dir))
    os.makedirs(target_dir, exist_ok=True)

    # download taxonomy first
    logging.info("Downloading taxonomy files")
    run_child(f"k2 download-taxonomy --db {target_dir}")
    # Download fungal and bacterial databases
    logging.info("Downloading fungal database")
    run_child(
        f"k2 download-library --library fungi --threads {n_proc} --db {target_dir} --no-masking"
    )
    logging.info("Downloading bacterial database")
    run_child(
        f"k2 download-library --library bacteria --threads {n_proc} --db {target_dir} --no-masking"
    )


def fetch_human(DB_DIR, n_proc):
    """Download human genome and annotations."""
    target_dir = os.path.join(DB_DIR, "human")
    logging.info("Database real path: %s" % os.path.realpath(target_dir))
    os.makedirs(target_dir, exist_ok=True)
    out_dir = os.path.join(target_dir, "fasta")
    os.makedirs(out_dir, exist_ok=True)
    base_name = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/"
    fasta_url = base_name + "GRCh38.primary_assembly.genome.fa.gz"
    version_name = find_file(
        base_name, r"gencode\.v[\d]+\.primary_assembly\.annotation\.gtf\.gz"
    )
    gtf_url = base_name + version_name
    logging.info("Downloading human annotation")
    ftp_down(gtf_url, os.path.join(out_dir, version_name))
    logging.info("Downloading human genome and bgzip compressing")
    fasta_path = os.path.join(out_dir, "GRCh38.fasta")
    if os.path.exists(fasta_path):
        os.remove(fasta_path)
    ftp_down(fasta_url, fasta_path)
    run_child(f"bgzip -@ {n_proc} -f {fasta_path}")


def fetch_bovine(DB_DIR, n_proc):
    """Download cow genome and annotations."""
    target_dir = os.path.join(DB_DIR, "bovine")
    logging.info("Database real path: %s" % os.path.realpath(target_dir))
    os.makedirs(target_dir, exist_ok=True)
    os.makedirs(os.path.join(target_dir, "fasta"), exist_ok=True)
    chromosomes = [f"chr{chrom}" for chrom in range(1, 30)]
    chromosomes.extend(["chrX"])  # chrY is missing
    logging.info("Downloading bovine genome")

    base_name = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/"
    version_name = find_file(base_name, r'href="([^"/]+/)"')
    new_base = (
        base_name
        + version_name
        + version_name.split("/")[0]
        + "_assembly_structure/"
    )
    local_file_name = os.path.join(target_dir, "fasta", "ref_Bos_taurus.fasta")
    if os.path.exists(local_file_name):
        os.remove(local_file_name)
    for chrom in chromosomes:
        logging.debug("Downloading bovine chromosome %s" % chrom)
        fasta_url = (
            new_base
            + f"Primary_Assembly/assembled_chromosomes/FASTA/{chrom}.fna.gz"
        )
        ftp_down(fasta_url, local_file_name)
        logging.debug("Downloaded bovine chromosome %s" % chrom)

    fasta_url = (
        new_base + "non-nuclear/assembled_chromosomes/FASTA/chrMT.fna.gz"
    )
    ftp_down(fasta_url, local_file_name)
    logging.debug("Downloaded bovine chromosome MT")

    fasta_url = (
        new_base
        + "Primary_Assembly/unplaced_scaffolds/FASTA/unplaced.scaf.fna.gz"
    )
    ftp_down(fasta_url, local_file_name)
    logging.debug("Downloaded bovine chromosome unplaced")

    run_child(f"bgzip -@ {n_proc} -f {local_file_name}")
    logging.info("Downloading gff annotation file")
    gff_url = (
        base_name
        + version_name
        + version_name.split("/")[0]
        + "_genomic.gff.gz"
    )
    ftp_down(
        gff_url,
        os.path.join(
            target_dir, "fasta", version_name.split("/")[0] + "_genomic.gff.gz"
        ),
    )


def main(args):
    """What the main does."""
    DB_DIR = os.path.expandvars(args.dbdir)
    logging.info("now in fetch_data")
    if args.viral:
        # print(args.no_db_compression)
        fetch_viral(
            DB_DIR, args.viral, n_proc, compression=not args.no_db_compression
        )
    if args.bact_fungal:
        fetch_bact_fungal(DB_DIR, n_proc)
    elif args.human:
        fetch_human(DB_DIR, n_proc)
    elif args.bovine:
        fetch_bovine(DB_DIR, n_proc)
