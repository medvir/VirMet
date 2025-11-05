#!/usr/bin/env python3
"""Update viral and bacterial database based on a new query to ncbi and a manually added list of GIs"""

import logging
import os
import re


from virmet.common import (
    n_proc,
    get_accs,
    run_child,
    viral_query,
    random_reduction,
)


def virupdate(
    DB_DIR, viral_type, picked=None, update_min_date=None, compression=True
):
    if viral_type == "n":
        db_type = "nuccore"
    elif viral_type == "p":
        db_type = "protein"
    else:
        raise ValueError(f'Invalid db_type value: "{db_type}".')
    target_dir = os.path.join(DB_DIR, f"viral_{db_type}")
    viral_database = os.path.join(target_dir, "viral_database_updated.fasta")
    viral_seqs_info = os.path.join(target_dir, "viral_seqs_info_updated.tsv")
    final_database = os.path.join(target_dir, "viral_database.fasta")
    final_viral_seqs = os.path.join(target_dir, "viral_seqs_info.tsv")

    # Set temporal file
    out_file = os.path.join(target_dir, "ncbi_dataset")

    if update_min_date:
        # this query downloads a new viral_seqs_info.tsv and parses the GI
        logging.info("Interrogating NCBI")
        ncbi_acc = viral_query(DB_DIR, viral_type, update_min_date)
        ncbi_file = os.path.join(target_dir, "ncbi_batch.txt")

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
        logging.info(
            "downloaded updated viral database info in %s" % target_dir
        )
        logging.info(
            "downloaded %d from a total of %d accession numbers"
            % (total_acc_output, len(ncbi_acc))
        )

    # sequences given manually by specifying file with ACCs
    if picked:
        cml_download = (
            "datasets download virus genome accession --inputfile %s --api-key $NCBI_API_KEY --filename %s.zip "
            % (picked, out_file)
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

    # Create new dmp file from viral_seqs_info.tsv, which contains Accn TaxId
    cml = "cut -f 1,2 %s > %s/viral_accn_taxid_updated.dmp" % (
        viral_seqs_info,
        target_dir,
    )
    run_child(cml)
    accs_1 = set(get_accs(viral_database))
    vir_acc = open("%s/viral_accn_taxid_updated.dmp" % target_dir)
    accs_2 = set([line.split()[0] for line in vir_acc])
    vir_acc.close()
    if accs_1 == accs_2:
        logging.info("taxonomy and fasta sequences match")
    else:
        logging.info("taxonomy and fasta sequences do not match")
        logging.info("fasta sequences: %s" % len(accs_1))
        logging.info("taxonomy information: %s" % len(accs_2))

    # Merge with current database and remove duplicates
    cml_merge = f"mv {final_database} {target_dir}/viral_database_outdated.fasta ; \
                cat {target_dir}/viral_database_outdated.fasta {viral_database} > {final_database} ; \
                mv {final_viral_seqs} {target_dir}/viral_seqs_info_outdated.tsv ; \
                cat {target_dir}/viral_seqs_info_outdated.tsv {viral_seqs_info} > {final_viral_seqs} ; \
                mv {target_dir}/viral_accn_taxid.dmp {target_dir}/viral_accn_taxid_outdated.dmp ; \
                cat {target_dir}/viral_accn_taxid_outdated.dmp {target_dir}/viral_accn_taxid_updated.dmp > {target_dir}/viral_accn_taxid.dmp; \
                rm -f {target_dir}/viral_database_original.fasta"
    run_child(cml_merge)

    rmdup_cmd = f"seqkit rmdup {final_database} --threads {n_proc} -i -o {target_dir}/viral_database_rmdup.fasta -D {target_dir}/duplicated_names.txt"
    run_child(rmdup_cmd)
    os.rename(final_database, "%s/viral_database_original.fasta" % target_dir)
    os.rename("%s/viral_database_rmdup.fasta" % target_dir, final_database)

    if compression:
        logging.info("Compress the database\n")
        random_reduction(viral_type, DB_DIR)


def main(args):
    DB_DIR = os.path.expandvars(args.dbdir)
    logging.info("now in update_db")
    logging.info("Database real path: %s" % os.path.realpath(DB_DIR))
    virupdate(
        DB_DIR,
        args.viral,
        args.picked,
        args.update_min_date,
        compression=not args.no_db_compression,
    )
