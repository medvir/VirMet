#!/usr/bin/env python

import glob
import subprocess
import sys
from collections import Counter


def run_child(cmd, exe="/bin/bash"):
    """use subrocess.check_output to run an external program with arguments"""
    try:
        output = subprocess.check_output(
            cmd,
            universal_newlines=True,
            shell=True,  # nosec B602: Required for shell piping.
            executable=exe,
            stderr=subprocess.STDOUT,
        )
    except subprocess.CalledProcessError as ee:
        print(
            "Execution of %s failed with returncode %d: %s"
            % (cmd, ee.returncode, ee.output)
        )
        output = None
    return output


def extract_reads(db_type):
    all_reads = Counter()
    for f in glob.glob("%s/good_*%s?.cram" % (sample_dir, db_type)):
        cml = "samtools view %s | cut -f 3" % f
        all_reads.update(run_child(cml).strip().split("\n"))

    print("In total", sum(all_reads.values()), db_type, "reads\n")

    print("|   reads   |    accn   |          organism         |")
    print("|----------:|:---------:|:--------------------------|")
    for accn, reads in all_reads.most_common(3):
        # cml = 'efetch -db nuccore -id %s -format docsum | xtract -pattern DocumentSummary -element Organism' % accn
        org = accn  #  run_child(cml).strip()
        print("| %s | %s | %s |" % (reads, accn, org))


sample_dir = sys.argv[1]

print("## Analysis on sample %s\n" % sample_dir)
# here the bacterial reads
print("\n### Bacterial reads\n")
extract_reads("bact")

# fungal reads
print("\n### Fungal reads\n")
extract_reads("fungi")
