"""Tests for wolfpack module."""

import os
import tempfile
import unittest

import pandas as pd

from virmet.common import n_proc, run_child
from virmet.wolfpack import get_nodes_names, hunter, viral_blast

virmet_dir = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)


class TestHunter(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.gettempdir()
        self.reads = os.path.join(
            virmet_dir, "VirMet", "data", "raw_reads.fastq.gz"
        )

    def test_hunter(self):
        s_dir = hunter(self.reads, self.tmpdir, n_proc)
        raw_reads = run_child("gunzip -c %s | wc -l" % self.reads)
        raw_reads = int(raw_reads.strip().split()[0]) / 4
        with open("%s/good.fastq" % s_dir) as f:
            good_n = sum(1 for line in f) / 4
        with open("%s/stats.tsv" % s_dir) as f:
            stats = dict(line.strip().split() for line in f)
        filtered_out = int(stats["low_entropy"]) + int(stats["low_quality"])
        self.assertEqual(good_n, int(stats["passing_filter"]))
        self.assertEqual(
            raw_reads, filtered_out + good_n + int(stats["trimmed_too_short"])
        )

    def tearDown(self):
        for f in [
            "%s/good.fastq" % self.tmpdir,
            "%s/stats.tsv" % self.tmpdir,
            "%s/fastp.json" % self.tmpdir,
            "%s/fastp.err" % self.tmpdir,
            "%s/fastp.html" % self.tmpdir,
        ]:
            os.remove(f)


class TestViralBlast(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.gettempdir()
        self.dbdir = os.path.join(virmet_dir, "VirMet", "data", "test_db")
        self.nodes, self.names = get_nodes_names(self.dbdir)
        self.reads = os.path.join(
            virmet_dir, "VirMet", "data", "hq_decont_reads.fastq"
        )

    def test_viral_blast(self):
        viral_blast(
            self.reads, 4, self.nodes, self.names, self.tmpdir, self.dbdir
        )

        orgs_file = os.path.join(os.path.split(self.reads)[0], "orgs_list.tsv")
        df_org_list = pd.read_csv(
            orgs_file, index_col="species", delimiter="\t"
        )

        MG212469_reads = df_org_list.loc[
            "Enterovirus coxsackiepol", "reads"
        ].sum()
        MG212469_reads_expected = 287
        self.assertEqual(MG212469_reads, MG212469_reads_expected)

    def tearDown(self):
        path_remove = os.path.split(self.reads)[0]
        for f in [
            "%s/blast_info.txt" % path_remove,
            "%s/hq_decont_reads.fasta" % path_remove,
            "%s/orgs_list.tsv" % path_remove,
            "%s/stats.tsv" % path_remove,
            "%s/unique.tsv" % path_remove,
        ]:
            os.remove(f)
