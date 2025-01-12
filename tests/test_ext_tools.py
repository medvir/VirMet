"""Tests for ext_tools module."""

import glob
import os
import tempfile
import unittest

from virmet.common import run_child


class TestToolsCallable(unittest.TestCase):
    """calls help of external tools just to check that they are installed"""

    def setUp(self):
        self.tmpdir = tempfile.gettempdir()
        self.genome_file = os.path.join(self.tmpdir, "HIV.fasta")

    def test_edirect(self):
        run_child(
            "efetch -db nuccore -id K03455 -format fasta > %s"
            % self.genome_file
        )
        self.assertTrue(os.path.isfile(self.genome_file))
        os.remove(self.genome_file)

    def test_bwa_index(self):
        run_child(
            f"efetch -db nuccore -id K03455 -format fasta > {self.genome_file}"
        )
        run_child(f"bwa index {self.genome_file} &> /dev/null")
        self.assertTrue(os.path.join(self.tmpdir, "HIV.bwt"))
        os.remove(self.genome_file)
        for f in glob.glob(f"{self.tmpdir}/HIV.*"):
            os.remove(f)

    def test_blast(self):
        log_file = os.path.join(self.tmpdir, "tmp.log")
        run_child(f"blastn -help > {log_file} 2>&1")
        with open(log_file) as f:
            self.assertGreater(len(f.read()), 6)
        os.remove(log_file)
