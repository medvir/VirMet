"""Tests for ext_tools module."""

import os
import tempfile
import unittest

from shutil import which
from virmet.common import run_child


class TestToolsCallable(unittest.TestCase):
    """check that external tools are installed"""

    def setUp(self):
        self.tmpdir = tempfile.gettempdir()

    def test_seqkit(self):
        self.assertIsNotNone(which("seqkit"))
    
    def test_samtools(self):
        self.assertIsNotNone(which("samtools"))

    def test_bwamem2(self):
        self.assertIsNotNone(which("bwa-mem2"))
    
    def test_fastp(self):
        self.assertIsNotNone(which("fastp"))
    
    def test_datasets(self):
        self.assertIsNotNone(which("datasets"))
        self.assertIsNotNone(which("dataformat"))

    def test_blast(self):
        self.assertIsNotNone(which("blastn"))
        self.assertIsNotNone(which("makeblastdb"))
        self.assertIsNotNone(which("blastdbcmd"))
    
    def test_bgzip(self):
        self.assertIsNotNone(which("bgzip"))
