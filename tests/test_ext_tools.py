#!/usr/bin/env python3
import os
import sys
import glob
import unittest
import tempfile
import subprocess

virmet_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.sys.path.insert(1, virmet_dir)
mod = __import__('virmet')
sys.modules["virmet"] = mod

from virmet.common import run_child


class TestToolsCallable(unittest.TestCase):
    '''calls help of external tools just to check that they are installed
    '''

    def setUp(self):
        self.genome_file = os.path.join(tempfile.gettempdir(), 'HIV.fasta')

    def test_edirect(self):
        run_child('efetch -db nuccore -id K03455 -format fasta > %s' % self.genome_file)
        self.assertTrue(os.path.isfile(self.genome_file))
        os.remove(self.genome_file)

    def test_bwa_index(self):
        run_child('efetch -db nuccore -id K03455 -format fasta > %s' % self.genome_file)
        run_child('bwa index %s &> /dev/null' % self.genome_file)
        self.assertTrue(os.path.join(tempfile.gettempdir(), 'HIV.bwt'))
        os.remove(self.genome_file)
        for f in glob.glob('%s/HIV.*' % tempfile.gettempdir()):
            os.remove(f)

    def test_blast(self):
        log_file = os.path.join(tempfile.gettempdir(), 'tmp.log')
        run_child('blastn -help > %s 2>&1' % log_file)
        with open(log_file) as f:
            l = sum(1 for _ in f)
        self.assertGreater(l, 6)
        os.remove(log_file)

    # def tearDown(self):
    #     os.remove(self.log_file)
    #     os.remove(self.genome_file)
