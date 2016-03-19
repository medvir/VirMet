#!/usr/bin/env python3.4
import os
import sys
import unittest
import subprocess

virmet_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.sys.path.insert(1, virmet_dir)
mod = __import__('virmet')
sys.modules["virmet"] = mod

from virmet.common import run_child


class Testbwa(unittest.TestCase):

    def setUp(self):
        self.remote_1 = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.primary_assembly.annotation.gtf.gz'
        self.remote_2 = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/_README.TXT'

    def test_help(self):
        out = run_child('bwa', 'index &> /tmp/tmp.log')
        l = sum(1 for l in open('/tmp/tmp.log'))
        self.assertGreater(l, 6)
        out = run_child('bwa', 'mem &> /tmp/tmp.log')
        l = sum(1 for l in open('/tmp/tmp.log'))
        self.assertGreater(l, 20)

    def tearDown(self):
        os.remove('/tmp/tmp.log')
