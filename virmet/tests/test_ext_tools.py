#!/usr/bin/env python3.4
import os
import sys
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
        self.log_file = os.path.join(tempfile.gettempdir(), 'tmp.log')

    def test_bwa(self):
        run_child('bwa', 'index > %s 2>&1' % self.log_file)
        with open(self.log_file) as f:
            l = sum(1 for _ in f)
        self.assertGreater(l, 6)

    def test_blast(self):
        run_child('blastn', '-H > %s 2>&1' % self.log_file)
        with open(self.log_file) as f:
            l = sum(1 for _ in f)
        self.assertGreater(l, 6)

    def tearDown(self):
        os.remove(self.log_file)
