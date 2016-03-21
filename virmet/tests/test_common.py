#!/usr/bin/env python3.4
import os
import sys
import unittest
import tempfile

virmet_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.sys.path.insert(1, virmet_dir)
mod = __import__('virmet')
sys.modules["virmet"] = mod

from virmet.common import run_child, ftp_down


def parse_file_line(line):
    ftl = line.strip().split(':')[1][1:]
    if ftl.startswith('gzip compressed') and 'extra field' in ftl:
        return 'bgzipped'
    elif ftl.startswith('gzip compressed') and 'extra field' not in ftl:
        return 'gzipped'
    elif ftl.startswith('ASCII'):
        return 'ascii'


class TestFTPDown(unittest.TestCase):

    def setUp(self):
        self.remote_1 = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.primary_assembly.annotation.gtf.gz'
        self.remote_2 = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/_README.TXT'

    def test_nodecompress(self):
        out_file = os.path.join(tempfile.gettempdir(), 'gtf.txt.gz')
        ftp_down(self.remote_1, out_file)
        ftl = run_child('file', out_file)
        os.remove(out_file)
        ft = parse_file_line(ftl)
        self.assertEqual(ft, 'gzipped')

    def test_decompress(self):
        out_file = os.path.join(tempfile.gettempdir(), 'gtf.txt')
        ftp_down(self.remote_1, out_file)
        ftl = run_child('file', out_file)
        os.remove(out_file)
        ft = parse_file_line(ftl)
        self.assertEqual(ft, 'ascii')

    def test_append(self):
        out_file = os.path.join(tempfile.gettempdir(), 'README.TXT')
        ftp_down(self.remote_2, out_file)
        with open(out_file) as f:
            n_lines_1 = sum(1 for _ in f)
        ftp_down(self.remote_2, out_file)
        with open(out_file) as f:
            n_lines_2 = sum(1 for _ in f)
        self.assertEqual(n_lines_2, 2 * n_lines_1)
        os.remove(out_file)
