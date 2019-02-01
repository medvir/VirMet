#!/usr/bin/env python3
import os
import io
import sys
import unittest
import tempfile
import pandas as pd
from pandas.util.testing import assert_frame_equal 

virmet_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.sys.path.insert(1, virmet_dir)
mod = __import__('virmet')
sys.modules["virmet"] = mod

from virmet.wolfpack import hunter, viral_blast, get_nodes_names
from virmet.common import run_child, DB_DIR


class TestHunter(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.gettempdir()
        self.reads = os.path.join(virmet_dir, 'VirMet', 'data', 'raw_reads.fastq.gz')

    def test_hunter(self):
        os.chdir(self.tmpdir)
        s_dir = hunter(self.reads)
        os.chdir(s_dir)
        raw_reads = run_child('gunzip -c %s | wc -l' % self.reads)
        raw_reads = int(raw_reads.strip().split()[0]) / 4
        with open('good.fastq') as f:
            good_n = sum(1 for l in f) / 4
        with open('bad.fastq') as f:
            bad_n = sum(1 for l in f) / 4
        with open('stats.tsv') as f:
            stats = dict(l.strip().split() for l in f)
        filtered_out = int(stats['low_entropy']) + int(stats['low_quality'])
        self.assertEqual(bad_n, filtered_out)
        self.assertEqual(good_n, int(stats['passing_filter']))
        self.assertEqual(raw_reads, bad_n + good_n + int(stats['trimmed_too_short']))

    def tearDown(self):
        os.chdir(self.tmpdir)
        for f in ['good.fastq', 'bad.fastq', 'prinseq.log', 'prinseq.err', 'stats.tsv']:
            os.remove(f)

    
class TestViralBlast(unittest.TestCase):
    
    def setUp(self):
        self.tmpdir =  tempfile.gettempdir()
        self.nodes, self.names = get_nodes_names(DB_DIR)
        self.reads = os.path.join(virmet_dir, 'VirMet', 'data', 'hq_decont_reads.fastq')
        			
    def test_viral_blast(self):
        os.chdir(self.tmpdir)
        viral_blast(self.reads, 4, self.nodes, self.names)
        
        self.orgs_file = os.path.join(self.tmpdir, 'orgs_list.tsv')
        df_org_list = pd.read_csv(self.orgs_file, index_col='species', delimiter="\t")
        
        MG212469_reads = df_org_list.loc['Enterovirus C','reads'].sum()
        MG212469_reads_expected = 1122
        self.assertEqual(MG212469_reads, MG212469_reads_expected)
        #assert_frame_equal(self.fixture, df_row)

        
