#!/usr/bin/env python3.4
import os
import sys
import unittest
import subprocess

def run_child(exe_name, arg_string):
    '''use subrocess to run an external program with arguments'''
    import subprocess

    if not arg_string.startswith(' '):
        arg_string = ' ' + arg_string

    try:
        retcode = subprocess.call(exe_name + arg_string, shell=True)
        if retcode > 0:
            sys.exit("Child %s %s terminated by signal %d" %
                     (exe_name, arg_string, retcode))
    except OSError as ee:
        sys.exit("Execution of %s failed: %s" % (exe_name, ee))

    return retcode

class TestSTAR(unittest.TestCase):

    def setUp(self):
        pass

    def testbasicpass(self):
        self.assertEqual('abc', 'abc')

    def testmakeindex(self):
        # download HIV genome
        run_child('efetch', '-db nuccore -id K03455 -format fasta > /tmp/genome.fasta')
        self.assertTrue(os.path.isfile('/tmp/genome.fasta'))
        # index
        run_child('STAR', '--runMode genomeGenerate --genomeSAindexNbases 6 \
        --genomeDir /tmp --genomeFastaFiles /tmp/genome.fasta')
        # download reads
        #run_child('wget', 'https://github.com/ozagordi/VirMet/blob/master/data/read1.fastq -O /tmp/read1.fastq')
        # align reads not possible since STAR hangs on small genomes
        #run_child('STAR', '--genomeDir /tmp --readFilesIn /tmp/read1.fastq \
        #--outFileNamePrefix /tmp/')
        self.assertTrue(os.path.isfile('/tmp/SA'))
        self.assertTrue(os.path.isfile('/tmp/SAindex'))
        self.assertTrue(os.path.isfile('/tmp/Genome'))

if __name__ == '__main__':
    unittest.main()
