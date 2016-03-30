#!/usr/bin/env python3.4

'''Runs on all samples of a MiSeq run or on a single fastq file'''
import os
import sys
import pandas as pd
from virmet.common import run_child

covpl_exe = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'scripts', 'covplot.R')

def main(args):
    ''''''
    outdir = args.outdir
    organism = args.organism

    try:
        os.chdir(outdir)
    except FileNotFoundError:
        sys.exit('Where is the output directory? Check the path.')
    for f in ['orgs_list.tsv', 'unique.tsv']:
        assert os.path.exists(f), 'file %s not found' % f

    orgs_list = pd.read_csv('orgs_list.tsv', sep='\t', header=0)
    # criterion is "startswith"
    criterion = orgs_list['organism'].map(lambda x: x.startswith(organism))
    matching_orgs = orgs_list[criterion]
    # organism matching that given on command line with most reads is:
    best_species = str(matching_orgs.iloc[0].organism)

    # blast results
    unique = pd.read_csv('unique.tsv', sep='\t', header=0)
    matching_reads = unique[unique['sscinames'] == best_species]

    best_seqids = matching_reads.groupby('sseqid').size().order(ascending=False)
    # TODO: upgrade for NCBI outphase of GI
    gi, dsc, acc = str(best_seqids.index.tolist()[0]).split('|')[1:4]

    # download single genome, index, align viral_reads
    cml = '-db nuccore -query \"%s[Accession]\" | efetch -format fasta > single.fasta' % acc
    run_child('esearch', cml)
    run_child('bwa', 'index single.fasta')
    run_child('bwa', 'mem single.fasta viral_reads.fastq.gz 2> /dev/null | samtools view -u - | samtools sort -O bam -T tmp -o single_sorted.bam -')
    run_child('samtools', 'index single_sorted.bam')
    run_child('samtools', 'depth -q 0 -Q 0 single_sorted.bam > depth.txt')
    image_name = organism.replace(' ', '_') + '_coverage.pdf'
    run_child('Rscript', '%s depth.txt %s %s' % (covpl_exe, acc, image_name))
