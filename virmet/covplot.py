#!/usr/bin/env python3.4

'''Runs on all samples of a MiSeq run or on a single fastq file'''
import os
import sys
import pandas as pd
from virmet.common import run_child

covpl_exe = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'scripts', 'covplot.R')


def best_species(orgs_file, organism):
    orgs_list = pd.read_csv(orgs_file, sep='\t', header=0)
    # assert decreasing sorted
    diff = orgs_list['reads'] - orgs_list['reads'].shift(1)
    assert (diff > 0).sum() == 0, diff
    # criterion is "startswith"
    criterion = orgs_list['organism'].map(lambda x: x.startswith(organism))
    matching_orgs = orgs_list[criterion]
    # organism matching that given on command line with most reads is the first
    # W.O. this assumes descending order of reads
    return str(matching_orgs.iloc[0].organism)


def main(args):
    ''''''
    outdir = args.outdir
    organism = args.organism

    assert os.path.isdir(outdir), 'Where is the output dir? Check the path.'

    org_file = os.path.join(outdir, 'orgs_list.tsv')
    best_spec = best_species(org_file, organism)

    # blast results
    blast_file = os.path.join(outdir, 'unique.tsv.gz')
    unique = pd.read_csv(blast_file, sep='\t', header=0, compression='gzip')
    matching_reads = unique[unique['sscinames'] == best_spec]
    best_seqids = matching_reads.groupby('sseqid').size().order(ascending=False)

    # TODO: upgrade for NCBI outphase of GI
    gi, dsc, acc = str(best_seqids.index.tolist()[0]).split('|')[1:4]

    os.chdir(outdir)
    # download single genome, index, align viral_reads
    cml = '-db nuccore -query \"%s[Accession]\" | efetch -format fasta > single.fasta' % acc
    run_child('esearch', cml)
    run_child('bwa', 'index single.fasta')
    run_child('bwa', 'mem single.fasta viral_reads.fastq.gz 2> /dev/null | samtools view -u - | samtools sort -O bam -T tmp -o single_sorted.bam -')
    run_child('samtools', 'index single_sorted.bam')
    run_child('samtools', 'depth -q 0 -Q 0 single_sorted.bam > depth.txt')
    image_name = organism.replace(' ', '_') + '_coverage.pdf'
    run_child('Rscript', '%s depth.txt %s %s' % (covpl_exe, acc, image_name))

    return best_species
