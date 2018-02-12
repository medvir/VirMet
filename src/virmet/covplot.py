#!/usr/bin/env python3

'''Runs on all samples of a MiSeq run or on a single fastq file'''
import os
import pandas as pd
from Bio import SeqIO
from virmet.common import run_child

covpl_exe = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'scripts', 'covplot.R')


def best_species(orgs_file, organism):
    orgs_list = pd.read_csv(orgs_file, sep='\t', header=0)
    # assert decreasing sorted
    diff = orgs_list['reads'] - orgs_list['reads'].shift(1)
    assert (diff > 0).sum() == 0, diff
    # criterion is "startswith"
    # criterion = orgs_list['sscinames'].map(lambda x: x.startswith(organism))
    criterion = orgs_list.iloc[:, 0].str.startswith(organism).fillna(False)
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
    best_seqids = matching_reads.groupby('sseqid').size().sort_values(ascending=False)

    dsc, acc = str(best_seqids.index.tolist()[0]).split('|')[:2]

    # copy single genome, index, align viral_reads
    os.chdir(outdir)
    organism = organism.replace(' ', '-')
    best_seq = [s for s in SeqIO.parse('/data/virmet_databases/viral_nuccore/viral_database.fasta', 'fasta') if acc in s.id]
    seq_len = len(best_seq[0])
    SeqIO.write(best_seq, 'single.fasta', 'fasta')
    run_child('bwa index single.fasta')
    bam_file = 'single_sorted_%s.bam' % organism
    run_child('bwa mem -t 2 single.fasta viral_reads.fastq.gz 2> /dev/null | samtools view -u - | samtools sort -O bam -T tmp -o %s -' % bam_file)
    run_child('samtools index %s' % bam_file)
    depth_file = 'depth_%s.txt' % organism
    run_child('samtools depth -a -q 0 -Q 0 %s > %s' % (bam_file, depth_file))
    image_name = organism.replace(' ', '_') + '_coverage.pdf'
    run_child('Rscript %s %s %s %s %s' % (covpl_exe, depth_file, acc, seq_len, image_name))

    return best_species

if __name__ == '__main__':
    import sys
    args = {'outdir': sys.argv[1], 'organism': sys.argv[2]}
    main(args)
