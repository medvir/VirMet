#!/usr/bin/env python3.4

'''Runs on all samples of a MiSeq run or on a single fastq file'''
import os
import sys
import glob
import logging
import pandas as pd
from virmet.common import run_child, single_process, DB_DIR

contaminant_db = ['/data/virmet_databases/human/bwa/humanGRCh38',
                  '/data/virmet_databases/bacteria/bwa/bact1',
                  '/data/virmet_databases/bacteria/bwa/bact2',
                  '/data/virmet_databases/bacteria/bwa/bact3',
                  '/data/virmet_databases/fungi/bwa/fungi1',
                  '/data/virmet_databases/bovine/bwa/bt_ref']
ref_map = {
    'humanGRCh38': '/data/virmet_databases/human/fasta/GRCh38.fasta.gz',
    'bact1': '/data/virmet_databases/bacteria/fasta/bact1.fasta.gz',
    'bact2': '/data/virmet_databases/bacteria/fasta/bact2.fasta.gz',
    'bact3': '/data/virmet_databases/bacteria/fasta/bact3.fasta.gz',
    'fungi1': '/data/virmet_databases/fungi/fasta/fungi1.fasta.gz',
    'bt_ref': '/data/virmet_databases/bovine/fasta/bt_ref_Bos_taurus_UMD_3.1.1.fasta.gz'
}

blast_cov_threshold = 75.
blast_ident_threshold = 75.


def hunter(fq_file):
    '''runs quality filter on a fastq file with seqtk and prinseq,
    simple parallelisation with xargs, returns output directory
    '''
    from virmet.common import prinseq_exe

    try:
        n_proc = min(os.cpu_count(), 16)
    except NotImplementedError:
        n_proc = 2

    logging.info('hunter will run on %s processors' % n_proc)
    if '_' in fq_file:
        s_dir = '_'.join(os.path.split(fq_file)[1].split('_')[:2])
        try:
            os.mkdir(s_dir)
        except FileExistsError:
            logging.info('entering %s already existing' % s_dir)
        os.chdir(s_dir)
        s_dir = os.getcwd()
    else:
        s_dir = os.getcwd()

    oh = open('stats.tsv', 'a')
    # count raw reads
    if fq_file.endswith('gz'):
        out1 = run_child('gunzip', '-c %s | wc -l' % fq_file)
    else:
        out1 = run_child('wc', '-l %s | cut -f 1 -d \" \"' % fq_file)
    n_reads = int(int(out1) / 4)
    oh.write('raw_reads\t%d\n' % n_reads)

    # trim and discard short reads, count
    logging.info('trimming with seqtk')
    cml = 'trimfq %s | seqtk seq -L 75 - > intermediate.fastq' % fq_file
    out1 = run_child('seqtk', cml)
    out1 = run_child('wc', '-l intermediate.fastq | cut -f 1 -d \" \"')
    n_reads = int(int(out1) / 4)
    oh.write('trimmed_long\t%d\n' % n_reads)

    # We want to split in n_proc processors, so each file has at most
    # (n_reads / n_proc) + 1 reads and 4 times as many lines
    # this fails if there are more cpus than reads!
    max_reads_per_file = int(n_reads / n_proc) + 1
    max_l = max_reads_per_file * 4
    # split and rename
    run_child('split', '-l %d intermediate.fastq splitted' % max_l)
    os.remove('intermediate.fastq')
    splitted = glob.glob('splitted*')
    n_splitted = len(splitted)
    for i, spf in enumerate(sorted(splitted)):
        os.rename(spf, 'splitted%0.2d.fastq' % i)  # W.O. max 100 files/cpus

    # filter with prinseq, parallelize with xargs
    logging.info('filtering with prinseq')
    cml = '-w 0 %d | xargs -P %d -I {} %s \
            -fastq splitted{}.fastq -lc_method entropy -lc_threshold 70 \
            -log prinseq{}.log -min_qual_mean 20 -ns_max_p 25 \
            -out_good ./good{} -out_bad ./bad{} > ./prinseq.err 2>&1' % (n_splitted - 1, n_splitted, prinseq_exe)
    run_child('seq', cml)

    logging.info('cleaning up')
    if len(glob.glob('good??.fastq')):
        run_child('cat', 'good??.fastq > good.fastq')
        run_child('rm', 'good??.fastq')

    if len(glob.glob('bad??.fastq')):
        run_child('cat', 'bad??.fastq > bad.fastq')
        run_child('rm', 'bad??.fastq')

    if len(glob.glob('prinseq??.log')):
        run_child('cat', 'prinseq??.log > prinseq.log')
        run_child('rm', 'prinseq??.log')
    run_child('rm', 'splitted*fastq')

    out1 = run_child('wc', '-l good.fastq | cut -f 1 -d \" \"')
    n_reads = int(int(out1) / 4)
    oh.write('high_quality\t%d\n' % n_reads)

    os.chdir(os.pardir)
    return s_dir


def victor(input_reads, contaminant):
    '''decontaminate reads by aligning against contaminants with bwa and removing
    reads with alignments
    '''
    import gzip
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    try:
        n_proc = min(os.cpu_count(), 16)
    except NotImplementedError:
        n_proc = 2

    # alignment with bwa
    rf_head = input_reads.split('.')[0]
    cont_name = os.path.split(contaminant)[1]
    sam_name = '%s_%s.sam' % (rf_head, cont_name)
    err_name = '%s_%s.err' % (rf_head, cont_name)
    cml = 'mem -t %d -R \'@RG\tID:foo\tSM:bar\tLB:library1\' -T 75 %s %s 2> \
    %s | samtools view -h -F 4 - > %s' % (n_proc, contaminant, input_reads, err_name, sam_name)
    run_child('bwa', cml)
    logging.info('running bwa %s %s on %d cores' % (cont_name, rf_head, n_proc))

    # reading sam file to remove reads with hits
    # test if an object is in set is way faster than in list
    mapped_reads = set(run_child('grep', '-v \"^@\" %s | cut -f 1' % sam_name).strip().split('\n'))
    # mapped_reads = set([l.split()[0] for l in open(sam_name)
    #                     if not l.startswith('@')])
    oh = open('stats.tsv', 'a')
    oh.write('matching_%s\t%d\n' % (cont_name, len(mapped_reads)))
    oh.close()
    clean_name = os.path.splitext(sam_name)[0] + '.fastq'

    output_handle = open(clean_name, 'w')
    logging.info('Cleaning reads in %s with alignments in %s' %
                 (input_reads, sam_name))
    logging.info('Writing to %s' % clean_name)
    if input_reads.endswith('.gz'):
        cont_handle = gzip.open(input_reads)
    else:
        cont_handle = open(input_reads)
    c = 0
    # Using FastqGeneralIterator allows fast performance
    for title, seq, qual in FastqGeneralIterator(cont_handle):
        if title.split()[0] not in mapped_reads:
            c += 1
            output_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
            if c % 100000 == 0:
                logging.info('written %d clean reads' % c)
    logging.info('written %d clean reads' % c)
    output_handle.close()

    return clean_name


def viral_blast(file_in, n_proc):
    '''runs blast against viral database, parallelise with xargs
    '''

    oh = open('stats.tsv', 'a')
    os.rename(file_in, 'hq_decont_reads.fastq')
    fasta_file = 'hq_decont_reads.fasta'
    run_child('seqtk', 'seq -A hq_decont_reads.fastq > %s' % fasta_file)
    tot_seqs = int(run_child('grep', '-c \"^>\" %s' % fasta_file).strip())
    oh.write('blasted_reads\t%d\n' % tot_seqs)
    max_n = (tot_seqs / n_proc) + 1

    # We want to split in n_proc processors, so each file has at most
    # (tot_seqs / n_proc) + 1 reads
    cml = "-v \"MAX_N=%d\" \'BEGIN {n_seq=0;} /^>/ \
    {if(n_seq %% %d == 0){file=sprintf(\"splitted_clean_%%d.fasta\", n_seq/%d);} \
    print >> file; n_seq++; next;} { print >> file; }' %s" % (max_n, max_n, max_n, fasta_file)
    run_child('awk', cml)

    # blast needs access to taxdb files to retrieve organism name
    os.environ['BLASTDB'] = DB_DIR

    xargs_thread = 0  # means on all available cores, caution
    cml = '0 %s | xargs -P %d -I {} blastn -task megablast \
           -query splitted_clean_{}.fasta -db %s \
           -out tmp_{}.tsv \
           -outfmt \'6 qseqid sseqid sscinames stitle pident qcovs score length mismatch gapopen qstart qend sstart send staxids\'' \
        % (n_proc - 1, xargs_thread, os.path.join('/data/virmet_databases', 'viral_nuccore/viral_db'))
    logging.info('running blast now')
    run_child('seq', cml)

    logging.info('parsing best HSP for each query sequence')
    qseqid = ''

    bh = open('unique.tsv', 'w')
    bh.write('qseqid\tsseqid\tsscinames\tstitle\tpident\tqcovs\tscore\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tstaxids\n')
    for i in range(n_proc):
        tmpf = 'tmp_%d.tsv' % i
        with open(tmpf) as f:
            for line in f:
                if line.split('\t')[0] != qseqid:
                    bh.write(line)
                    qseqid = line.split('\t')[0]
        os.remove(tmpf)
        os.remove('splitted_clean_%d.fasta' % i)
    bh.close()

    logging.info('filtering and grouping by scientific name')
    hits = pd.read_csv('unique.tsv', index_col='qseqid',  # delim_whitespace=True)
                       delimiter="\t")
    logging.info('found %d hits' % hits.shape[0])

    oh.write('blast_hits\t%s\n' % hits.shape[0])
    # select according to identity and coverage
    good_hits = hits[(hits.pident > blast_ident_threshold) & (hits.qcovs > blast_cov_threshold)]
    # define a column with genbank id only
    # good_hits['sid'] = [int(s.split('_')[1]) for s in good_hits.sseqid]
    # good_hits = good_hits.join(taxonomy, on='sid')
    matched_reads = good_hits.shape[0]
    logging.info('%d hits passing coverage and identity filter' % matched_reads)
    oh.write('viral_good_reads\t%s\n' % matched_reads)
    org_count = good_hits.groupby('sscinames').size()
    org_count.order(ascending=False).to_csv('orgs_list.csv', header=True)
    unknown_reads = tot_seqs - matched_reads
    oh.write('unknown_reads\t%s\n' % unknown_reads)
    oh.close()


def cleaning_up():
    '''sift reads into viral/unknown, compresses and removes files
    '''
    import multiprocessing as mp
    from Bio.SeqIO.QualityIO import FastqGeneralIterator

    # selects reads with coverage and identity higher than 75
    df = pd.read_csv('unique.tsv', sep='\t')
    viral_ids = set(df[(df.qcovs > blast_cov_threshold) & (df.pident > blast_ident_threshold)].qseqid)
    viral_c = 0
    undet_c = 0
    all_reads = 'hq_decont_reads.fastq'
    all_handle = open(all_reads)
    undet_handle = open('undetermined_reads.fastq', 'w')
    viral_handle = open('viral_reads.fastq', 'w')
    # Using FastqGeneralIterator allows fast performance
    for title, seq, qual in FastqGeneralIterator(all_handle):
        if title.split()[0] not in viral_ids:
            undet_c += 1
            undet_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
            if undet_c % 100000 == 0:
                logging.debug('written %d undet reads' % undet_c)
        else:
            viral_c += 1
            viral_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
            if viral_c % 10000 == 0:
                logging.debug('written %d viral reads' % viral_c)
    undet_handle.close()
    viral_handle.close()
    logging.info('written %d undet reads' % undet_c)
    logging.info('written %d viral reads' % viral_c)

    run_child('gzip', '-f viral_reads.fastq')
    run_child('gzip', '-f undetermined_reads.fastq')
    os.remove(all_reads)

    cmls = []
    for samfile in glob.glob('*.sam'):
        stem = os.path.splitext(samfile)[0]
        cont = stem.split('_')[-1]
        if cont == 'ref':  # hack because _ in bovine file name
            cont = 'bt_ref'
        cml = 'sort -O bam -l 0 -T /tmp %s | \
        samtools view -T %s -C -o %s.cram -' % (samfile, ref_map[cont], stem)
        cmls.append(('samtools', cml))

    # run in parallel
    pool = mp.Pool()
    results = pool.map(single_process, cmls)
    for r in results:
        logging.info(r)

    for samfile in glob.glob('*.sam'):
        os.remove(samfile)


def main(args):
    ''''''

    if args.run:
        miseq_dir = args.run.rstrip('/')
        run_name = os.path.split(miseq_dir)[1]
        if run_name.startswith('1'):
            try:
                run_date, machine_name = run_name.split('_')[:2]
                logging.info('running on run %s from machine %s' % (run_name, machine_name))
            except ValueError:
                logging.info('running on directory %s' % miseq_dir)
            bc_dir = os.path.join(miseq_dir, 'Data/Intensities/BaseCalls/')
        else:
            bc_dir = miseq_dir

        rel_fastq_files = glob.glob('%s/*_S*.fastq*' % bc_dir)
        samples_to_run = [os.path.split(fq)[1].split('_')[1] for fq in rel_fastq_files]
        logging.info('samples to run: %s' % ' '.join(samples_to_run))
        all_fastq_files = [os.path.abspath(f) for f in rel_fastq_files]
    elif args.file:
        all_fastq_files = [os.path.abspath(args.file)]
        logging.info('running on a single file %s' % all_fastq_files[0])
        run_name = args.file.split('.')[0]

    out_dir = 'virmet_output_%s' % run_name
    try:
        os.mkdir(out_dir)
    except OSError:
        logging.error('directory %s exists' % out_dir)
    os.chdir(out_dir)

    # run hunter on all fastq files
    s_dirs = []
    for fq in all_fastq_files:
        logging.info('running hunter on %s' % fq)
        sd = hunter(fq)
        s_dirs.append(sd)

    # run mapping against contaminants to remove
    cont_reads = 'good.fastq'  # first run on good.fastq
    for cont in contaminant_db:
        for sample_dir in s_dirs:
            os.chdir(sample_dir)
            decont_reads = victor(input_reads=cont_reads, contaminant=cont)

            os.chdir(os.pardir)
        cont_reads = decont_reads  # decontaminated reads are input for next round (equal across samples)

    logging.info('blasting against viral database')
    file_to_blast = cont_reads  # last output of victor is input for blast
    try:
        n_proc = min(os.cpu_count(), 12)
    except NotImplementedError:
        n_proc = 2
    logging.info('%d cores that will be used' % n_proc)

    for sample_dir in s_dirs:
        os.chdir(sample_dir)
        viral_blast(file_to_blast, n_proc)
        logging.info('sample %s blasted' % sample_dir)
        os.chdir(os.pardir)

    logging.info('summarising and cleaning up')
    for sample_dir in s_dirs:
        os.chdir(sample_dir)
        cleaning_up()
        os.chdir(os.pardir)
