#!/opt/python2.7/bin/python2.7

'''Decontaminate reads by aligning against contaminants with STAR'''

import sys
import os

import logging
import logging.handlers

nthreads = 12  # not too many processors, huge RAM requirements
multimapNmax = 30  # reads mapping more times are not reported

loc_dir = {'human': '/data/databases/Homo_sapiens',
           'bact1': '/data/databases/Bacteria/db1',
           'bact2': '/data/databases/Bacteria/db2',
           'bact3': '/data/databases/Bacteria/db3',
           'bos': '/data/databases/Bos_taurus',
           'dog': '/data/databases/Canis_lupus',
           'hiv': '/Users/ozagordi/References/HIV-STAR',
           'pcv': '/Users/ozagordi/References/PCV1-STAR'
           }

cont_desc = {'human': 'Homo sapiens genome',
             'bact1': 'Bacterial database, part 1',
             'bact2': 'Bacterial database, part 2',
             'bact3': 'Bacterial database, part 3',
             'bos': 'Bos taurus genome',
             'dog': 'Domestic dog genome'
             }

def run_child(exe_name, arg_string):
    '''use subrocess to run an external program with arguments'''
    import subprocess

    if not arg_string.startswith(' '):
        arg_string = ' ' + arg_string

    try:
        viclog.debug(exe_name + arg_string)
    except:
        print >> sys.stderr, exe_name + arg_string
    try:
        retcode = subprocess.call(exe_name + arg_string, shell=True)
        if retcode > 0:
            try:
                viclog.error("Child %s terminated by signal" % exe_name, retcode)
            except:
                pass
            print "Child %s terminated by signal" % exe_name, retcode
        else:
            try:
                viclog.debug("Child %s returned %i" % (exe_name, retcode))
            except:
                pass
            print "Child %s returned %i" % (exe_name, retcode)
    except OSError as ee:
        #viclog.error("Execution of %s failed:" % exe_name, ee)
        print "Execution of %s failed:" % exe_name, ee
    return retcode


def parse_com_line():
    '''Standard option parsing'''
    import argparse
    args = None

    # default action is 'store'
    parser = argparse.ArgumentParser(description='Decontaminate fastq file \
    according to a list of contaminant organisms',
                                     epilog='%(prog) -l for a list of contaminants')
    parser.add_argument('-r', '--readfile', dest='readfile',
                        help='fastq/fasta file with reads')
    parser.add_argument('-d', '--contaminant', dest='contaminant',
                        help='contaminant organism')
    parser.add_argument('-o', '--output', dest='output',
                        help='output file name stem')
    parser.add_argument('-c', '--coverage', dest='coverage', default=0.90,
                        help='coverage threshold required for a match')
    parser.add_argument('-i', '--identity', dest='identity', default=0.95,
                        help='coverage threshold required for a match')
    parser.add_argument('-l', '--list', dest='list', action='store_true',
                        help='list the available contaminants')
    
    args = parser.parse_args()
    if args.list:
        for k, v in cont_desc:
            print k, v
        sys.exit()
    return args


def remove_reads(sam_file, contaminated_reads):
    '''Remove reads found in sam file from the original fastq file.
    Returns the name of file with decontaminated reads'''

    import gzip
    from Bio.SeqIO.QualityIO import FastqGeneralIterator

    # test if an object is in set is way faster than in list
    mapped_reads = set([l.split()[0] for l in open(sam_file)
                        if not l.startswith('@')])

    clean_name = '.'.join(sam_file.split('.')[:-1]) + '.fastq'
    output_handle = open(clean_name, 'w')
    try:
        viclog.info('Cleaning reads in %s with alignments in %s' %
                    (contaminated_reads, sam_file))
        viclog.info('Writing to %s' % clean_name)
    except:
        print >> sys.stderr, 'Cleaning reads in %s with alignments in %s' % (contaminated_reads, sam_file)
        print >> sys.stderr, 'Writing to %s' % clean_name

    if contaminated_reads.endswith('.gz'):
        cont_handle = gzip.open(contaminated_reads)
    else:
        cont_handle = open(contaminated_reads)

    c = 0
    # Using FastqGeneralIterator allows fast performance
    for title, seq, qual in FastqGeneralIterator(cont_handle):
        if title.split()[0] not in mapped_reads:
            c += 1
            output_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
            if c % 100000 == 0:
                print >> sys.stderr, 'written %d clean reads' % c

    output_handle.close()

    return clean_name


def run_star(readfile, contaminant, overlap, identity):
    '''Runs the aligner STAR and returns the name of sam file'''

    exe_name = 'STAR'
    arg_string = '--genomeDir %s' % loc_dir[contaminant]
    arg_string += ' --readFilesIn %s' % readfile
    if readfile.endswith('.gz'):
        arg_string += ' --readFilesCommand zcat'
    arg_string += ' --runThreadN %d' % nthreads
    arg_string += ' --outFilterScoreMinOverLread %f' % overlap
    arg_string += ' --outFilterMatchNminOverLread %f' % identity
    arg_string += ' --outFilterMultimapNmax %d' % multimapNmax
    arg_string += ' --genomeLoad LoadAndKeep'
    run_child(exe_name, arg_string)

    rf_head = os.path.split(readfile)[1]
    if readfile.endswith('.gz'):
        rf_head = '.'.join(rf_head.split('.')[:-2])
    else:
        rf_head = '.'.join(rf_head.split('.')[:-1])

    sam_name = '%s_%s.sam' % (rf_head, contaminant)
    os.rename('Aligned.out.sam', sam_name)
    sj_name = 'SJ_%s.tab' % contaminant
    os.rename('SJ.out.tab', sj_name)
    os.rename('Log.final.out', 'Log_%s.final.out' % contaminant)
    try:
        viclog.info('Alignment written to %s' % sam_name)
    except:
        print >> sys.stderr, 'Alignment written to %s' % sam_name
    return sam_name


def main(input_reads='', contaminant='', coverage=0.90, identity=0.90,
         remove=False):

    # Make a global logging object.
    viclog = logging.getLogger(__name__)

    # set logging level
    viclog.setLevel(logging.DEBUG)
    # This handler writes everything to a file.
    LOG_FILENAME = './victor.log'
    hl = logging.handlers.RotatingFileHandler(LOG_FILENAME, 'w',
                                              maxBytes=100000, backupCount=5)
    f = logging.Formatter("%(levelname)s %(asctime)s %(funcName)s\
                          %(lineno)d %(message)s")
    hl.setFormatter(f)
    viclog.addHandler(hl)

    if remove:
        exe_name = 'STAR'
        arg_string = '--genomeDir %s --genomeLoad Remove' % loc_dir[contaminant] 
        run_child(exe_name, arg_string)
        return 'Genome from %s removed' % contaminant

    sam_cont = run_star(input_reads, contaminant, coverage,
                        identity)
    clean_name = remove_reads(sam_cont, input_reads)

    return clean_name


if __name__ == '__main__':
    cml_args = parse_com_line()

    # Make a global logging object.
    viclog = logging.getLogger(__name__)

    # set logging level
    viclog.setLevel(logging.DEBUG)
    # This handler writes everything to a file.
    LOG_FILENAME = './victor.log'
    hl = logging.handlers.RotatingFileHandler(LOG_FILENAME, 'w',
                                              maxBytes=100000, backupCount=5)
    f = logging.Formatter("%(levelname)s %(asctime)s %(funcName)s\
                          %(lineno)d %(message)s")
    hl.setFormatter(f)
    viclog.addHandler(hl)
    viclog.info(' '.join(sys.argv))

    input_reads = cml_args.readfile
    # for cont in cml_args.contaminants:
    viclog.info('Aligning %s againts %s' % (input_reads, cml_args.contaminant))
    sam_cont = run_star(input_reads, cml_args.contaminant, cml_args.coverage,
                        cml_args.identity)
    clean_name = remove_reads(sam_cont, input_reads)
    os.rename(input_reads, '%s_clean.fastq' % cml_args.output)
    print >> sys.stderr, 'Clean reads now in %s_clean.fastq' % cml_args.output
