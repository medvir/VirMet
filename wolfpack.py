#!/opt/python2.7/bin/python2.7
'''Runs hunter on all samples of a single MiSeq run'''
import sys
import os
import os.path
import glob
import victor

contaminants = ['human', 'bact1', 'bact2', 'bact3', 'bos', 'dog']

dn = os.path.abspath(os.path.dirname(__file__))
hunter_exe = os.path.join(dn, 'hunter.sh')
viral_blast_exe = os.path.join(dn, 'viral_blast.sh')
summarise_exe = os.path.join(dn, 'summarise.sh')

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


miseq_dir = sys.argv[1].rstrip('/')
run_name = os.path.split(miseq_dir)[1]
if run_name.startswith('1'):
    run_date, machine_name = run_name.split('_')[:2]
    assert run_date.startswith('1') and len(run_date) == 6
    assert machine_name.startswith('M0')
    bc_dir = os.path.join(miseq_dir, 'Data/Intensities/BaseCalls/')
else:
    bc_dir = miseq_dir
all_fastq_files = glob.glob('%s/*_S*.fastq*' % bc_dir)

samples_to_run = set([])
for fq in all_fastq_files:
    sample_n = os.path.split(fq)[1].split('_')[1]
    samples_to_run.add(sample_n)

print 'samples to run: ', ' '.join(samples_to_run)

fq_to_run = []
for smp in samples_to_run:
    fgz_patt = '%s/*_%s_*fastq.gz' % (bc_dir, smp)
    full_f = glob.glob(fgz_patt)
    if full_f == []:
        full_f = glob.glob('%s/*_%s_*fastq' % (bc_dir, smp))
    fq_to_run.append(os.path.abspath(full_f[0]))

#out_dir = os.path.join('/data/VirMetResults', run_name)
out_dir = os.path.join('./', 'wolf_out_%s' % run_name)
try:
    os.mkdir(out_dir)
except OSError:
    print 'directory %s exists' % out_dir
os.chdir(out_dir)


# Run hunter stopping at quality

for fq in fq_to_run:
    sample_n = os.path.split(fq)[1].split('_')[1]
    try:
        os.mkdir(sample_n)
    except OSError:
        print 'directory %s exists' % sample_n
    os.chdir(sample_n)
    print sample_n
    run_child(hunter_exe, ' %s quality > hunter.log' % fq)
    # hunter finished
    os.chdir('../')


# Run decontamination with STAR loading each genome first

# initial input file is done on good.fastq, then good_human.fastq and so on
cont_reads = 'good.fastq'

# contaminant loop
for cont in contaminants:
    # sample loop
    for fq in fq_to_run:
        sample_n = os.path.split(fq)[1].split('_')[1]
        os.chdir(sample_n)
        out1 = victor.main(input_reads=cont_reads, contaminant=cont)
        os.chdir('../')

    # saves reads file name for next STAR filter: same for all samples
    try:
        cont_reads = out1
    except NameError:
        print >> sys.stderr, 'Did not decontaminate'

    # at the end of mapping runs, remove genomve from memory
    rmout = victor.main(contaminant=cont, remove=True)
    os.remove('Log.progress.out')
    os.remove('Log.out')
    os.remove('Aligned.out.sam')
    print >> sys.stderr, rmout


# Now run viral_blast (blast against viral database) and summarise/clean
# first change name
for fq in fq_to_run:
    sample_n = os.path.split(fq)[1].split('_')[1]
    os.chdir(sample_n)
    print >> sys.stderr, os.getcwd()
    os.rename(cont_reads, 'clean_filtered_reads.fastq')
    run_child(viral_blast_exe, 'clean_filtered_reads.fastq')
    run_child(summarise_exe, fq)
    os.chdir('../')
