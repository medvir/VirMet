#!/opt/python2.7/bin/python2.7
'''Runs hunter on all samples of a single MiSeq run'''
import sys
import os
import glob

hunter_exe = '/home/ozagordi/Dropbox/Software/VirMet/hunter.sh'

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
assert run_name.startswith('1'), 'Check your run directory'
run_date, machine_name = run_name.split('_')[:2]

assert run_date.startswith('1') and len(run_date) == 6
assert machine_name.startswith('M0')

bc_dir = os.path.join(miseq_dir, 'Data/Intensities/BaseCalls/')

all_fastq_files = glob.glob('%s/*_S*.fastq*' % bc_dir)

samples_to_run = set([])
for fq in all_fastq_files:
    sample_n = os.path.split(fq)[1].split('_')[1]
    samples_to_run.add(sample_n)
    
fq_to_run = []
for smp in samples_to_run:
    full_f = glob.glob('%s/*_%s_*fastq.gz' % (bc_dir, smp))
    if full_f == []:
        full_f = glob.glob('%s/*_%s_*fastq' % (bc_dir, smp))
    fq_to_run.append(full_f[0])


out_dir = os.path.join('/VirMetResults', run_name)
try:
    os.mkdir(out_dir)
except OSError:
    print 'directory %s exists' % out_dir
os.chdir(out_dir)

for fq in fq_to_run:
    sample_n = os.path.split(fq)[1].split('_')[1]
    try:
        os.mkdir(sample_n)
    except OSError:
        print 'directory %s exists' % sample_n
    os.chdir(sample_n)
    # run hunter here
    run_child(hunter_exe, '%s > hunter.log' % fq)
    # hunter finished
    print os.getcwd()
    os.chdir('../')
