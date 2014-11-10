#!/bin/bash
#
# This program takes a fastq file as input, then it runs several tools
# to filter, decontaminate and then blast reads to viral database.
# Finally, it writes files of taxonomy and summary statistics
# It can be called with option 'quality' (as from wolfpack.py), and it
# stops after quality filtering

VERSION=0.3.3
LOCKFILE="pipeline_version.lock"
FILEIN=$1
RUNTYPE=$2
# allowed runtype:
# - quality (stops after seqtk and prinseq)
# - if runtype is not set, full hunter is run
echo $FILEIN $RUNTYPE

DEC_OUT_NAME=clean_filtered_reads
FASTAFILE=clean_filtered_reads.fasta

#OWNDIR=$(dirname $(readlink -f "$BASH_SOURCE"))
OWNDIR=$(dirname $0)

KMER=11  # default=11
GLOBAL=1  # 0:local, 1:global
MATCH=75
THRESHDECONCOV=90
THRESHDECONID=94

prinseq() { /usr/local/bin/prinseq-lite.pl "$@"; }
deconseq() { $OWNDIR/deconseq.pl "$@"; }
tax_orgs() { $OWNDIR/tax_orgs.py "$@"; }
victor() { $OWNDIR/victor.py "$@"; }
sift_reads() { $OWNDIR/sift_reads.py "$@"; }

if [[ "$RUNTYPE" != "" && "$RUNTYPE" != "quality" ]]; then
		echo 'RUNTYPE can only be quality or undefined'
		exit
fi

if [ -f $LOCKFILE ] ; then
	echo 'pipeline lock file exists'
	exit
else
	echo "$VERSION" > $LOCKFILE
fi

echo 'Pipeline version' $VERSION
echo 'Analysing file' $1

# Counts the reads correctly even if the file is gzipped
if [[ $FILEIN == *fastq ]];
	then
	echo 'here'
	NREADS=`wc -l $FILEIN | cut -f 1 -d " "`
elif [[ $FILEIN == *fastq.gz ]];
	then
	echo 'there'
	NREADS=`gunzip -c $FILEIN | wc -l | cut -f 1 -d " "`
fi

let "NREADS /= 4"
echo 'Reads to analyze:' $NREADS
	
echo `date`
echo 'cleaning with seqtk'
seqtk trimfq $FILEIN | seqtk seq -L 75 - > intermediate.fastq
echo ''

INTREADS=`wc -l intermediate.fastq | sed 's/^ *//' | cut -f 1 -d " "`
let "INTREADS /= 4"
echo 'Reads passing seqtk:' $INTREADS
# We want to split in 16 processors, so each file has at most
# (NREADS / 16) + 1 reads and 4 times as many lines
let "MAX_READS_PER_FILE = INTREADS / 16"
let "MAX_READS_PER_FILE += 1"
let "MAX_L = MAX_READS_PER_FILE * 4"
echo 'Max lines per file is:' $MAX_L

#split -l $MAX_L -d intermediate.fastq splitted
split -l $MAX_L intermediate.fastq splitted
find . -name "splitted*" | xargs -I % mv % %.fastq
rm intermediate.fastq

c=0
for SPF in $(ls splitted*fastq)
do
	cc=`printf '%0.2d' $c`
	mv $SPF splitted${cc}.fastq
	let "c += 1"
done

echo `date`
echo 'cleaning with prinseq'
if [[ $HOSTNAME == "virologymc17.local" ]];
	then
	XARGS_THREAD=4
elif [[ $HOSTNAME == "virologysrv04.uzh.ch" ]];
	then
	XARGS_THREAD=0
fi
seq -w 0 15 | xargs -P $XARGS_THREAD -I {} /usr/local/bin/prinseq-lite.pl \
	    -fastq splitted{}.fastq -lc_method entropy -lc_threshold 70 \
        -log prinseq{}.log -min_qual_mean 20 -ns_max_p 25 \
		-out_good ./good{} -out_bad ./bad{}
echo ''

cat good??.fastq > good.fastq
cat bad??.fastq > bad.fastq
cat prinseq??.log > prinseq.log
seq -w 0 15 | xargs -I {} rm splitted{}.fastq good{}.fastq bad{}.fastq \
	prinseq{}.log

# wolfpack stops here
if [ -z "$RUNTYPE" ] ; then
	echo "Continuing to decontamination"
elif [ "$RUNTYPE" = "quality" ] ; then
	echo "Exiting after quality filter"
	exit
fi

echo `date`
echo 'decontaminating from human, bacterial, bovine with victor'
victor -r good.fastq -d human -d bact1 -d bact2 -d bact3 -d bos -d dog \
	-o ${DEC_OUT_NAME} &> victor.err.log
echo ''

echo `date`
echo `blasting against viral database and extracting best matches`
viral_blast.sh ${DEC_OUT_NAME}

echo `date`
echo 'summary statistics and clean up'
summarise.sh 

echo ''
echo -e "\033[1;31m===========================================================\033[0m"
echo ''
