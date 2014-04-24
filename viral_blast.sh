#!/bin/bash

OWNDIR=$(dirname $0)
tax_orgs() { $OWNDIR/tax_orgs.py "$@"; }

FILEIN=$1
FASTAFILE=clean_filtered_reads.fasta


seqret -auto $FILEIN fasta::$FASTAFILE

export BLASTDB="/data/databases"
echo `date`
echo 'splitting input file'


FASTAREADS=`grep -c ">" $FASTAFILE`
echo 'Reads to blast:' $FASTAREADS
# We want to split in 16 processors, so each file has at most
# (FASTAREADS / 16) + 1 reads
let "MAX_N = FASTAREADS / 16"
let "MAX_N += 1"

awk -v "MAX_N=$MAX_N" 'BEGIN {n_seq=0;} /^>/ {if(n_seq%MAX_N==0){file=sprintf("splitted_clean_%d.fasta", n_seq/MAX_N);} print >> file; n_seq++; next;} { print >> file; }' $FASTAFILE


echo `date`
echo 'running blast in parallel'
if [[ $HOSTNAME == "virologymc17.local" ]];
	then
	XARGS_THREAD=4
elif [[ $HOSTNAME == "virologysrv04.uzh.ch" ]];
	then
	XARGS_THREAD=0
fi
seq 0 7 | xargs -P $XARGS_THREAD -I {} blastn -task megablast \
	    -query splitted_clean_{}.fasta -db /data/databases/viral_db \
		-perc_identity 75 -max_target_seqs 1 -out tmp_{}.tsv \
		-outfmt '6 qseqid sseqid sscinames stitle pident qcovs score length mismatch gapopen qstart qend sstart send staxids' 
echo ''

cat tmp_*.tsv > tmp.tsv
seq 0 7 | xargs -I {} rm splitted_clean_{}.fasta tmp_{}.tsv

echo -e 'qseqid\tsseqid\tsscinames\tstitle\tpident\tqcovs\tscore\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tstaxids' > unique.tsv
# save the header
cp unique.tsv results.tsv

sort -k1,1 -u tmp.tsv >> unique.tsv

cat tmp.tsv >> results.tsv
rm tmp.tsv

echo `date`
echo 'listing organisms'
tax_orgs unique.tsv

