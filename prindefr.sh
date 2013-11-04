#!/bin/bash
FILEIN=$1
DEC_OUT_NAME=decon_out
FASTAFILE=processed.fasta

OWNDIR=$(dirname $(readlink -f "$BASH_SOURCE"))

DB="/data/databases/rins_viral_database.fasta"
KMER=11  # default=11
GLOBAL=1  # 0:local, 1:global
MATCH=75


prinseq() { /usr/local/bin/prinseq-lite.pl "$@"; }
deconseq() { $OWNDIR/deconseq.pl "$@"; }
tax_orgs() { $OWNDIR/tax_orgs.py "$@"; }

echo 'Analysing file' $1

# number of reads to analyse and lines for the splitting into 16 processes
# integer division so that it stays divisible by 4
NREADS=`wc -l $FILEIN | cut -f 1 -d " "`
let "NREADS /= 4"
echo 'Reads to analyze:' $NREADS
let "MAX_L = NREADS / 16"
let "MAX_L *= 4"
#echo $MAX_L

echo `date`
echo 'cleaning with seqtk'
seqtk trimfq $FILEIN | seqtk seq -L 75 - > intermediate.fastq
echo ''

split -l $MAX_L -d intermediate.fastq splitted
find . -name "splitted*" | xargs -I % mv % %.fastq
rm intermediate.fastq


echo `date`
echo 'cleaning with prinseq'
seq -w 0 15 | xargs -P 0 -I {} /usr/local/bin/prinseq-lite.pl \
	    -fastq splitted{}.fastq -lc_method entropy -lc_threshold 70 \
        -log prinseq{}.log -min_qual_mean 20 -ns_max_p 25 \
		-out_good ./good{} -out_bad ./bad{}
echo ''

cat good??.fastq > good.fastq
cat bad??.fastq > bad.fastq
cat prinseq??.log > prinseq.log
seq -w 0 15 | xargs -I {} rm splitted{}.fastq good{}.fastq bad{}.fastq \
	prinseq{}.log


echo `date`
echo 'decontaminating from human, bacterial and bovine with deconseq'
deconseq -f good.fastq -dbs hsref,bact,bos -id $DEC_OUT_NAME -keep_tmp_files \
         -c 90 -i 94 &> deconseq.log
echo ''

seqret ${DEC_OUT_NAME}_clean.fq fasta::$FASTAFILE

echo `date`
echo 'running blast'
blastn -query $FASTAFILE -db /data/databases/rins_viral_blast_db \
-perc_identity 75 -max_target_seqs 1 -num_threads 24 \
-outfmt  '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send' > results.tsv

echo `date`
echo 'extracting unique high-scoring segment pairs (HSP)'
echo 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send' > unique.tsv
sort -k1,1 -u results.tsv >> unique.tsv

echo `date`
echo 'listing organisms'
tax_orgs unique.tsv

echo `summary statistics`
echo 'total_reads,passing_quality,from_human,from_bacteria,from_bos_taurus,clean,matching_viral_db' > stats.csv

PASS_READS=`wc -l good.fastq | cut -f 1 -d " "`
let "PASS_READS /= 4"

H_READS=`cut -f 1 decon_out_hsr*tsv | sort | uniq | wc -l | cut -f 1 -d " "`
BAC_READS=`cut -f 1 decon_out_bac*tsv | sort | uniq | wc -l | cut -f 1 -d " "`
BOS_READS=`cut -f 1 decon_out_bos*tsv | sort | uniq | wc -l | cut -f 1 -d " "`

CLEAN_READS=`grep -c ">" processed.fasta`
VIR_READS=`wc -l unique.tsv | cut -f 1 -d " "`
let "VIR_READS -= 1"
echo $NREADS,$PASS_READS,$H_READS,$BAC_READS,$BOS_READS,$CLEAN_READS,$VIR_READS >> stats.csv

echo 'cleaning and zipping'
rm good.fastq
gzip bad.fastq
gzip decon_out_cont.fq
gzip decon_out_clean.fq

echo ''
echo -e "\033[1;31m===========================================================\033[0m"
echo ''
