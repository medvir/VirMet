#!/bin/bash
OWNDIR=$(dirname $0)
sift_reads() { $OWNDIR/sift_reads.py "$@"; }
FILEIN=$1
echo 'Analysing file' $FILEIN

# Counts the reads correctly even if the file is gzipped
if [[ $FILEIN == *fastq ]];
	then
	echo 'counting original reads'
	NREADS=`wc -l $FILEIN | cut -f 1 -d " "`
elif [[ $FILEIN == *fastq.gz ]];
	then
	echo 'counting reads gzipped'
	NREADS=`gunzip -c $FILEIN | wc -l | cut -f 1 -d " "`
fi

let "NREADS /= 4"
echo 'found' $NREADS

PASS_READS=`wc -l good.fastq | cut -f 1 -d " "`
let "PASS_READS /= 4"
echo $PASS_READS 'of which passed the quality filter'

echo 'writing stat files'
H_READS=$(grep -v '^@\w\w' good_human.sam | cut -f 1 | sort -u | wc -l)
BAC1_READS=$(grep -v '^@\w\w' good_human_bact1.sam | cut -f 1 | sort -u | wc -l)
BAC2_READS=$(grep -v '^@\w\w' good_human_bact1_bact2.sam | cut -f 1 | sort -u | wc -l)
BAC3_READS=$(grep -v '^@\w\w' good_human_bact1_bact2_bact3.sam | cut -f 1 | sort -u | wc -l)
BOS_READS=$(grep -v '^@\w\w' good_human_bact1_bact2_bact3_bos.sam | cut -f 1 | sort -u | wc -l)
DOG_READS=$(grep -v '^@\w\w' good_human_bact1_bact2_bact3_bos_dog.sam | cut -f 1 | sort -u | wc -l)
let "BAC_READS = BAC1_READS + BAC2_READS + BAC3_READS"

CLEAN_READS=`grep -c ">" clean_filtered_reads.fasta`
VIR_READS=`wc -l unique.tsv | cut -f 1 -d " "`
let "VIR_READS -= 1"
echo 'total_reads,passing_quality,from_human,from_bacteria,from_bos_taurus,from_canis,clean,matching_viral_db' > stats.csv
echo $NREADS,$PASS_READS,$H_READS,$BAC_READS,$BOS_READS,$DOG_READS,$CLEAN_READS,$VIR_READS >> stats.csv

let "DISC_READS = NREADS - PASS_READS"
let "UNKN_READS = CLEAN_READS - VIR_READS"
echo 'id,reads,category' > hier_stats.csv

echo discarded,$DISC_READS,low quality >> hier_stats.csv
echo human,$H_READS,contaminated >> hier_stats.csv
echo bacterial,$BAC_READS,contaminated >> hier_stats.csv
echo bovine,$BOS_READS,contaminated >> hier_stats.csv
echo canine,$DOG_READS,contaminated >> hier_stats.csv
echo unknown,$UNKN_READS,clean >> hier_stats.csv
echo viral,$VIR_READS,clean >> hier_stats.csv

echo 'treemap plot with R'
Rscript $OWNDIR/treemap.R hier_stats.csv
mv Rplots.pdf read_origin_treemap.pdf

echo 'sifting, cleaning and zipping'
sift_reads clean_filtered_reads.fastq results.tsv
gzip viral_reads.fastq
gzip undetermined_reads.fastq
rm good.fastq
rm bad.fastq good_*.fastq
rm clean_filtered_reads.fastq
rm splitted_clean*fasta

for SAMFILE in good_human good_human_bact1 good_human_bact1_bact2 \
	good_human_bact1_bact2_bact3 good_human_bact1_bact2_bact3_bos \
		good_human_bact1_bact2_bact3_bos_dog
do
	samtools view -Su ${SAMFILE}.sam | samtools sort - ${SAMFILE}_sorted
	rm ${SAMFILE}.sam
done
