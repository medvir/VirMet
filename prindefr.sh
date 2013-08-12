#!/bin/bash
FILEIN=$1
DEC_OUT_NAME=decon_out
FASTAFILE=processed.fasta


DB="/data/databases/rins_viral_database.fasta"
KMER=11  # default=11
GLOBAL=1  # 0:local, 1:global
MATCH=75


prinseq() { /usr/local/bin/prinseq-lite.pl "$@"; }
deconseq() { /home/ozagordi/Dropbox/Software/VirMet/deconseq.pl "$@"; }
list_organisms() { /home/ozagordi/Dropbox/Software/VirMet/list_organisms.py\
	 "$@"; }

tax_orgs() { /home/ozagordi/Dropbox/Software/VirMet/tax_orgs.py "$@"; }

# print the number of reads to analyze
echo 'Reads to analyze'
echo `wc -l $FILEIN | cut -f 1 -d " "` / 4 | bc

echo `date`
echo 'cleaning with seqtk'
seqtk trimfq $FILEIN | seqtk seq -L 75 - > intermediate.fastq
echo ''

echo `date`
echo 'cleaning with prinseq'
prinseq -fastq intermediate.fastq -lc_method entropy -lc_threshold 70 \
        -log prinseq.log -min_qual_mean 20 -ns_max_p 25 \
		-out_good ./good -out_bad ./bad
rm intermediate.fastq
echo ''

echo `date`
echo 'decontaminating from human, bacterial and bovine with deconseq'
deconseq -f good.fastq -dbs hsref,bact,bos -id $DEC_OUT_NAME -keep_tmp_files \
         -c 80 -i 80 &> deconseq.log
echo ''

# convert fastq to fasta
#CONV_STR="from Bio import SeqIO; rs = SeqIO.convert('${DEC_OUT_NAME}_clean.fq',\
#	 'fastq', '$FASTAFILE', 'fasta'); print str(rs) + ' reads converted'"
#python -c "$CONV_STR"

seqret ${DEC_OUT_NAME}_clean.fq fasta::$FASTAFILE

#echo `date`
#echo 'running fr-hit'
#fr-hit -T 28 -r 1 -k $KMER -g $GLOBAL -m $MATCH -o \
#	clean-r_1-k_${KMER}-m_${MATCH}-global_${GLOBAL}.tsv -d $DB -a $FASTAFILE \
#		&> fr-hit.log
#echo ''

echo `date`
echo 'running blast'
echo 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send' > results.tsv
blastn -query $FASTAFILE -db /data/databases/rins_viral_blast_db \
-perc_identity 75 -max_target_seqs 1 -num_threads 24 \
-outfmt  '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send' >> results.tsv

echo `date`
echo 'listing organisms'
tax_orgs results.tsv
#list_organisms clean-r_1-k_${KMER}-m_${MATCH}-global_${GLOBAL}.tsv > \
#	 orgs_list.csv 2> all_reads_id_orgs.log
echo ''
echo -e "\033[1;31m===========================================================\033[0m"
echo ''
