#!/bin/bash
FILEIN=$1
DEC_OUT_NAME=decon_out
FASTAFILE=processed.fasta


DB="/data/databases/rins_viral_database.fasta"
KMER=11  # default=11
GLOBAL=1  # 0:local, 1:global
MATCH=75


prinseq() { /usr/local/bin/prinseq-lite.pl "$@"; }
#fr-hit() { /home/ozagordi/Projects/3rd_party_SW/fr-hit-v0.7-x86_64/fr-hit "$@"; }
deconseq() { /home/ozagordi/Dropbox/Software/VirMet/deconseq.pl "$@"; }
list_organisms() { /home/ozagordi/Dropbox/Software/VirMet/list_organisms.py "$@"; }

# print the number of reads to analyze
echo 'Reads to analyze'
echo `wc -l $FILEIN | cut -f 1 -d " "` / 4 | bc

# cleaning with prinseq
prinseq -fastq $FILEIN -derep 1 -exact_only -lc_method entropy -lc_threshold 70 -graph_data -graph_stats ld,gc,qd,de -log prinseq.log -min_len 75 -min_qual_mean 20 -ns_max_p 25 -out_good ./good -out_bad ./bad

# decontaminate from human and bacterial with deconseq
deconseq -f good.fastq -dbs hsref,bact -id $DEC_OUT_NAME -keep_tmp_files &> deconseq.log

# convert fastq to fasta
CONV_STR="from Bio import SeqIO; rs = SeqIO.convert('${DEC_OUT_NAME}_clean.fq', 'fastq', '$FASTAFILE', 'fasta'); print str(rs) + ' reads converted'"
python -c "$CONV_STR"

# run fr-hit
fr-hit -T 28 -r 1 -k $KMER -g $GLOBAL -m $MATCH -o clean-r_1-k_${KMER}-m_${MATCH}-global_${GLOBAL}.tsv -d $DB -a $FASTAFILE &> fr-hit.log

# list organisms in a csv file
list_organisms clean-r_1-k_${KMER}-m_${MATCH}-global_${GLOBAL}.tsv > orgs_list.csv 2> all_reads_id_orgs.log
