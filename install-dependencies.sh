#!/bin/bash

#  NCBI edirect tools
wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz -O /tmp/edirect.tar.gz
tar -xvf /tmp/edirect.tar.gz

# prinseq
wget http://downloads.sourceforge.net/project/prinseq/standalone/prinseq-lite-0.20.4.tar.gz -O /tmp/prinseq-lite-0.20.4.tar.gz
tar -xvf /tmp/prinseq-lite-0.20.4.tar.gz

# STAR
wget https://github.com/alexdobin/STAR/archive/STAR_2.4.2a.tar.gz -O /tmp/STAR_2.4.2a.tar.gz
tar -xvf /tmp/STAR_2.4.2a.tar.gz
cd STAR-STAR_2.4.2a && make STAR

# pip
wget https://bootstrap.pypa.io/get-pip.py -O /tmp/get-pip.py
sudo python3 /tmp/get-pip.py

# nose
pip install nose
