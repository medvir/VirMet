#!/bin/bash

#  NCBI edirect tools
cd /tmp
perl -MNet::FTP -e \
  '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1); $ftp->login;
   $ftp->binary; $ftp->get("/entrez/entrezdirect/edirect.zip");'
unzip -u -q edirect.zip
rm edirect.zip
export PATH=$PATH:/tmp/edirect
./edirect/setup.sh

# prinseq
wget http://downloads.sourceforge.net/project/prinseq/standalone/prinseq-lite-0.20.4.tar.gz \
-O /tmp/prinseq-lite-0.20.4.tar.gz
tar -xvf /tmp/prinseq-lite-0.20.4.tar.gz
sudo install -p tmp/prinseq-lite-0.20.4/prinseq-lite.pl /usr/local/bin

# samtools 1.3
wget https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2 \
-O /tmp/samtools-1.3.tar.bz2
tar xvfj /tmp/samtools-1.3.tar.bz2
cd /tmp/samtools-1.3
make
sudo make prefix=/usr/local install

# NCBI blast+ 2.3.0

cd /tmp
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.3.0+-x64-linux.tar.gz
tar xvzfp ncbi-blast-2.3.0+-x64-linux.tar.gz
