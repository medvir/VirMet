#!/bin/bash

set -e

#  NCBI edirect tools
cd $HOME
if [ ! -d "$HOME/edirect/" ]; then
    perl -MNet::FTP -e \
      '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1); $ftp->login;
       $ftp->binary; $ftp->get("/entrez/entrezdirect/edirect.zip");';
    unzip -u -q edirect.zip;
    rm edirect.zip;
    ./edirect/setup.sh;
else
    echo 'Using cached directory';
fi
# cd edirect
# install -p econtact edirutil efilter elink entrez-phrase-search eproxy \
# espell ftp-cp join-into-groups-of sort-uniq-count-rank xtract xtract.Linux \
# eaddress edirect.pl efetch einfo enotify epost esearch esummary ftp-ls nquire \
# reorder-columns setup-deps.pl sort-uniq-count word-at-a-time xtract.pl /usr/local/bin

# prinseq
cd /tmp
wget http://downloads.sourceforge.net/project/prinseq/standalone/prinseq-lite-0.20.4.tar.gz \
-O /tmp/prinseq-lite-0.20.4.tar.gz
tar -xvf /tmp/prinseq-lite-0.20.4.tar.gz
#install -p /tmp/prinseq-lite-0.20.4/prinseq-lite.pl /usr/local/bin

# samtools 1.3
wget https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2 \
-O /tmp/samtools-1.3.tar.bz2
tar xvfj /tmp/samtools-1.3.tar.bz2
cd /tmp/samtools-1.3 && make && make prefix=$HOME/samtools-1.3 install

# NCBI blast+ 2.3.0
cd /tmp
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.3.0+-x64-linux.tar.gz
tar xzfp ncbi-blast-2.3.0+-x64-linux.tar.gz
install -p ./ncbi-blast-2.3.0+/bin/* /usr/local/bin
