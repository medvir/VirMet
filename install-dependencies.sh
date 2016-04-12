#!/bin/bash
set -e

#  NCBI edirect tools
cd $HOME
#if [ ! -d "$HOME/edirect/" ]; then
    perl -MNet::FTP -e \
      '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1); $ftp->login;
       $ftp->binary; $ftp->get("/entrez/entrezdirect/edirect.zip");';
    unzip -u -q edirect.zip;
    rm edirect.zip;
    ./edirect/setup.sh;
#else
#    echo 'Using cache for edirect';
#fi

# prinseq
cd $HOME
#if [ ! -d "$HOME/prinseq/" ]; then
    mkdir prinseq;
    wget http://downloads.sourceforge.net/project/prinseq/standalone/prinseq-lite-0.20.4.tar.gz \
    -O /tmp/prinseq-lite-0.20.4.tar.gz;
    tar -xvf /tmp/prinseq-lite-0.20.4.tar.gz;
    cp /tmp/prinseq-lite-0.20.4/prinseq-lite.pl $HOME/prinseq;
#else
#    echo "Using cache for prinseq"
#fi

# samtools 1.3
#if [ ! -d "$HOME/samtools-1.3/" ]; then
    cd /tmp
    mkdir $HOME/samtools-1.3;
    wget https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2 \
    -O /tmp/samtools-1.3.tar.bz2;
    tar xvfj /tmp/samtools-1.3.tar.bz2;
    cd /tmp/samtools-1.3 && make && make prefix=$HOME/samtools-1.3 install;
#else
#    echo "Using cache for samtools"
#fi

# NCBI blast+ 2.3.0
#if [ ! -d "$HOME/ncbi-blast-2.3.0+/" ]; then
    cd $HOME;
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.3.0+-x64-linux.tar.gz;
    tar xzfp ncbi-blast-2.3.0+-x64-linux.tar.gz;
#else
#    echo "Using cache for NCBI blast+"
#fi
echo "TMPDIR is " $TMPDIR
