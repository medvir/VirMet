#!/bin/bash
set -e

# seqtk
if [ ! -e "$HOME/seqtk/seqtk" ]; then
    echo "Compiling and installing seqtk"
    if [ ! -d "$HOME/seqtk" ]; then
        mkdir $HOME/seqtk
    fi
    cd /tmp
    wget https://github.com/lh3/seqtk/archive/v1.1.tar.gz \
    -O seqtk-1.1.tar.gz;
    tar xvfz seqtk-1.1.tar.gz;
    cd seqtk-1.1 && make && install -v /tmp/seqtk-1.1/seqtk $HOME/seqtk/;
else
   echo "Using cache for seqtk"
fi

# NCBI edirect tools
if [ ! -e "$HOME/edirect/efetch" ]; then
    echo "Installing edirect"
    cd $HOME
    perl -MNet::FTP -e \
      '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1); $ftp->login;
       $ftp->binary; $ftp->get("/entrez/entrezdirect/edirect.zip");';
    unzip -u -q edirect.zip;
    rm edirect.zip;
    ./edirect/setup.sh;
    ./edirect/efetch -version;
else
    echo 'Using cache for edirect';
fi

# prinseq
if [ ! -e "$HOME/prinseq/prinseq-lite.pl" ]; then
    if [ ! -d "$HOME/prinseq" ]; then
        mkdir $HOME/prinseq
    fi
    cd /tmp
    wget http://downloads.sourceforge.net/project/prinseq/standalone/prinseq-lite-0.20.4.tar.gz \
    -O /tmp/prinseq-lite-0.20.4.tar.gz;
    tar -xvf prinseq-lite-0.20.4.tar.gz;
    install -v prinseq-lite-0.20.4/prinseq-lite.pl $HOME/prinseq;
else
    echo "Using cache for prinseq"
fi

# samtools 1.3
if [ ! -e "$HOME/samtools-1.3/bin" ]; then
    cd /tmp
    wget https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2 \
    -O samtools-1.3.tar.bz2;
    tar xvfj samtools-1.3.tar.bz2;
    cd samtools-1.3 && make && make prefix=$HOME/samtools-1.3 install;
else
   echo "Using cache for samtools"
fi

# NCBI blast+ 2.3.0
if [ ! -d "$HOME/ncbi-blast-2.3.0+/bin" ]; then
    cd $HOME;
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.3.0+-x64-linux.tar.gz;
    tar xzfp ncbi-blast-2.3.0+-x64-linux.tar.gz;
else
   echo "Using cache for NCBI blast+"
fi
