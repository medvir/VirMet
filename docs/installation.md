## Installation

The classic `python setup.py install` should work, provided the user has the
necessary permission.

VirMet relies on a number of third-party tools used to access databases, trim,
convert, filter and map reads.

The tools VirMet dependends on are:

- bwa
- samtools 1.3
- tabix
- seqtk
- prinseq-lite
- edirect command line tools
- blast+ 2.3.0

Moreover, the following languages are used

- python (3.x, it's 2016...) with pandas and Biopython
- R (for `covplot` only)


### Commands to install dependencies on Ubuntu
On a Ubuntu 14.04 the following commands should provide a system wide
installation of the tools mentioned above plus R, although on Travis a slightly
different strategy is used.

    # system wide configuration available as Ubuntu packages
    sudo apt-get update -qq
    sudo apt-get install -qq -y build-essential ftp golang unzip \
    bwa tabix seqtk libwww-perl r-base

    #  NCBI edirect tools
    cd /tmp
    perl -MNet::FTP -e \
      '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1); $ftp->login;
       $ftp->binary; $ftp->get("/entrez/entrezdirect/edirect.zip");'
    unzip -u -q edirect.zip
    rm edirect.zip
    export PATH=$PATH:/tmp/edirect
    ./edirect/setup.sh
    cd edirect
    sudo install -p econtact edirutil efilter elink entrez-phrase-search eproxy \
    espell ftp-cp join-into-groups-of sort-uniq-count-rank xtract xtract.Linux \
    eaddress edirect.pl efetch einfo enotify epost esearch esummary ftp-ls nquire \
    reorder-columns setup-deps.pl sort-uniq-count word-at-a-time xtract.pl /usr/local/bin

    # prinseq
    cd /tmp
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
    tar xzfp ncbi-blast-2.3.0+-x64-linux.tar.gz
    sudo install -p ./ncbi-blast-2.3.0+/bin/* /usr/local/bin

    # edirect looks for uname in the wrong place
    sudo ln -sv /bin/uname /usr/bin/uname

    export PATH=/usr/local/bin:$PATH

Then, one needs python 3 (VirMet was mainly developed and tested on 3.4, but
any 3.x should work), together with [pandas](http://pandas.pydata.org) and
[Biopython](http://biopython.org/wiki/Main_Page). Go to the respective
installation pages and choose your favourite method. For continuous
integration on Travis we used conda (see [`.travis.yml`](./.travis.yml)).
Finally, R packages ggplot2 and ggthemes are needed to run `covplot`.
To install these, type `install.packages(c('ggplot2', 'ggthemes'))` inside the
R console.
