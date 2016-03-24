VirMet
------

[![Build Status](https://travis-ci.org/ozagordi/VirMet.svg?branch=master)](https://travis-ci.org/ozagordi/VirMet)
(caution: coverage of the test is very low)

A set of tools for viral metagenomics.

virmet is now called with a command subcommand
syntax: `virmet fetch --viral n`, for example, downloads the bacterial
database. Other available subcommands so far are

- `fetch`               download genomes
- `update`              update viral/bacterial database
- `index`               index genomes
- `wolfpack`            analyze a Miseq run

A short help is obtained with `virmet subcommand -h`.

### Installation
VirMet relies on a number of third-party tools used to access databases, trim,
convert, filter and map reads. One can refer to the files [`.travis.yml`](./.travis.yml)
and [`install-dependencies.sh`](./install-dependencies.sh) for details or
further down in this README.
The dependencies are:

- bwa
- samtools 1.3
- tabix
- seqtk
- prinseq-lite
- edirect command line tools
- blast+ 2.3.0
- python (3.x, it's 2016...) with pandas and Biopython


### Preparation: fetching databases

After installation, one needs to populate the database directory. By default
this will be `\data\virmet_databases` and will occupy about 60 GB. In order to
populate this, use the subcommand `fetch`, for example as follows

    virmet fetch --viral n  # this downloads viral sequences, nucleotide only
    virmet fetch --human
    virmet fetch --bact
    virmet fetch --fungal
    virmet fetch --bovine

This might take long. If it's taking too long, you might want to include the
above commands in a `down.sh` file and run them overnight.

Only the option `--viral` takes an argument: `n` for nucleotide and `p` for
protein viral database. Currently only nucleotide sequences are used, while the
protein ones are foreseen as useful in discovery of novel viral sequences
(in a future version).


### Preparation: indexing databases

The subcommand `virmet index` can take multiple arguments, so you can run

    virmet index --viral n --human --bact --fungal --bovine

and wait for the indexing to finish.


### Running a virus scan

This can be run on a single file or on a directory. It will try to guess from
the naming scheme if it is a Miseq output directory (_i.e._ with
`Data/Intensities/BaseCalls/` structure) and analyze all fastq files in there.
The extension must be `.fastq` or `.fastq.gz`. It will then run a filtering
step based on quality, length and entropy (in short: reads with a lot of
repeats will be discarded), followed by a decontamination step where reads
of human/bacterial/bovine/fungal origin will be discarded. Finally, remaining
reads are _blasted_ against the viral database. The list of organisms with the
count of reads is in files `orgs_list.csv` in the output directory
(naming is `virmet_output_...`). For example, if we have a directory named
`exp_01` with files

    exp_01/AR-1_S1_L001_R1_001.fastq.gz
    exp_01/AR-2_S2_L001_R1_001.fastq.gz
    exp_01/AR-3_S3_L001_R1_001.fastq.gz
    exp_01/AR-4_S4_L001_R1_001.fastq.gz

we could run

    virmet wolfpack --dir exp_01

and, after some time, find the results in `virmet_output_exp01`. Many files are
present, the most important ones being `orgs_list.csv` and `stats.tsv`. The
first lists the viral organisms found with a count of reads that could be
matched to them.

    [user@host test_virmet]$ cat virmet_output_test_dir_150123/3-1-65_S5/orgs_list.tsv
    organism	reads
    Human adenovirus 7	126
    Human poliovirus 1 strain Sabin	45
    Human poliovirus 1 Mahoney	29
    Human adenovirus 3+11p	19
    Human adenovirus 16	1

The second file is a summary of all reads analyzed for this sample and how many
were passing a specific step of the pipeline or matching a specific database.

    [user@host test_virmet]$ cat virmet_output_exp01/AR-1_S1/stats.tsv
    raw_reads       6250
    trimmed_too_short       462
    low_entropy     1905
    low_quality     0
    passing_filter  3883
    matching_humanGRCh38    3463
    matching_bact1  0
    matching_bact2  0
    matching_bact3  0
    matching_fungi1 0
    matching_bt_ref 0
    reads_to_blast  420
    viral_reads     257
    undetermined_reads      163


### Updating the database

More and more sequences are uploaded to NCBI database every month. The figure
shows the number of viral sequences with _complete genome_ in the title
that are submitted every month to NCBI.

![Code used to create the figure is [here](https://gist.github.com/ozagordi/c1e1c4158ab4e94e4683)](./viral_genomes.png "NCBI complete viral genomes per month")

VirMet provides a simple way to update the viral database with the subcommand
`update` as in the example

    virmet update --viral n

Similarly, new bacterial sequences can be added as well as fungal. Bacterial
and bovine database, since they consist of a single organism, don't need to be
updated so often.

Don't forget to index the database again once it has been updated.

#### Adding sequences manually

By adding the switch `--picked file_with_ids` users can add sequences by
writing their ids in a file, one per line.

----

### Details on the installation
In essence, on a Ubuntu 14.04 one can run the following commands for a system wide
configuration.

    # system wide configuration available as Ubuntu packages
    sudo apt-get update -qq
    sudo apt-get install -qq -y build-essential ftp golang unzip \
    bwa tabix seqtk libwww-perl

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

    #
    sudo ln -sv /bin/uname /usr/bin/uname
    export PATH=/usr/local/bin:$PATH

Then, one needs python 3 (VirMet was mainly developed and tested on 3.4, but
any 3.x should work), together with [pandas](http://pandas.pydata.org) and
[Biopython](http://biopython.org/wiki/Main_Page). Go to the respective
installation pages and choose your favourite method. For continuous
integration on Travis we used conda (see [`.travis.yml`](./.travis.yml)).
