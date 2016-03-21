VirMet
------

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
and [`install-dependencies.sh`](./install-dependencies.sh) for details. In
essence, on a Ubuntu 14.04 one can run the following commands for a system wide
configuration.

    # system wide configuration available as Ubuntu packages
    sudo apt-get update
    sudo apt-get install build-essential ncurses-dev cmake \
        bwa tabix ncbi-blast+ libwww-perl

    #  NCBI edirect tools
    cd /usr/local/
    perl -MNet::FTP -e \
       '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1); $ftp->login;
        $ftp->binary; $ftp->get("/entrez/entrezdirect/edirect.tar.gz");'
     tar xfz edirect.tar.gz
     rm edirect.tar.gz
     ./edirect/setup.sh

    # prinseq
    wget http://downloads.sourceforge.net/project/prinseq/standalone/prinseq-lite-0.20.4.tar.gz -O /usr/local/prinseq-lite-0.20.4.tar.gz
    tar -xvf /usr/local/prinseq-lite-0.20.4.tar.gz

    # samtools 1.3
    wget https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2 -O /tmp/samtools-1.3.tar.bz2
    tar xvfj /tmp/samtools-1.3.tar.bz2
    cd /tmp/samtools-1.3
    make
    sudo make prefix=/usr/local install

    export PATH=$PATH:/usr/local/edirect
    sudo ln -sv PATH=$PATH:/usr/local/prinseq-lite-0.20.4/prinseq-lite.pl /usr/local/bin/prinseq-lite.pl
    export PATH=/usr/local/bin:$PATH

Then, one needs python 3 (VirMet was mainly developed and tested on 3.4, but
any 3.x should work), together with [pandas](http://pandas.pydata.org) and
[Biopython](http://biopython.org/wiki/Main_Page). Go to the respective
installation pages and choose your favourite method. For continuous
integration on Travis we used conda (see [`.travis.yml`](./.travis.yml)).


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

    [user@host test_virmet]$ cat virmet_output_exp01/AR-1_S1/orgs_list.csv
    sscinames,0
    Human immunodeficiency virus 1,118
    Lactobacillus phage LF1,78
    Torque teno virus,43
    Human adenovirus 7,24

The second file is a summary of all reads analyzed for this sample and how many
were passing a specific step of the pipeline or matching a specific database.

    [user@host test_virmet]$ cat virmet_output_exp01/AR-1_S1/stats.tsv
    raw_reads	2721226
    trimmed_long	2425400
    high_quality	2161625
    matching_humanGRCh38	493252
    matching_bact1	430976
    matching_bact2	443481
    matching_bact3	44645
    matching_fungi1	1424
    matching_bt_ref	4536
    blasted_reads	743311
