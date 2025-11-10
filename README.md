VirMet
------

VirMet is a software designed to help users running viral metagenomics 
(mNGS) experiments.

For full Documentation, please [Read the Docs](https://medvir.github.io/VirMet/).

VirMet is called with a command-subcommand syntax. All the possible subcommands are:

* [`fetch`](https://medvir.github.io/VirMet/Preparation/): download all databases
* [`update`](https://medvir.github.io/VirMet/Preparation/): update viral database
* [`index`](https://medvir.github.io/VirMet/Preparation/): index all genomes
* [`wolfpack`](https://medvir.github.io/VirMet/Wolfpack/): analyze a Miseq run or file
* [`covplot`](https://medvir.github.io/VirMet/Covplot/): plot coverage for a specific organism

Some help can be obtained with `virmet <subcommand> -h` or simply `virmet -h`:

```
virmet -h
usage: virmet <command> [options]

positional arguments:
  {fetch,update,index,wolfpack,covplot}
                        available sub-commands
    fetch               download databases
    update              update viral database
    index               index genomes
    wolfpack            analyze a Miseq run
    covplot             create coverage plot

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
```

Run `virmet subcommand -h` for more help.

### Installation

#### Bioconda

VirMet is available through [Bioconda](https://bioconda.github.io), a channel
for the [conda](http://conda.pydata.org/docs/intro.html) package manager. Once
conda is [installed](https://bioconda.github.io/#install-conda) and the
[channels](https://bioconda.github.io/#set-up-channels) are set up,
`conda install virmet` installs the package with all its dependencies.

### Preparation

VirMet contains programs to download and index the genome sequences,
instructions [here](https://medvir.github.io/VirMet/Preparation/).

### Running a virus scan

Once the `fetch` and `index` subcommands have ben run and the databases are downloaded and indexed, users can
use [`Wolfpack`](https://medvir.github.io/VirMet/Wolfpack/), the main subcommand of VirMet.

A very simple way to use Wolfpack is the following:
`virmet wolfpack --run path_to_run_directory`

With that, sequencing reads are filtered (quality-control), decontaminated, and finally blasted against a (large)
set of viral sequences. 

Virmet Wolfpack provides several outputs that can be found in the
output folder `virmet_output_RUN_NAME`, located at the current working directory.

For example, if users have a directory named `Exp01` with files:

```
exp_01/ABC-DNA_S9.fastq.gz
exp_01/DEF-RNA_S11.fastq.gz
```

and run ```virmet wolfpack --run Exp01```, their results will
appear into `virmet_output_Exp01`.

The most important output files are:

* **Orgs_species_found.csv**: table showing the viral organisms identified per
FASTQ file as well as the read counts matching each organism.
It looks as follows:

```
species                accn        reads  stitle	           ssciname	              covered_region  seq_len  sample	    run
Tunavirus T1	       MK213796.1  455	  Phage T1             Escherichia phage T1	  30152	          48836	   ABC-DNA_S9	Exp01
Varicellovirus a5	   KY559403.2  79	  Alphaherpesvirus 5   Alphaherpesvirus 5	  321	          137741   ABC-DNA_S9	Exp01
Lentivirus humimdef1   MT222957.1  13	  HIV-1 isolate X 	   HIV-1	              141	          9372	   DEF-RNA_S11  Exp01
Emesvirus zinderi	   EF204940.1  6265	  Phage MS2 isolate Y  Escherichia phage MS2  2836	          3569	   DEF-RNA_S11  Exp01
```

with each column meaning:

`species`: scientific name of the species corresponding to the database sequence.  
`accn`: accession number of the viral species corresponding to the database sequence.   
`reads`: number of reads assigned to this specific sequence.  
`stitle`: title of the sequence in the database (fasta header).  
`ssciname`: scientific name of the sequence.  
`covered_region`: number of nucleotides covered by at least one read.  
`seq_len`: length of the sequence.  
`sample`: name of the FASTQ file that lead to the results.  
`run`: name of the sequencing run or main folder.  

* **Run_reads_summary.tsv**: summary of all reads analyzed per sample, showing the
number of reads passing each step of the pipeline (e.g., QC, decontamination),
and the number of reads matching the human, bovine, bacterial, fungal and
viral database.
It looks as follows:

```
category	          reads	    sample	     run
raw_reads	          2656448	ABC-DNA_S9	 Exp01
trimmed_too_short	  46	    ABC-DNA_S9	 Exp01
low_entropy	          25	    ABC-DNA_S9	 Exp01
low_quality	          1081	    ABC-DNA_S9	 Exp01
passing_filter	      2655296	ABC-DNA_S9	 Exp01
matching_humanGRCh38  339020	ABC-DNA_S9	 Exp01
matching_bt_ref	      1319509	ABC-DNA_S9	 Exp01
matching_bacteria	  46024	    ABC-DNA_S9	 Exp01
matching_fungi	      10	    ABC-DNA_S9	 Exp01
matching_other_cells  7403	    ABC-DNA_S9	 Exp01
reads_to_blast	      943330	ABC-DNA_S9	 Exp01
viral_reads	          70579	    ABC-DNA_S9	 Exp01
undetermined_reads	  872751	ABC-DNA_S9	 Exp01
raw_reads	          2202164	DEF-RNA_S11	 Exp01
trimmed_too_short	  259441	DEF-RNA_S11	 Exp01
low_entropy	          18	    DEF-RNA_S11	 Exp01
low_quality	          1277	    DEF-RNA_S11	 Exp01
passing_filter	      1941428	DEF-RNA_S11	 Exp01
matching_humanGRCh38  526874	DEF-RNA_S11	 Exp01
matching_bt_ref	      819466	DEF-RNA_S11	 Exp01
matching_bacteria	  8721	    DEF-RNA_S11	 Exp01
matching_fungi	      21	    DEF-RNA_S11	 Exp01
matching_other_cells  299	    DEF-RNA_S11	 Exp01
reads_to_blast	      586047	DEF-RNA_S11	 Exp01
viral_reads	          168247	DEF-RNA_S11	 Exp01
undetermined_reads	  417800	DEF-RNA_S11	 Exp01
```

In addition to these main outputs, Wolfpack creates inside
`virmet_output_RUN_NAME` a subdirectory for each sample (FASTQ file).
There, users can find these additional files: 

* Unique.tsv.gz
* Viral_reads.fastq.gz
* Undetermined_reads.fastq.gz
* Orgs_list.tsv
* Stats.tsv
* Fastp.html
* **Folders named as viral organisms**
* Other.err

Files `orgs_list.tsv` and `stats.tsv` report the main output of the tool for
each sample. They are the same as `Orgs_species_found.csv` and `Run_reads_summary.tsv`,
respectively, but containing information only of one sample. Therefore, the
last two columns (_sample_ and _run_) are missing.

The `unique.tsv.gz` reports all BLAST hits to the viral database. These are
not the same as in `orgs_list.tsv` because VirMet further filters the
BLAST hits (`unique.tsv.gz`) to ensure that only those with ≥75% pident and
≥75% qcov are ultimately considered viral reads for diagnosis. 

As the names say, `viral_reads.fastq.gz` and `undetermined_reads.fastq.gz`
contain, respectively, reads identified as of viral origin and reads not
matching any of the considered genomes.

Finally, `fastp.html` shows the quality-control statistics and information
about the QC-filtering step.

In the _decontamination_ step, reads are aligned against the human genome first,
those matching are discarded while those not matching are aligned against
the bovine genome, and so on. In each step, some files ending with _.err_
are generated, and they can be used for inspection, if needed. However, they
do not provide valuable information for the users (unless there is any error)
and can be removed if desired.

Besides all these files, Wolfpack automatically creates (unless disabled with
`--nocovplot`) a few folders with the names of viral organisms. Such folders
contain the coverage plots of the viral assignments, named `Viral_organism_coverage.pdf`.
It is recommended to manually have a look at them before a final viral diagnosis.

Overall, users can expect the following structure for the outputs:

```
virmet_output_Exp01/
│
├── ABC-DNA_S9
│   ├── Virus_X/
│   │    └── Virus_X_coverage.pdf
│   ├── unique.tsv.gz, viral_reads.fastq.gz, undetermined_reads.fastq.gz
│   ├── orgs_list.tsv, stats.tsv, fastp.html
│   └── Others.err
│ 
├── DEF-RNA_S11
│   ├── Virus_Y/
│   │   └── Virus_Y_coverage.pdf
│   ├── Virus_Z/
│   │   └── Virus_Z_coverage.pdf
│   ├── unique.tsv.gz, viral_reads.fastq.gz, undetermined_reads.fastq.gz
│   ├── orgs_list.tsv, stats.tsv, fastp.html
│   └── Others.err
│
├── run_reads_summary.tsv
└── orgs_species_found.tsv
```

Please, see [VirMet Documentation](https://medvir.github.io/VirMet/) for a more extensive
explanation on how to use `fetch`, `index`, `wolfpack` and `covplot` subcommands.