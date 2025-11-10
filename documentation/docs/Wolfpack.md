# Running a virus scan: `wolfpack`

## Wolfpack

Wolfpack is the main subcommand of VirMet and it allows to perform
mNGS analyses of the given samples (virus scan).

Information on how to use it can be obtained with `-h`:

```
virmet wolfpack -h
usage: virmet wolfpack [options]

Options:
  -h, --help       show this help message and exit
  --run RUN        Miseq run directory
  --file FILE      single fastq file
  --dbdir [DBDIR]  path to find and use the Virmet database. Default: /data/virmet_databases/
  --nocovplot      do not make the covplots. Default: make them
  --noctrls        do not analyse ntc- or Undetermined samples. Default: analyse them
```

Wolfpack can be run on a single file (`--file`) or on a directory (`--run`).
If a directory is provided, it analyses all the FASTQ files in there.
Please, note that Wolfpack works if the structure of the `--run` directory is
any of these 2:

```
# Structure 1
run_directory/
├── FASTQ_file_1
├── FASTQ_file_2
├── ...
└── FASTQ_file_n

# Structure 2
run_directory/
└── Data/
    └── Intensities/
        └── BaseCalls/
            ├── FASTQ_file_1
            ├── FASTQ_file_2
            ├── ...
            └── FASTQ_file_n
```

Importantly, the extension of the FASTQ files must be
`.fastq` or `.fastq.gz`.

To run Wolfpack, users have the option to specify `--dbdir`, which is
the path to the VirMet database and must include the viral, human, bovine, 
bacterial and fungal (bact_fungi) databases.
See how to [prepare for a VirMet run](./Preparation.md) before running
Wolfpack.
If `-dbdir` is not specified, VirMet will assume that it is located in
`/data/virmet_databases/`.

Users also have the option to run the mNGS analyses on all the FASTQ files
(the default behaviour), which usually include negative controls and 
undertermined reads, or to skip those with `--noctrls`.
In this case, all files named `ntc-*` or `Undetermined*` will be ignored.

To carry out the virus scan, Wolfpack will first run a filtering step (QC)
based on quality, length and entropy of the given reads (named quality-control).
Then, it will perform a decontamination step, where reads of human, bovine,
bacterial and fungal origin will be discarded.
Finally, remaining reads will be _blasted_ against the viral database and
the viruses present in the mNGS samples will be identified.

Before providing the results, Wolfpack will (by default) automatically check
which of the identified viral organisms fulfill these criteria:

* ≥3 sequencing reads match this specific viral organism.
* The weighted mean coverage score for this viral organism is ≥10%.
* The viral organism is not a phage nor a non-human virus (host is human).

For all the viral organisms meeting these criteria, Wolfpack will generate
[Coverage plots](./Covplot.md), where users can see what fraction of the genome is covered
by the sequenncing reads. This will provide them with further evidence
supporting (or not) the presence of the organism in the sample.

If users want to avoid this last step and are not interested in getting the
coverage plots, they can disable such option with `--nocovplot`.

## Outputs

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

More information about the [Coverage plots](./Covplot.md) is provided in the corresponding
section.

## Hot run

A virus scan on a full MiSeq run typically lasts a few hours, most of which are spent
in the decontamination phase. Sometimes, after a run is completed, we would like to
run it again with a new viral database. In these cases, `wolfpack` would run skipping
the previous phases to save time. It relies on the presence of intermediate files that,
if present, signals the pipeline that a specific step must be skipped.

These are the rules (must be intended for each sample):

* If `fastp.html` exists, skip quality filtering.
* If `good_humanGRCh38.err` exists, skip the human reads decontamination.
* If `good_humanGRCh38_bt_ref.err` exists, skip the bovine reads decontamination.
* if `good_humanGRCh38_bt_ref_bact_fungi.err.err` exists, skip the 
decontamination bacterial and fungal database.

Blasting against viral database will always be performed. If both `viral_reads.fastq.gz`
and `undetermined_reads.fastq.gz` exist, their content will be copied into a file,
they will be removed, and this new file will be blasted against the viral database.

Therefore, if users change the viral database after a run has already been analyzed,
running `virmet wolfpack` again will skip the quality filtering and decontamination steps
and proceed directly with the viral assignment.

**Enjoy with your viral classification!**