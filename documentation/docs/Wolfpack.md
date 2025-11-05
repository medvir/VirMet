# Running a virus scan: `wolfpack`

## Wolfpack

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

    virmet wolfpack --run exp_01

and, after some time, find the results in `virmet_output_exp01`. Many files are
present, the most important ones being `orgs_list.csv` and `stats.tsv`. The
first lists the viral organisms found with a count of reads that could be
matched to them.

    [user@host test_virmet]$ cat virmet_output_exp_01/AR-1_S1/orgs_list.tsv
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


### Additional files

At the end of a run a directory for each sample (fastq file analyzed) is
created containing the following files:

    good_humanGRCh38_bact1_bact2_bact3_fungi1_bt_ref.cram
    good_humanGRCh38_bact1_bact2_bact3_fungi1_bt_ref.err
    ...
    good_humanGRCh38_bact1.cram
    good_humanGRCh38_bact1.err
    good_humanGRCh38.cram
    good_humanGRCh38.err

    orgs_list.tsv
    prinseq.err
    prinseq.log
    stats.tsv
    undetermined_reads.fastq.gz
    unique.tsv.gz
    viral_reads.fastq.gz

Files `orgs_list.tsv` and `stats.tsv` report the main output of the tool as
reported above, while `unique.tsv.gz` reports blast hits to viral database.

As the names say, `viral_reads.fastq.gz` and `undetermined_reads.fastq.gz`
contain, respectively, reads identified as of viral origin and reads not
matching any of the considered genomes.

`prinseq.err` and `prinseq.log` are, respectively, the standard error and log
file of prinseq, used to filter reads. By inspecting this log file, VirMet
determines how many reads were discarded because of low entropy or low quality.

In the _decontamination_ step, reads are aligned against human genome first,
those matching are discarded while those not matching are aligned against
the first set of bacterial genomes, and so on.
File `good_humanGRCh38.cram` is the alignment of high quality reads (good) to
human genome, saved in CRAM format. File `good_humanGRCh38_bact1.cram` contains
the alignment to bacterial genomes in set bact1 of high quality reads (good) minus
those that were identified as matching human genome, and so on. File ending in
`err` contain the standard error of the conversion bam -> cram.


### Hot run

A virus scan on a full MiSeq run typically lasts a few hours, many of which are spent
in the decontamination phase. Sometimes, after a run is completed, we would like to
run it again with a new viral database. In these cases, `wolfpack` would run skipping
the previous phases to save time. It relies on the presence of intermediate files that,
if present, signals the pipeline that a specific step must be skipped.

These are the rules (must be intended for each sample):

- if both `prinseq.log` and `prinseq.err` exist, skip quality filtering,
- if `good_humanGRCh38.err` exists, skip the human reads decontamination,
- if `good_humanGRCh38_bact1.err` exists, skip the decontamination against first bacterial
  database, and so on.

Blasting against viral database will always be performed. If both `viral_reads.fastq.gz`
and `undetermined_reads.fastq.gz` exist, their content will be copied into a file,
they will be removed, and this new file will be blasted against the viral database.

In short, if we change the viral database after a run has already been analyzed, simply
running `virmet wolfpack` again will skip the quality filtering and go straight to
blast against viral database.