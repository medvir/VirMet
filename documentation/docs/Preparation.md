# Preparing for a real run: `fetch` or `update`, and `index`

## Fetch

After installation, one needs to populate the database directory. By default
this will be `/data/virmet_databases_update`, but it can be modified by the users.
In order to populate this database, use the subcommand `fetch` as follows:

```
virmet fetch --viral n
virmet fetch --human
virmet fetch --bact_fungal
virmet fetch --bovine
```

This might take long, because databases are big.
Therefore, we strongly recommend running it inside a `tmux` session.

By default, the download of the viral database performs a
compression step. Users can avoid it with `--no_db_compression`.
However, we highly discourage to do that because subsequent analyses
will be much slower and the senstivity improvement during the 
viral classification process is minimal.

Remarkably, only the option `--viral` takes an argument: `n` for nucleotide and `p` for
the protein viral database. However, please note that only nucleotide sequences
are used currently by [Wolfpack](./Wolfpack.md), while the protein sequences 
are planned for future expansions of the tool.

Further information on the `fetch` subcommand can be obtained using `-h`:

```
virmet fetch -h
usage: virmet <command> [options]

options:
  -h, --help           show this help message and exit
  --viral {n,p}        viral [nucleic acids/proteins]
  --human              human
  --bact_fungal        bacterial and fungal(RefSeq)
  --bovine             bovine (Bos taurus)
  --no_db_compression  do not compress the viral database
  --dbdir [DBDIR]      path to store the new Virmet database
```

## Update

More and more genomes are uploaded to NCBI database every month.

VirMet provides a simple way to update the viral database without the need to download
all the genomes again. This can be done with the subcommand `update` as follows:

`virmet update --viral n --update_min_date YYYY/MM/DD`

Information on how to use the `update` subcomand can be obtained with `-h`:

```
virmet update -h
usage: virmet <command> [options] 

options:
  -h, --help            show this help message and exit
  --viral {n,p}         update viral [n]ucleic/[p]rotein
  --picked              PICKED       
                        file with additional sequences, one ACC per line
  --update_min_date     UPDATE_MIN_DATE
                        update viral database with sequences produced after date YYYY/MM/DD
  --no_db_compression   do not compress the viral database
  --dbdir [DBDIR]       path to store the updated Virmet database
```

In addition to choosing a specific date from which to look for reference sequences (`--update_min_date`),
users can also provide a *.txt* file with specific viral genomes to include in the database (`--picked`).
For that, they need to provide all the Accession Numbers (ACC) of these genomes
as they appear in the NCBI RefSeq database.
The .txt file should contain one Accession Number per row.

As mentioned for the `fetch` subcommand, `update` also allows users to skip
the compression step (`--no_db_compression`) of the viral database, although this is
highly discouraged. Users can also change the default database path (`--dbdir`) if 
their whole VirMet database is in another location.

Important: don't forget to index the viral database again once it has been updated.

## Index

After dowloading or updating the databases, it's always needed to index them.
The subcommand `virmet index` is used for that, and it can take multiple arguments.
Therefore, if users need to index the whole VirMet database, they can run: 

`virmet index --viral n --human --bact_fungal --bovine`

and wait for the indexing to finish.

Alternatively, if they need to index only the viral database, they can do it as follows:
`virmet index --viral n`

Further information on how to use the `index` subcommand can be obtained with `-h`:

```
virmet index -h
usage: virmet <command> [options]

options:
  -h, --help       show this help message and exit
  --viral {n,p}    make blast index of viral database
  --human          make bwa index of human database
  --bact_fungal    build kraken2 bacterial and fungal database
  --bovine         make bwa index of bovine database
  --dbdir [DBDIR]  path to store the indexed Virmet database
```

## Database structure

If everything works as expected, your database directory should have the following structure:

```
virmet_databases_update/
├── viral_nuccore/
│   ├── viral_database.fasta
│   ├── viral_accn_taxid.dmp
│   └── viral_seqs_info.tsv
├── human/
│   ├── fasta/
│   │   └── GRCh38.fasta.gz
│   └── bwa/
│       └── bwa_files
├── bovine/
│   ├── fasta/
│   │   └── ref_Bos_taurus.fasta.gz
│   └── bwa/
│       └── bwa_files
├── bact_fungi/
│   ├── library/
│   │   ├── bacteria
│   │   │   └── library.fna
│   │   └── fungi
│   │       └── library.fna
│   └── taxonomy/
│       └── taxdump.tar.gz
├── names.dmp.gz
└── nodes.dmp.gz
```
Please, note that the name of the database doesn't need to be necessarily
`virmet_databases_update`. This is the default name.

If you read until here, you are ready to use [VirMet Wolfpack](./Wolfpack.md)!

**Enjoy analysing your mNGS samples!**