# Preparing for a real run: `fetch` or `upload`, and `index`

## Fetch

After installation, one needs to populate the database directory. By default
this will be `/data/virmet_databases` and will occupy about 60 GB. In order to
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

## Update

More and more sequences are uploaded to NCBI database every month.

VirMet provides a simple way to update the viral database without the need to download
all the genomes again. This can be done with the subcommand
`update` as in the example

    virmet update --viral n

Similarly, new bacterial sequences can be added as well as fungal. Human
and bovine database, since they consist of a single organism, don't need to be
updated so often.

Don't forget to index the database again once it has been updated.

### Adding sequences manually

By adding the switch `--picked file_with_ids` users can add sequences by
writing their ids in a file, one per line.

## Index

After dowloading or updating the databases, it's always needed to index them.
The subcommand `virmet index` is used for that, and it can take multiple arguments, 
so you can run

    virmet index --viral n --human --bact --fungal --bovine

and wait for the indexing to finish.