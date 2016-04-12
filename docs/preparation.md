## `fetch` and `index`: preparing for a real run

### Fetching databases

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


### Indexing databases

The subcommand `virmet index` can take multiple arguments, so you can run

    virmet index --viral n --human --bact --fungal --bovine

and wait for the indexing to finish.
