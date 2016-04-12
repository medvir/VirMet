### Updating the database

More and more sequences are uploaded to NCBI database every month. The figure
shows the number of viral sequences with _complete genome_ in the title
that are submitted every month to NCBI ([code](https://gist.github.com/ozagordi/c1e1c4158ab4e94e4683)).

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
