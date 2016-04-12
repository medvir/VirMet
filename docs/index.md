# VirMet: a set of tools for viral metagenomics

VirMet is a software suite designed to help users running viral metagenomics
experiments: unspecific massively parallel sequencing with the aim of
discovering and characterizing the viruses present in biological samples.



virmet is now called with a command subcommand
syntax: `virmet fetch --viral n`, for example, downloads the bacterial
database. Other available subcommands so far are

- `fetch`               download genomes
- `update`              update viral/bacterial database
- `index`               index genomes
- `wolfpack`            analyze a Miseq run

A short help is obtained with `virmet subcommand -h`.
