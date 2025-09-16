# Welcome to VirMet

VirMet is a software suite designed to help users running viral metagenomics (mNGS) experiments: unspecific massively parallel NGS sequencing with the aim of discovering and characterizing the virus fraction of biological samples.

To see and download the code, visit our [GitHub](https://github.com/medvir/VirMet/).

VirMet is called with a command-subcommand syntax. All the possible subcommands are:

* [`fetch`](./Preparation.md): download genomes
* [`update`](./Preparation.md): update viral/bacterial database
* [`index`](./Preparation.md): index genomes
* [`wolfpack`](./Wolfpack.md): analyze a Miseq run
* [`covplot`](./Covplot.md): plot coverage for a specific organism

Some help can be obtained with `virmet <subcommand> -h`.

As an example, if you want to download the viral nucleotide database, you can use it as follows:
`virmet fetch --viral n`.

**Enjoy using VirMet!**

![](assets/logo.svg){: style="display:block; margin-left:auto; margin-right:auto; width:50%;" }
