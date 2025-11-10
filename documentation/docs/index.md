# Welcome to VirMet

VirMet is a software suite designed to help users running viral metagenomics 
(mNGS) experiments: unspecific massively parallel NGS sequencing with the aim 
of discovering and characterizing the virus fraction of biological samples.

To see and download the code, visit our [GitHub](https://github.com/medvir/VirMet/).
Information on how to install VirMet can be found in the [Installation](./Installation.md) section.

VirMet is called with a command-subcommand syntax. All the possible subcommands are:

* [`fetch`](./Preparation.md): download genomes
* [`update`](./Preparation.md): update viral database
* [`index`](./Preparation.md): index genomes
* [`wolfpack`](./Wolfpack.md): analyze a Miseq run or file
* [`covplot`](./Covplot.md): plot coverage for a specific organism

Some help can be obtained with `virmet <subcommand> -h` or simply `virmet -h`:

```
virmet -h
usage: virmet <command> [options]

positional arguments:
  {fetch,update,index,wolfpack,covplot}
                        available sub-commands
    fetch               download genomes
    update              update viral database
    index               index genomes
    wolfpack            analyze a Miseq run
    covplot             create coverage plot

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

Run `virmet subcommand -h` for more help
```

**Enjoy using VirMet!**

![](assets/logo.svg){: style="display:block; margin-left:auto; margin-right:auto; width:50%;" }
