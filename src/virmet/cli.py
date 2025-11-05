#!/usr/bin/env python3
"""
Module that contains the command line app.

Why does this file exist, and why not put this in __main__?

  You might be tempted to import things from __main__ later, but that will cause
  problems: the code will get executed twice:

  - When you run `python -mminvar` python will execute
    ``__main__.py`` as a script. That means there won't be any
    ``minvar.__main__`` in ``sys.modules``.
  - When you import __main__ it will get executed again (as a module) because
    there's no ``minvar.__main__`` in ``sys.modules``.

  Also see (1) from http://click.pocoo.org/5/setuptools/#setuptools-integration
"""

import argparse
import logging
import os
import shutil
import sys

from virmet import fetch, index, update, wolfpack, covplot
from virmet.__init__ import __version__


def wolfpack_run(args):
    """Default function for command line parser."""
    od = wolfpack.main(args)
    try:
        om = shutil.move("virmet.log", od)
        sys.exit("logging file is in %s" % om)
    except FileNotFoundError:
        sys.exit("now in %s, file not found" % os.getcwd())
    except shutil.Error:
        sys.exit("logging file not moved after hot run")


def main():
    """Parse command line, run default functions."""
    # parse command line
    # create the top-level parser
    parser = argparse.ArgumentParser(
        usage="%(prog)s <command> [options]",
        epilog="Run `virmet subcommand -h` for more help",
    )
    parser.add_argument(
        "-v", "--version", action="version", version=__version__
    )

    subparsers = parser.add_subparsers(help="available sub-commands")

    # create the parser for command "fetch"
    parser_fetch = subparsers.add_parser("fetch", help="download genomes")
    parser_fetch.add_argument(
        "--viral", choices="np", help="viral [nucleic acids/proteins]"
    )
    parser_fetch.add_argument("--human", help="human", action="store_true")
    parser_fetch.add_argument(
        "--bact_fungal",
        help="bacterial and fungal(RefSeq)",
        action="store_true",
    )
    parser_fetch.add_argument(
        "--bovine", help="bovine (Bos taurus)", action="store_true"
    )
    parser_fetch.add_argument(
        "--no_db_compression",
        help="do not compress the viral database",
        action="store_true",
        default=False,
    )
    parser_fetch.add_argument(
        "--dbdir",
        type=str,
        nargs="?",
        default="/data/virmet_databases_update/",
        help="path to store the new Virmet database",
    )
    parser_fetch.set_defaults(func=fetch.main)

    # create the parser for command "update"
    parser_update = subparsers.add_parser(
        "update", help="update viral database"
    )
    parser_update.add_argument(
        "--viral",
        choices="np",
        help="update viral [n]ucleic/[p]rotein",
        default=False,
    )
    parser_update.add_argument(
        "--picked",
        help="file with additional sequences, one GI per line",
        default=None,
    )
    parser_update.add_argument(
        "--update_min_date",
        help="update viral [n]ucleic/[p]rotein with sequences produced after the date YYYY/MM/DD",
        default=None,
    )
    parser_update.add_argument(
        "--no_db_compression",
        help="do not compress the viral database",
        action="store_true",
        default=False,
    )
    parser_update.add_argument(
        "--dbdir",
        type=str,
        nargs="?",
        default="/data/virmet_databases_update/",
        help="path to store the updated Virmet database",
    )
    parser_update.set_defaults(func=update.main)

    # create the parser for command "index"
    parser_index = subparsers.add_parser("index", help="index genomes")
    parser_index.add_argument(
        "--viral", choices="np", help="make blast index of viral database"
    )
    parser_index.add_argument(
        "--human", action="store_true", help="make bwa index of human database"
    )
    parser_index.add_argument(
        "--bact_fungal",
        action="store_true",
        help="build kraken2 bacterial and fungal database",
    )
    parser_index.add_argument(
        "--bovine",
        action="store_true",
        help="make bwa index of bovine database",
    )
    parser_index.add_argument(
        "--dbdir",
        type=str,
        nargs="?",
        default="/data/virmet_databases_update/",
        help="path to store the indexed Virmet database",
    )
    parser_index.set_defaults(func=index.main)

    # create the parser for command "wolfpack"
    parser_wolf = subparsers.add_parser("wolfpack", help="analyze a Miseq run")
    parser_wolf.add_argument("--run", type=str, help="Miseq run directory")
    parser_wolf.add_argument("--file", type=str, help="single fastq file")
    parser_wolf.add_argument(
        "--dbdir",
        type=str,
        nargs="?",
        default="/data/virmet_databases",
        help="path to find and use the Virmet database",
    )
    parser_wolf.add_argument(
        "--nocovplot",
        action="store_true",
        help="do not make the covplots. Default: make them",
    )
    parser_wolf.add_argument(
        "--noctrls",
        action="store_true",
        help="do not analyse ntc- or Undetermined samples. Default: analyse them",
    )
    parser_wolf.set_defaults(func=wolfpack_run)

    # create the parser for command "covplot"
    parser_covplot = subparsers.add_parser(
        "covplot", help="create coverage plot"
    )
    parser_covplot.add_argument(
        "--outdir",
        type=str,
        nargs="?",
        default="./",
        help="path to store the coverage plots",
    )
    parser_covplot.add_argument(
        "--organism", default=None, help="ssciname of the organism of interest"
    )
    parser_covplot.add_argument(
        "--dbdir",
        type=str,
        nargs="?",
        default="/data/virmet_databases/",
        help="path to find and use the Virmet database",
    )
    parser_covplot.set_defaults(func=covplot.main)

    # exit so that log file is not written
    if len(sys.argv) == 1 or sys.argv[1] == "-h" or sys.argv[1] == "--help":
        # or sys.argv[1] == '-v' or sys.argv[1] == '--version':
        parser.print_help()
        sys.exit()

    # logging configuration
    logging.basicConfig(
        filename="virmet.log",
        level=logging.DEBUG,
        format="%(levelname)s %(asctime)s %(filename)s: %(funcName)s() %(lineno)d: \t%(message)s",
        datefmt="%Y/%m/%d %H:%M:%S",
    )

    logging.info(" ".join(sys.argv))
    logging.info("VirMet version:%s" % __version__)
    # parse the args
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
