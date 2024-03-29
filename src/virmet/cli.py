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
import os
import sys

from pkg_resources import (get_distribution, DistributionNotFound)

try:
    __version__ = get_distribution('virmet').version
except DistributionNotFound:
    # package is not installed
    pass

# manipulate path to import functions
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(1, parent_dir)
mod = __import__('virmet')
sys.modules["virmet"] = mod

# each subcommand takes one of these functions as default


def tidytable_run(args):
    """Default function for command line parser."""
    from virmet import tidytable
    tidytable.main(args)


def covplot_run(args):
    """Default function for command line parser."""
    from virmet import covplot
    covplot.main(args)


def wolfpack_run(args):
    """Default function for command line parser."""
    from shutil import move, Error
    from virmet import wolfpack
    od = wolfpack.main(args)
    try:
        om = move('virmet.log', od)
        sys.exit('logging file is in %s' % om)
    except FileNotFoundError:
        sys.exit('now in %s, file not found' % os.getcwd())
    except Error:
        sys.exit('logging file not moved after hot run')


def index_db(args):
    """Default function for command line parser."""
    from virmet import index
    index.main(args)


def update_db(args):
    """Default function for command line parser."""
    from virmet import update
    update.main(args)


def fetch_db(args):
    """Default function for command line parser."""
    from virmet import fetch
    fetch.main(args)


def main():
    """Parse command line, run default functions."""
    import argparse
    # parse command line
    # create the top-level parser
    parser = argparse.ArgumentParser(usage='%(prog)s <command> [options]',
                                     epilog="Run `virmet subcommand -h` for more help")
    parser.add_argument('-v', '--version', action='version', version=__version__)

    subparsers = parser.add_subparsers(help='available sub-commands')

    # create the parser for command "fetch"
    parser_fetch = subparsers.add_parser('fetch', help='download genomes')
    parser_fetch.add_argument('--viral', choices='np',
                              help='viral [nucleic acids/proteins]')
    parser_fetch.add_argument('--human', help='human', action='store_true')
    parser_fetch.add_argument('--bact', help='bacterial (RefSeq)', action='store_true')
    parser_fetch.add_argument('--fungal', help='fungal (RefSeq)', action='store_true')
    parser_fetch.add_argument('--bovine', help='bovine (Bos taurus)', action='store_true')
    parser_fetch.add_argument('--no_db_compression', help='do not compress the viral database', action='store_true', default=False)
    parser_fetch.set_defaults(func=fetch_db)

    # create the parser for command "update"
    parser_update = subparsers.add_parser('update', help='update viral/bacterial database')
    parser_update.add_argument('--viral', choices='np', help='update viral [n]ucleic/[p]rotein', default=False)
    parser_update.add_argument('--bact', help='update bacterial database', action='store_true', default=False)
    parser_update.add_argument('--fungal', help='update fungal database', action='store_true', default=False)
    parser_update.add_argument('--picked', help='file with additional sequences, one GI per line', default=None)
    parser_update.add_argument('--update_min_date', help='update viral [n]ucleic/[p]rotein with sequences produced after the date YYYY/MM/DD (e.g. 2022/09/01)', default=None)
    parser_update.set_defaults(func=update_db)

    # create the parser for command "index"
    parser_index = subparsers.add_parser('index', help='index genomes')
    parser_index.add_argument('--viral', choices='np',
                              help='make blast index of viral database')
    parser_index.add_argument('--human', action='store_true',
                              help='make bwa index of human database')
    parser_index.add_argument('--bact', action='store_true',
                              help='make bwa index of bacterial database')
    parser_index.add_argument('--fungal', action='store_true',
                              help='make bwa index of fungal database')
    parser_index.add_argument('--bovine', action='store_true',
                              help='make bwa index of bovine database')
    parser_index.set_defaults(func=index_db)

    # create the parser for command "wolfpack"
    parser_wolf = subparsers.add_parser('wolfpack', help='analyze a Miseq run')
    parser_wolf.add_argument('--run', type=str, help='Miseq run directory')
    parser_wolf.add_argument('--file', type=str, help='single fastq file')
    parser_wolf.set_defaults(func=wolfpack_run)

    # create the parser for command "tidytable"
    parser_tidy = subparsers.add_parser('tidytable', help='make tables summarising the whole run')
    parser_tidy.add_argument('--outdir', type=str, help='path to run results directory (virmet_output_...)')
    parser_tidy.set_defaults(func=tidytable_run)

    # create the parser for command "covplot"
    parser_cov = subparsers.add_parser('covplot', help='plot coverage for a specific organism')
    parser_cov.add_argument('--outdir', type=str, help='path to sample results directory')
    parser_cov.add_argument('--organism', type=str, help='name of the organism as reported in orgs_list.tsv file')
    parser_cov.set_defaults(func=covplot_run)

    # exit so that log file is not written
    if len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
        # or sys.argv[1] == '-v' or sys.argv[1] == '--version':
        parser.print_help()
        sys.exit()

    # logging configuration
    import logging
    import logging.handlers
    logging.basicConfig(filename='virmet.log', level=logging.DEBUG,
                        format='%(levelname)s %(asctime)s %(filename)s: %(funcName)s() %(lineno)d: \t%(message)s',
                        datefmt='%Y/%m/%d %H:%M:%S')

    logging.info(' '.join(sys.argv))
    logging.info('VirMet version:%s', __version__)
    # parse the args
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":  # and __package__ is None:
    main()
