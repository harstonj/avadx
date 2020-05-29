from __future__ import print_function
import os
import sys
import shutil
import argparse
import datetime
import builtins as __builtin__
from pathlib import Path


def print(*args, **kwargs):
    if not QUIET:
        # '{:23s}'.format(f'[{datetime.datetime.now().timestamp()}]')
        return __builtin__.print(*args, **kwargs)
    else:
        return __builtin__.print(end='')


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-a', '--action', required=True)
    parser.add_argument('-f', '--inputfile', action='append',
                        help="input file(s)")
    parser.add_argument('-o', '--outputfolder',
                        help="path to base output folder; default: INPUTFILE_out")
    parser.add_argument('-t', '--threads', default=1, type=int,
                        help="number of threads; default: 1")
    parser.add_argument('-c', '--cpu', type=int,
                        help="max cpus per thread; default: all available")
    parser.add_argument('-p', '--preserve', action='store_true',
                        help="if flag is set intermediate results are kept")
    parser.add_argument('-u', '--update', type=str,
                        help="update databases")
    parser.add_argument('-D', '--createdb', type=str, nargs="*", metavar="arg",
                        help="create new reference database: <db_name> <db_sequences.fasta> ")
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help="set verbosity level; default: log level INFO")
    parser.add_argument('-q', '--quiet', action='store_true',
                        help="if flag is set console output is logged to file")

    return parser.parse_args()

def main():
    args = parse_arguments()
    pass

if __name__ == "__main__":
    main()
