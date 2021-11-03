"""
Create and initialize a new analysis directory.
"""
import logging
import os
import os.path
import sys
from pathlib import Path

import pkg_resources

logger = logging.getLogger(__name__)

ACCEPTED_FILE_EXT = (".bedpe", "bed")

def add_arguments(parser):
    parser.add_argument(
        "-i", "--input_bedpe", type=list, action='append', nargs='+',
        help="input bedpe file(s) separated by spaces, to perform the filteration and intersection "
             "detected in the same directory"
    )
    parser.add_argument(
        "-o", "--out_directory", type=Path,
        help="New analysis directory to create"
    )


def main(args):
    init(args.input_bedpe, args.out_directory)
    config_file = pkg_resources.resource_filename("lsvtool","/cli/parameters.config")
    os.system("cp {} {}/parameters.config".format(config_file, args.out_directory))
    
def init(bedpe: list, directory: Path):
    bedpe = [''.join([str(elem) for elem in i]) for i in bedpe[0]]

    if " " in str(directory):
        logger.error("The name of the analysis directory must not contain spaces")
        sys.exit(1)

    create_and_populate_analysis_directory(bedpe, directory)

    logger.info(f"Directory {directory} initialized.")
    logger.info(f"To start the analysis: 'cd {directory} && lsvtool run' ")


def create_and_populate_analysis_directory(bedpe: list, directory: Path):
    try:
        directory.mkdir()
    except OSError as e:
        logger.error(e)
        sys.exit(1)

    for bedpe_file in bedpe:
        print(bedpe_file)
        bedpe_file = Path(bedpe_file)
        fail_if_inaccessible(bedpe_file)
        create_symlink(bedpe_file, directory, bedpe_file)


def fail_if_inaccessible(path):
    try:
        with path.open():
            pass
    except OSError as e:
        logger.error("Could not open %r: %s", path, e)
        sys.exit(1)


def create_symlink(readspath, dirname, target):
    if not os.path.isabs(readspath):
        src = os.path.relpath(readspath, dirname)
    else:
        src = readspath

    # Check if file has the correct extension
    if not str(readspath).endswith(ACCEPTED_FILE_EXT):
        raise FileNotFoundError(f"File {readspath} is not accepted input.")

    os.symlink(src, os.path.join(dirname, target))


