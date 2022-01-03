"""
Create and initialize a new analysis directory.
"""
import logging
import os
import os.path
import sys
from pathlib import Path
from typing import List, Any
import shutil

import pkg_resources

logger = logging.getLogger(__name__)

CONFIG_FILE_NAME = "parameters.config"
ACCEPTED_FILE_EXT = (".vcf", ".vcf.gz")

def add_arguments(parser):
    parser.add_argument(
        "-i", "--input-vcfs", type=Path, nargs='*', required=True,
        help="input vcf file(s) separated by spaces, to perform the filteration and intersection "
             "detected in the same directory"
    )
    parser.add_argument(
        "-o", "--out_directory", type=Path, required=True,
        help="New analysis directory to create"
    )
    parser.add_argument(
        "-s", "--sample-names", nargs="*", 
        help="sample names to use for each input file separated by spaces. Otherwise input file "
             "names are used."
    )


def main(args):
    init(vcfs = args.input_vcfs, directory = args.out_directory, sample_names = args.sample_names)


def init(vcfs: List[Path], directory: Path, sample_names: List[str] = None):
    if " " in str(directory):
        logger.error("The name of the analysis directory must not contain spaces")
        sys.exit(1)

    if len(vcfs) < 2:
        logger.error("Please provide more than one VCF.")
        sys.exit(1)
    
    if sample_names is None:
        sample_names = [vcf.name.replace(".vcf.gz", "").replace(".vcf", "") for vcf in vcfs]
    else:
        assert len(sample_names) == len(vcfs), "Provide one sample name for each input file!"
    
    create_and_populate_analysis_directory(vcfs, directory, sample_names)

    logger.info(f"Directory {directory} initialized.")
    logger.info(f"To start the analysis: 'cd {directory} && lsvtool run' ")


def create_and_populate_analysis_directory(vcfs: List[Path], directory: Path, sample_names: List[str]):
    try:
        directory.mkdir()
    except OSError as e:
        logger.error(e)
        sys.exit(1)

    # Copy default configs
    config_file = pkg_resources.resource_filename("lsvtool", CONFIG_FILE_NAME)
    shutil.copy(config_file, directory)

    for vcf, name in zip(vcfs, sample_names):
        fail_if_inaccessible(vcf)
        create_symlink(vcf, directory, name)


def fail_if_inaccessible(path: Path):
    try:
        with path.open():
            pass
    except OSError as e:
        logger.error("Could not open %r: %s", path, e)
        sys.exit(1)


def create_symlink(source: Path, dirname: Path, target: str):
    if not os.path.isabs(source):
        src = os.path.relpath(source, dirname)
    else:
        src = source

    # Check if file has the correct extension
    if not str(source).endswith(ACCEPTED_FILE_EXT):
        raise FileNotFoundError(f"File {source} is not accepted input.")

    ext = ".vcf.gz" if source.match("*.vcf.gz") else ".vcf"
    dest = dirname / f"{target}{ext}"

    logger.info(f"Creating symlink for file {source} --> {dest}")
    try:
        os.symlink(src, dest)
    except FileExistsError as e:
        logger.error("%s: Try to supply a custom names using -s/--sample-names", e)
        sys.exit(1) 


