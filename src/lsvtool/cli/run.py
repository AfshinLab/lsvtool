"""
Run the pipeline.

This is a small wrapper around Snakemake that sets some default parameters.

To run full pipeline use:

    $ lsvtool run

Arguments for snakemake are passed as:

    $ lsvtool run [snakemake_args [snakemake_args ...]]]

For example, to do a dry-run (not executing anything, only showing what would be done) use:

    $ lsvtool run -n

For info about arguments related to Snakemake run '$ snakemake -h' or look at the official
documentation at https://snakemake.readthedocs.io/en/stable/executing/cli.html.
"""

# Snakemake wrapping parially based on:
#  - http://ivory.idyll.org/blog/2020-improved-workflows-as-applications.html
#  - https://github.com/dib-lab/charcoal/blob/latest/charcoal/__main__.py

import logging
import sys
import os
import subprocess
from typing import List
import pkg_resources

from snakemake.utils import available_cpu_count

logger = logging.getLogger(__name__)


def add_arguments(parser):
    parser.add_argument(
        '-c', '--cores', metavar='N', type=int, default=available_cpu_count(),
        help='Run on at most N CPU cores in parallel. '
        'Default: %(default)s (all available cores).'
    )

    # This argument will not capture any arguments due to nargs=-1. Instead parse_known_args()
    # is used in __main__.py to add any arguments not captured here to snakemake_args.
    smk_args = parser.add_argument_group("snakemake arguments")
    smk_args.add_argument(
        'snakemake_args', nargs=-1,
        help="Arguments passed to snakemake. For info about snakemake options run "
             "'snakemake --help'."
    )


def main(args):
    try:
        run(cores=args.cores,
            snakemake_args=args.snakemake_args)
    except subprocess.CalledProcessError as e:
        print(f'Error in snakemake invocation: {e}', file=sys.stderr)
        sys.exit(e.returncode)
    sys.exit(0)


def run(
    cores: int = 4,
    snakefile: str = "lsvtool.smk",
    workdir=None,
    snakemake_args: List[str] = None,
):   
    snakefile_path = pkg_resources.resource_filename("lsvtool", "lsvtool.smk")

    cmd = ["snakemake", "-s", str(snakefile_path), "--cores", str(cores)]

    # Set defaults
    cmd += ["--printshellcmds"]

    if workdir is not None:
        cmd += ["--directory", str(workdir)]

    if snakemake_args is not None:
        cmd += snakemake_args

    logger.debug(f"Command: {' '.join(cmd)}")
    subprocess.check_call(cmd)