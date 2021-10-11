"""
Start the analysis
"""

import logging
import sys
import pkg_resources
from snakemake import snakemake
from snakemake.utils import available_cpu_count

logger = logging.getLogger(__name__)


class SnakemakeError(Exception):
    pass


def add_arguments(parser):
    # Positionals
    parser.add_argument(
        "targets", nargs="*", metavar="<TARGETS>",
        help="File(s) to create excluding paths). If omitted, the full pipeline is run."
    )
    # Options
    parser.add_argument(
        "-n", "--dryrun", default=False, action="store_true",
        help="Perform dry run of pipeline. DEFAULT: False."
    )
    parser.add_argument(
        "--dag", default=False, action="store_true",
        help="Print the dag in the graphviz dot language. DEFAULT: False. To det "
             "output to pdf file, pipe output into dot as follows: '$ dbsrep run --dag "
             "| dot -Tpdf > dag.pdf'"
    )
    parser.add_argument(
        "-c", "--cores", metavar="<int>", type=int,
        default=available_cpu_count(),
        help="Maximum number of cores to run in parallel. DEFAULT: Use as many as "
             "available."
    )
    parser.add_argument(
        '--keepgoing', '-k', default=False, action='store_true',
        help='If one job fails, finish the others.')
    parser.add_argument(
        '--unlock', default=False, action='store_true',
        help='Remove a lock on the working directory.'
    )
    parser.add_argument(
        "--dir",
        help="Path to analysis directory. DEFAULT: CWD"
    )


def main(args):
    targets = args.targets if args.targets else None
    try:
        run(args.dryrun, args.cores, args.keepgoing, args.unlock, args.dag, targets)
    except SnakemakeError:
        sys.exit(1)
    sys.exit(0)


def run(
        dryrun: bool = False,
        cores: int = 4,
        keepgoing: bool = False,
        unlock: bool = False,
        printdag: bool = False,
        targets=None,
        workdir=None,
):
    # snakemake sets up its own logging, and this cannot be easily changed
    # (setting keep_logger=True crashes), so remove our own log handler
    # for now
    logger.root.handlers = []
    snakefile_path = pkg_resources.resource_filename("lsvtool", "lsvtool.smk")
    success = snakemake(snakefile_path,
                        snakemakepath="snakemake",  # Needed in snakemake 3.9.0
                        dryrun=dryrun,
                        printdag=printdag,
                        quiet=False if not printdag else True,
                        cores=cores,
                        keepgoing=keepgoing,
                        unlock=unlock,
                        printshellcmds=True,
                        targets=targets,
                        log_handler=[print_log_on_error],
                        workdir=workdir)
    if not success:
        raise SnakemakeError()


def print_log_on_error(msg):
    """Prints logs of failed rules in case of error"""
    if msg["level"] == "job_error" and msg["log"]:
        for log in msg["log"]:
            head = f"=== Output from log: '{log}' ==="
            print(head)
            if log.exists:
                with open(log) as f:
                    print(f.read().strip())
            print("-"*len(head))