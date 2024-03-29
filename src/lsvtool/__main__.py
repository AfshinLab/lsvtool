
"""
lsvtool is a code for filtering and intersecting Large Structural Variants
from different technologies (NAIBR, LinkedSV, 10X and Pacbio)
"""

import sys
import logging
import pkgutil
import importlib
from argparse import ArgumentParser

import lsvtool.cli as cli_package

logger = logging.getLogger(__name__)


def main(commandline_arguments=None) -> int:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(module)s - %(levelname)s: %(message)s",
                        datefmt='%Y-%m-%d %H:%M:%S')
    parser = ArgumentParser(description=__doc__, prog="lsvtool")
    parser.add_argument("--version", action="version", version="%(prog)s 0.1")
    parser.add_argument("--debug", action="store_true", default=False,
                        help="Print debug messages")
    subparsers = parser.add_subparsers()

    # Import each module that implements a subcommand and add a subparser for it.
    # Each subcommand is implemented as a module in the cli subpackage.
    # It needs to implement an add_arguments() and a main() function.
    modules = pkgutil.iter_modules(cli_package.__path__)
    for _, module_name, _ in modules:
        module = importlib.import_module("." + module_name, cli_package.__name__)
        help_message = module.__doc__.strip().split("\n", maxsplit=1)[0]
        subparser = subparsers.add_parser(
            module_name, help=help_message, description=module.__doc__
        )
        subparser.set_defaults(module=module)
        module.add_arguments(subparser)

    # Module 'run' needs to accept addition arguments
    args, extra_args = parser.parse_known_args(commandline_arguments)

    # For module 'run' extra_args are added to existing snakemake_args in namespace
    if hasattr(args, "snakemake_args"):
        args.snakemake_args += extra_args

    if args.debug:
        root = logging.getLogger()
        root.setLevel(logging.DEBUG)

    if not hasattr(args, "module"):
        parser.error("Please provide the name of a subcommand to run")
    else:
        module = args.module
        del args.module
        module_name = module.__name__.split('.')[-1]

        # Re-parse extra arguments if module is not "run" to raise the expected error
        if module_name != "run" and extra_args:
            parser.parse_args(extra_args)

        # Print settings for module
        sys.stderr.write(f"SETTINGS FOR: {module_name}\n")
        for object_variable, value in vars(args).items():
            sys.stderr.write(f" {object_variable}: {value}\n")

        module.main(args)

    return 0


if __name__ == "__main__":
    sys.exit(main())
