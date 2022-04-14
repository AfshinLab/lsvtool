import sys
import numpy as np
from collections import Counter
from pathlib import Path


class Summary(Counter):

    def print_stats(self, name=None, value_width=15, print_to=sys.stderr):
        """
        Prints stats in nice table with two column for the key and value pairs in summary
        :param name: name of script for header e.g. '__name__'
        :param value_width: width for values column in table
        :param print_to: Where to direct output. Default: stderr
        """
        # Get widths for formatting
        max_name_width = max(map(len, self.keys()), default=10)
        width = value_width + max_name_width + 1

        # Header
        print("="*width, file=print_to)
        print(f"STATS SUMMARY - {name}", file=print_to)
        print("-"*width, file=print_to)

        # Print stats in columns
        for name, value in self.items():
            value_str = str(value)
            if isinstance(value, (int, np.integer)):
                value_str = f"{value:>{value_width},}"
            elif isinstance(value, (float, np.float)):
                value_str = f"{value:>{value_width+4},.3f}"

            print(f"{name:<{max_name_width}} {value_str}", file=print_to)
        print("="*width, file=print_to)


def check_path(path_str):
    if path_str == "" or path_str is None:
        return None

    path = Path(path_str)
    if not path.exists():
        raise FileNotFoundError(f"Path {path_str} does not exist!")

    if not path.is_absolute():
        return path.resolve()
    return path
