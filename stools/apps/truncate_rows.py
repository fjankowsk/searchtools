#
#   2026 Fabian Jankowski
#   Truncate PSRFITS search mode data.
#

import argparse
import os.path
import shutil
import sys

from astropy.io import fits


class CustomFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter
):
    pass


def parse_args():
    """
    Parse the commandline arguments.

    Returns
    -------
    args: populated namespace
        The commandline arguments.
    """

    _epilog = """\
This program truncates PSRFITS search mode data to a given number of sub:nrows of samples, specified on the command line.

This is useful for time-aligned data that differ in row count due to subtle offsets in backend process or thread synchronisation when a stop command is received. The program truncates all input files to the same number of rows so that they can be processed. A common task is to truncate all frequency lanes to the same nrows to splice them together.
"""

    parser = argparse.ArgumentParser(
        description="Truncate PSRFITS search mode data.",
        epilog=_epilog,
        formatter_class=CustomFormatter,
    )

    parser.add_argument(
        "files", type=str, nargs="+", help="Names of search mode data files to process."
    )

    parser.add_argument(
        "--nrows",
        type=int,
        required=True,
        metavar=("value"),
        help="The number of the data rows to retain in the output files.",
    )

    args = parser.parse_args()

    return args


def check_args(args):
    """
    Sanity check the commandline arguments.

    Parameters
    ----------
    args: populated namespace
        The commandline arguments.
    """

    # files
    for item in args.files:
        if not os.path.isfile(item):
            print(f"File does not exist: {item}")
            sys.exit(1)

    # nrows
    if not args.nrows >= 1:
        print(f"The number of rows is invalid: {args.nrows}")
        sys.exit(1)


#
# MAIN
#


def main():
    # handle command line arguments
    args = parse_args()
    check_args(args)

    print(args.files)
    print(f"Retaining nrows: {args.nrows}")

    for item in args.files:
        print(f"Processing: {item}")

        _root, _extension = os.path.splitext(item)
        outname = f"{_root}_truncated{_extension}"

        shutil.copy(item, outname)

        with fits.open(outname, mode="update") as hdul:
            data = hdul[1].data
            # truncate
            data = data[0 : args.nrows]
            hdul[1].data = data

    print("All done.")


if __name__ == "__main__":
    main()
