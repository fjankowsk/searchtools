#
#   2026 Fabian Jankowski
#   Fix multi-chunk PSRFITS search mode data.
#

import argparse
import os.path
import shutil
import sys

from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy import units as u


def parse_args():
    """
    Parse the commandline arguments.

    Returns
    -------
    args: populated namespace
        The commandline arguments.
    """

    parser = argparse.ArgumentParser(
        description="Fix multi-chunk PSRFITS search mode data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "files", type=str, nargs="+", help="Names of search mode data files to process."
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


#
# MAIN
#


def main():
    # handle command line arguments
    args = parse_args()
    check_args(args)

    print(args.files)

    SECPERDAY = 24 * 60 * 60

    for item in args.files:
        print(f"Processing: {item}")

        # check if required
        with fits.open(item, mode="readonly") as hdul:
            header1 = hdul[1].header
            if not header1["NSUBOFFS"] > 0:
                print(f"NSUBOFFS values is not positive. Skipping file: {item}")
                continue

        _root, _extension = os.path.splitext(item)
        outname = f"{_root}_fixed{_extension}"

        shutil.copy(item, outname)

        with fits.open(outname, mode="update") as hdul:
            header0 = hdul[0].header
            _mjd = header0["STT_IMJD"] + (
                header0["STT_SMJD"] + header0["STT_OFFS"]
            ) / float(SECPERDAY)
            start = Time(_mjd, format="mjd")
            print(f"Original start time: {start.mjd}, {start.iso}")

            header1 = hdul[1].header
            _offset = header1["NSUBOFFS"] * header1["NSBLK"] * header1["TBIN"]
            offset = TimeDelta(_offset * u.second)
            print(f"NSUBOFFS time offset: {offset}")

            correct_start = start + offset
            print(f"Corrected start time: {correct_start.mjd}, {correct_start.iso}")

            # update start time and reset nsuboffs
            mjd_integer = int((correct_start.jd1 - 2400000.5) + correct_start.jd2)
            mjd_fraction = (
                correct_start.jd1 - 2400000.5 + correct_start.jd2
            ) - mjd_integer
            smjd_integer = int(mjd_fraction * SECPERDAY)
            smjd_fraction = (mjd_fraction * SECPERDAY) - smjd_integer

            hdul[0].header["STT_IMJD"] = mjd_integer
            hdul[0].header["STT_SMJD"] = smjd_integer
            hdul[0].header["STT_OFFS"] = smjd_fraction

            hdul[1].header["NSUBOFFS"] = 0

    print("All done.")


if __name__ == "__main__":
    main()
