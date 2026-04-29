#
#   2026 Fabian Jankowski
#   Fix multi-chunk PSRFITS search mode data.
#

import argparse
from decimal import Decimal
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

    SECPERDAY = Decimal("86400")

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
            _mjd = (
                Decimal(header0["STT_IMJD"])
                + (Decimal(header0["STT_SMJD"]) + Decimal(header0["STT_OFFS"]))
                / SECPERDAY
            )
            start = Time(_mjd, format="mjd")
            print(f"Original start time: {start.mjd}, {start.iso}")

            header1 = hdul[1].header
            _offset_subints = Decimal(header1["NSUBOFFS"])
            _offset_samples = _offset_subints * Decimal(header1["NSBLK"])
            _offset_time = _offset_samples * Decimal(header1["TBIN"])
            offset = TimeDelta(_offset_time * u.second)
            print(
                f"NSUBOFFS offset: {_offset_subints} subints, {_offset_samples} samples, {offset.sec} s"
            )

            correct_start = start + offset
            print(f"Corrected start time: {correct_start.mjd}, {correct_start.iso}")

            # update start time and reset nsuboffs
            hdul[0].header["DATE-OBS"] = correct_start.isot

            _correct_mjd = correct_start.to_value("mjd", subfmt="decimal")

            mjd_integer = int(_correct_mjd)
            mjd_fraction = _correct_mjd - mjd_integer
            smjd_integer = int(mjd_fraction * SECPERDAY)
            smjd_fraction = mjd_fraction * SECPERDAY - smjd_integer

            hdul[0].header["STT_IMJD"] = mjd_integer
            hdul[0].header["STT_SMJD"] = smjd_integer
            # XXX: this reduces the precision to double, which is ok
            # ideally, we want to retain the full decimal precision
            # we need to use the fits free card for this
            # unfortunately, it does not work at the moment
            hdul[0].header["STT_OFFS"] = float(smjd_fraction)

            hdul[1].header["NSUBOFFS"] = 0

    print("All done.")


if __name__ == "__main__":
    main()
