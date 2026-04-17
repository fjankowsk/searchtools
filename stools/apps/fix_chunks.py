import os.path
import shutil

from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy import units as u

#
# MAIN
#


def main():
    SECPERDAY = 24 * 60 * 60

    files = ["merged.fits"]

    for item in files:
        print(f"Processing: {item}")

        _root, _extension = os.path.splitext(item)
        outname = f"{_root}_fixed{_extension}"

        shutil.copy(item, outname)

        with fits.open(outname, mode="update") as hdul:
            header0 = hdul[0].header
            _mjd = header0["STT_IMJD"] + (
                header0["STT_SMJD"] + header0["STT_OFFS"]
            ) / float(SECPERDAY)
            start = Time(_mjd, format="mjd")
            print(start.mjd)
            print(start.iso)

            header1 = hdul[1].header
            _offset = header1["NSUBOFFS"] * header1["NSBLK"] * header1["TBIN"]
            offset = TimeDelta(_offset * u.second)
            print(offset)

            correct_start = start + offset
            print(correct_start.mjd)
            print(correct_start.iso)

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
