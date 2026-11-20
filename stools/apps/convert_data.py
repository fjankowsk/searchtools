#
#   2023 - 2026 Fabian Jankowski
#   Convert and clean search mode data.
#

import argparse

from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from mtcutils.core import automatic_channel_mask_fj, normalise, zdot
import numpy as np
import pandas as pd
from tqdm import tqdm
from your import Your
from your.formats.filwriter import make_sigproc_object

from stools.cc2d import cross_correlate_2d


def parse_args():
    """
    Parse the commandline arguments.

    Returns
    -------
    args: populated namespace
        The commandline arguments.
    """

    parser = argparse.ArgumentParser(
        description="Convert PSRFITS to SIGPROC filterbank data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("file", type=str, help="Filename of PSRFITS file to process.")

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        dest="outfile",
        default="converted.fil",
        metavar=("filename"),
        help="Name of the output SIGPROC filerbank file.",
    )

    parser.add_argument(
        "--plot",
        dest="plot",
        action="store_true",
        default=False,
        help="Display the data at each step.",
    )

    parser.add_argument(
        "--rfi",
        dest="rfi",
        action="store_true",
        default=False,
        help="Perform RFI excision.",
    )

    parser.add_argument(
        "--rmedian",
        dest="rmedian",
        action="store_true",
        default=False,
        help="Remove a running median from the data in each channel to flatten the baseline and bandpass.",
    )

    parser.add_argument(
        "--2dcorr",
        dest="twodcorr",
        action="store_true",
        default=False,
        help="Remove the 6 minute analogue pointing jumps by cross-correlating with a 2D template.",
    )

    parser.add_argument(
        "--tscrunch",
        type=int,
        dest="tscrunch",
        default=None,
        help="Decimate the data in time by this factor.",
    )

    args = parser.parse_args()

    return args


def running_median(t_data, params):
    """
    Calculate a running median across the data.

    Parameters
    ----------
    t_data: ~np.array
        The data to be analysed.
    params: dict
        Additional parameters for the processing.

    Returns
    -------
    rmedian: ~np.array
        The running median.
    """

    data = np.copy(t_data)

    if "window" in params:
        window = params["window"]
    else:
        window = 20

    if "win_type" in params:
        win_type = params["win_type"]
    else:
        win_type = None

    df = pd.DataFrame(data, columns=["data"])
    rmedian = (
        df["data"]
        .rolling(center=True, min_periods=1, window=window, win_type=win_type)
        .median()
        .to_numpy()
    )

    return rmedian


def plot_data_1d(data, title):
    """
    Display the 1D data.
    """

    fig = plt.figure()
    ax = fig.add_subplot()

    ax.plot(data, color="black")

    ax.set_title(title)

    fig.tight_layout()


def plot_data_2d(data, title, showmax=False):
    """
    Display the 2D data.
    """

    fig = plt.figure()
    ax = fig.add_subplot()

    im = ax.imshow(data, aspect="auto", cmap="Greys", origin="upper")

    if showmax:
        argmax = np.argmax(data)
        argmax = np.unravel_index(argmax, data.shape)
        ax.scatter(argmax[1], argmax[0], marker="*", color="tab:red", zorder=7)

    cbar = plt.colorbar(im, ax=ax, location="right")

    ax.set_xlabel("Time sample")
    ax.set_ylabel("Frequency channel")
    ax.set_title(title)

    fig.tight_layout()


def decimate_2d(data, f0, f1):
    """
    Decimate the 2D data in its dimensions.

    Parameters
    ----------
    data: ~np.array
        The input data.
    f0: int
        The factor by which to decimate dimension 0.
    f1: int
        The factor by which to decimate dimension 1.

    Returns
    -------
    decimated: ~np.array
        The decimated data.

    Raises
    ------
    ValueError
        If the data dimensions cannot be exactly divided by f0 or f1.
    """

    shape = data.shape

    for axis, item in zip([0, 1], [f0, f1]):
        if shape[axis] % item > 0:
            raise ValueError(
                f"Scrunch factor must divide data shape exactly: {shape[axis]}, {item}"
            )

    decimated = data.reshape(shape[0] // f0, f0, shape[1] // f1, f1).mean(axis=(1, 3))

    return decimated


def convert_to_integer(data, nbit=8):
    """
    Convert the float data to unsigned integers for storing in SIGPROC
    filterbank files or PSRFITS ones.

    Note that SIGPROC does not implement scales and offset storage.

    Parameters
    ----------
    data: ~np.array of float
        The input data.
    nbit: int
        The number of bits in the output data.

    Returns
    -------
    data_int: ~np.array of int
        The output data in integer form.

    Raises
    ------
    NotImplementedError
        If the requested bit depth is not supported.
    """

    mapping = {8: np.uint8, 16: np.uint16, 32: np.uint32}

    if nbit not in mapping:
        raise NotImplementedError(
            "Requested bit depth not implemented: {0}".format(nbit)
        )

    # normalise the data before integer conversion
    # we follow the psrfits convention here
    # https://www.atnf.csiro.au/research/pulsar/psrfits_definition/PsrfitsDocumentation.html
    # real_value = (data_value - zero_offset) * data_scale + data_offset
    # => data_value = (real_value - data_offset) / data_scale + zero_offset
    zero_offset = 2 ** (nbit - 1) - 0.5
    data_offset = np.mean(data)
    # add some extra dynamic range
    data_scale = (np.max(data) - data_offset) / float(zero_offset - 0.5)

    print(
        "Zero offset, data offset, data scale: {0}, {1}, {2}".format(
            zero_offset, data_offset, data_scale
        )
    )
    if data_scale > 0:
        scaled_data = (data - data_offset) / data_scale + zero_offset
    else:
        scaled_data = data

    print(
        "Scaled data: {0}, {1}, {2}".format(
            np.min(scaled_data), np.mean(scaled_data), np.max(scaled_data)
        )
    )

    assert np.all(np.isfinite(scaled_data))

    # convert to unsigned integer
    int_type = mapping[nbit]
    min_int = 0
    max_int = 2**nbit - 1
    data_int = np.clip(np.rint(scaled_data), min_int, max_int).astype(int_type)

    assert np.all(np.isfinite(data_int))

    print(
        "Integer data: {0}, {1}, {2}".format(
            np.min(data_int), np.mean(data_int), np.max(data_int)
        )
    )

    return data_int


# use the inverse dispersion constant rounded to three digits that
# is in common use in pulsar astronomy to be consistent throughout
# software tools (see page 129 of Manchester and Taylor 1977)
KDM = 1.0 / 2.41e-4


def get_dm_smearing(f_lo, f_hi, dm):
    """
    Compute the intra-channel dispersive smearing.

    This returns the *full* DM smearing width from -dm_sweep to +dm_sweep,
    and not just a single side. In other words, one must use the range
    -0.5 * get_dm_smearing to +0.5 * get_dm_smearing.

    Parameters
    ----------
    f_lo: float
        The low frequency edge of the channel in MHz.
    f_hi: float
        The high frequency edge of the channel in MHz.
    dm: float
        The dispersion measure of the source in pc/cm^3.

    Returns
    -------
    dt: float
        The intra-channel dispersive smearing in seconds.
    """

    dt = KDM * (f_lo ** (-2) - f_hi ** (-2)) * dm

    return dt


def get_impulse_template(robj, fch1, foff, shape, idx_middle):
    """
    Construct a 2D impulse template in FT order.

    Parameters
    ----------
    robj: ~your.Your
        The input object.
    fch1: float
        The frequency of the first (highest) channel.
    foff: float
        The channel bandwidth foff < 0.
    shape: tuple
        The shape of the output array.
    idx_middle: int
        The location index of the impulse. This is the fiducial point.

    Returns
    -------
    impulse: ~np.array
        The impulse template in FT order.
    """

    nchan = shape[0]
    nsamp = shape[1]

    assert foff < 0
    assert 0 <= idx_middle <= nsamp - 1

    # coherent dedispersion dm
    coherent_dm = robj.header["CHAN_DM"]

    # assemble the impulse template
    # boxcar width of the impulse plus twice the dm smearing time
    f_hi = fch1
    f_lo = fch1 + foff
    impulse_time_max = 0.2 + get_dm_smearing(f_lo, f_hi, coherent_dm)
    try:
        impulse_time_max = impulse_time_max[0]
    except IndexError:
        print("Impulse time max is float already.")

    print(f"Maximum impulse time: {impulse_time_max:.2f} s")

    impulse_nsamp = int(np.ceil(impulse_time_max / robj.your_header.tsamp))
    assert impulse_nsamp <= nsamp

    impulse = np.zeros(shape=shape)

    # boxcar impulse
    # the boxcar impulse was ~40 ms wide in june 2023
    # and about 200 ms in may 2024
    width_boxcar = 40.0e-3  # in seconds
    nsamp_boxcar = int(np.ceil(width_boxcar / robj.your_header.tsamp))
    nsamp_boxcar_half = int(np.ceil(0.5 * nsamp_boxcar))
    idx_start = idx_middle - nsamp_boxcar_half
    idx_end = idx_middle + nsamp_boxcar_half
    impulse[:, idx_start:idx_end] = 1.0

    # dm sweep
    fchans = fch1 + foff * np.arange(nchan)
    for ichan, f_hi in fchans:
        f_lo = f_hi + foff
        sweep_end = 0.5 * get_dm_smearing(f_lo, f_hi, coherent_dm)
        delta_sweep = int(np.ceil(sweep_end / robj.your_header.tsamp))

        # right side
        idx_end = idx_middle + nsamp_boxcar_half + delta_sweep
        if idx_end > nsamp - 1:
            idx_end = nsamp - 1
        impulse[ichan, idx_middle:idx_end] = 1.0

        # left side
        idx_start = idx_middle - nsamp_boxcar_half - delta_sweep
        if idx_start < 0:
            idx_start = 0
        impulse[ichan, idx_start:idx_middle] = 1.0

    return impulse


def main():
    args = parse_args()
    plt.ion()

    nbit = 8

    robj = Your(args.file)
    print(robj.your_header)
    print(f"File format: {robj.format}")

    if robj.your_header.ra_deg is not None and robj.your_header.dec_deg is not None:
        position = SkyCoord(
            robj.your_header.ra_deg, robj.your_header.dec_deg, unit="deg"
        )
        ra_hms = position.ra.hms
        dec_dms = position.dec.dms
        src_raj = float(
            f"{int(ra_hms[0]):02d}{abs(int(ra_hms[1])):02d}{abs(ra_hms[2]):07.4f}"
        )
        src_dej = float(
            f"{int(dec_dms[0]):02d}{abs(int(dec_dms[1])):02d}{abs(dec_dms[2]):07.4f}"
        )
    else:
        src_raj = 0.0
        src_dej = 0.0

    # fch1 must be the highest frequency channel in sigproc files
    # flip accordingly
    if robj.foff < 0:
        fch1 = robj.your_header.fch1
        foff = robj.your_header.foff
    else:
        fch1 = robj.freqs[-1]
        foff = -1.0 * abs(robj.your_header.foff)

    assert foff < 0
    assert fch1 == np.max(robj.chan_freqs)

    if robj.format == "fil":
        telescope_id = robj.telescope_id
        machine_id = robj.machine_id
    else:
        if "nenufar" in robj.telescope.lower():
            telescope_id = 40  # nenufar
            machine_id = 40  # nenufar-luppi
        else:
            telescope_id = 0  # fake
            machine_id = 0  # fake

    if args.tscrunch is not None:
        tsamp = robj.your_header.tsamp * args.tscrunch
    else:
        tsamp = robj.your_header.tsamp

    sigproc_obj = make_sigproc_object(
        rawdatafile=args.outfile,
        source_name=robj.source_name,
        # this is the *centre* frequency of the first filterbank channel
        fch1=fch1,  # MHz
        foff=foff,  # MHz
        nchans=robj.your_header.nchans,
        tsamp=tsamp,  # seconds
        tstart=robj.your_header.tstart,  # MJD
        src_raj=src_raj,  # HHMMSS.SS
        src_dej=src_dej,  # DDMMSS.SS
        machine_id=machine_id,  # backend
        nbeams=1,
        ibeam=0,
        nbits=nbit,
        nifs=1,  # total intensity
        barycentric=0,
        pulsarcentric=0,
        telescope_id=telescope_id,  # telescope
        data_type=1,  # filterbank
        az_start=-1,
        za_start=-1,
    )

    sigproc_obj.write_header(args.outfile)

    # data processing gulp size
    gulp_time = 50.0  # in seconds
    gulp = int(np.ceil(gulp_time / robj.your_header.tsamp))
    # use an even number for convenience
    if gulp % 2 > 0:
        gulp -= 1
    assert gulp % 2 == 0

    if args.tscrunch is not None:
        while gulp % args.tscrunch > 0:
            gulp -= 1

        assert gulp % args.tscrunch == 0
    print(
        "Gulp size: {0} samples, {1:.1f} seconds".format(
            gulp, gulp * robj.your_header.tsamp
        )
    )
    istart = 0
    nchan = robj.your_header.nchans
    nspectra = robj.nspectra()
    pbar = tqdm(total=nspectra)

    # running median window size for baseline flattening
    if args.rmedian:
        rmedian_window = 20.0  # in seconds
        mparams = {"window": int(np.ceil(rmedian_window / robj.your_header.tsamp))}
        print(
            "Running median window size: {0} samples, {1:.1f} seconds".format(
                mparams["window"], mparams["window"] * robj.your_header.tsamp
            )
        )

    # coherent dedispersion dm
    coherent_dm = robj.header["CHAN_DM"]
    print(f"Data were coherently dedispersed at the DM: {coherent_dm} pc/cm3")

    # assemble the impulse template for the 6-min analogue re-pointing jumps
    impulse = get_impulse_template(
        robj, fch1, foff, shape=(nchan, gulp), idx_middle=gulp // 2
    )
    impulse -= np.mean(impulse)
    impulse /= np.linalg.norm(impulse)

    while istart < nspectra:
        # always produce total intensity data
        data = robj.get_data(nstart=istart, nsamp=gulp, pol=0, npoln=1)
        nsamp = data.shape[0]

        # data in TF order
        assert data.shape[0] == nsamp
        assert data.shape[1] == nchan

        # transpose to FT order
        data = data.T
        assert data.shape[0] == nchan
        assert data.shape[1] == nsamp

        # flip frequency band
        if robj.foff > 0:
            data = np.flip(data, axis=0)

        if args.plot:
            plot_data_2d(data, "loaded")

        # normalise
        data, _, _ = normalise(data)

        if args.plot:
            plot_data_2d(data, "normalised")

        if args.twodcorr:
            # find the 6-min analogue re-pointing jumps by 2d correlating
            # with a matched boxcar filter in the frequency - time domain
            _chunk = data - np.mean(data)
            _chunk /= np.linalg.norm(_chunk)
            corr = cross_correlate_2d(_chunk, impulse, mode="same")

            print("Maximum cross-correlation: {0:.4f}".format(np.max(corr)))
            argmax = np.argmax(corr)
            argmax = np.unravel_index(argmax, corr.shape)
            argmax_sample = argmax[1]
            print(f"Time sample of the maximum cross-correlation: {argmax_sample}")

            if args.plot:
                plot_data_2d(data, "data")
                plot_data_2d(impulse, "template")
                plot_data_2d(corr, "correlation", showmax=True)

            if np.max(corr) > 0.1:
                # remove the sweep
                zap_array = get_impulse_template(
                    fch1, foff, shape=(nchan, nsamp), idx_middle=argmax_sample
                )
                zap_array += 1.0
                _mask = np.isclose(zap_array, 2.0)
                zap_array[_mask] = 0.0

                data = data * zap_array

                if args.plot:
                    plot_data_2d(zap_array, "zap-array")
                    plot_data_2d(data, "post-zap")

        if args.rfi:
            # iqrm & acc1 masks
            mask = automatic_channel_mask_fj(data)
            print(
                "Channels masked based on stddev (via IQRM) and acc1: {} / {} ({:.2%})".format(
                    mask.sum(), nchan, mask.sum() / nchan
                )
            )
            data[mask] = 0.0

            if args.plot:
                plot_data_2d(data, "masked")

            # zdot
            data = zdot(data)

            if args.plot:
                plot_data_2d(data, "zdot")

        if args.rmedian:
            # remove a running median per channel from the data
            # this flattens the baseline and bandpass further
            for ichan in range(nchan):
                data[ichan, :] -= running_median(data[ichan, :], mparams)

            if args.plot:
                plot_data_2d(data, "rmedian")

        # normalise again after rfi excision
        data, _, _ = normalise(data)

        # decimate the FT data in time
        if args.tscrunch is not None:
            if nsamp % args.tscrunch > 0:
                print(f"Last gulp detected of sample size: {nsamp}")
                last_gulp = nsamp
                while last_gulp % args.tscrunch > 0:
                    last_gulp -= 1

                assert last_gulp % args.tscrunch == 0
                print(f"Last gulp truncated to sample size: {last_gulp}")
                data = data[:, 0 : last_gulp - 1]

            data = decimate_2d(data, f0=1, f1=args.tscrunch)

        # convert to integers
        data_int = convert_to_integer(data, nbit=nbit)

        if args.plot:
            plot_data_2d(data_int, "integer")
            input("Press any key.")

        # the function expects TF ordered data
        data_int = data_int.T
        sigproc_obj.append_spectra(data_int, args.outfile)

        istart += gulp
        pbar.update(istart - pbar.n)

    pbar.update(nspectra - pbar.n)
    pbar.close()

    your_obj = Your(args.outfile)
    print(your_obj.your_header)

    assert np.isclose(robj.your_header.ra_deg, your_obj.your_header.ra_deg)
    assert np.isclose(robj.your_header.dec_deg, your_obj.your_header.dec_deg)

    print("All done.")


if __name__ == "__main__":
    main()
