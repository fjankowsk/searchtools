#
#   2023 - 2026 Fabian Jankowski
#   Convert and clean search mode data.
#

import argparse
import sys

from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from mtcutils.core import automatic_channel_mask_fj, normalise, zdot
import numpy as np
import pandas as pd
from scipy.fft import rfft, irfft, next_fast_len
from tqdm import tqdm
from your import Your
from your.formats.filwriter import make_sigproc_object


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
        "-m",
        "--mask",
        action="append",
        type=float,
        dest="mask",
        nargs=2,
        metavar=("start", "end"),
        default=None,
        help="Mask the static frequency range given in MHz from start to end inclusive. Can use several times to specify multiple frequency ranges.",
    )

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
        "--tscrunch",
        type=int,
        dest="tscrunch",
        default=None,
        help="Decimate the data in time by this factor.",
    )

    parser.add_argument(
        "--2dcorr",
        dest="twodcorr",
        action="store_true",
        default=False,
        help="Remove the NenuFAR 6 minute analogue re-pointing jumps by cross-correlating with a 2D impulse template.",
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

    # sanity check freqeuncy mask
    if args.mask is not None:
        for item in args.mask:
            if item[0] > 0 and item[1] > 0 and item[0] < item[1]:
                pass
            else:
                print(f"Mask frequencies are invalid: {item[0]}, {item[1]}")
                sys.exit(1)


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


def convert_to_integer(data, nbit=8, data_offset=None, data_scale=None):
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
    data_offset: float or None
        The data offset to use. If None, determine from the data.
    data_scale: float or None
        The data scale to use. If None, determine from the data.

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

    if data_offset is None:
        data_offset = np.mean(data)

    if data_scale is None:
        # add some extra dynamic range
        data_scale = (np.max(data) - data_offset) / (zero_offset - 0.5)

    print(
        f"Zero offset, data offset, data scale: {zero_offset}, {data_offset}, {data_scale}"
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


def gaussian_window(shape, sigma):
    """
    2D Gaussian window.
    """

    h, w = shape
    y = np.arange(h) - h // 2
    x = np.arange(w) - w // 2
    Y, X = np.meshgrid(y, x, indexing="ij")

    return np.exp(-(X**2 + Y**2) / (2 * sigma**2))


def get_impulse_template(robj, fch1, foff, shape, idx_middle, window):
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
    window: bool
        Whether to apply a Gaussian window to the template.

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

    fchans = fch1 + foff * np.arange(nchan)

    # assemble the impulse template
    # which is just the intrachannel dm smearing
    impulse = np.zeros(shape=shape)

    # this is the dispersive sweep in the lowest frequency channel
    f_hi = np.min(fchans)
    f_lo = f_hi + foff
    impulse_time_max = get_dm_smearing(f_lo, f_hi, coherent_dm)
    print(f"Maximum impulse time: {impulse_time_max:.2f} s")

    impulse_nsamp = int(np.ceil(impulse_time_max / robj.your_header.tsamp))
    assert impulse_nsamp <= nsamp

    # dm sweep
    for ichan, f_hi in enumerate(fchans):
        f_lo = f_hi + foff
        sweep_end = 0.5 * get_dm_smearing(f_lo, f_hi, coherent_dm)
        delta_sweep = int(np.ceil(sweep_end / robj.your_header.tsamp))

        # left side
        idx_start = idx_middle - delta_sweep
        if idx_start < 0:
            idx_start = 0

        # right side
        idx_end = idx_middle + delta_sweep
        if idx_end > nsamp - 1:
            idx_end = nsamp - 1

        assert idx_start < idx_end
        impulse[ichan, idx_start:idx_end] = 1.0

    if window:
        sigma = int(np.ceil(0.5 * impulse_time_max / robj.your_header.tsamp))
        impulse *= gaussian_window(shape, sigma=sigma)

    return impulse


def time_shift_detection(data, template):
    """
    Per-channel 1D time correlation, summed over frequency.
    Handles edge clipping via overlap normalization.

    Parameters
    ----------
    data: (T, F)
        The real input data array in TF order.
    template: (T_t, F)
        The template real array in TF order. It may be zero-padded and will be auto-trimmed.

    Returns
    -------
    R_norm: (T,)
        The normalized cross correlation.
    """

    T = data.shape[0]
    T_t_full = template.shape[0]
    assert data.shape[1] == template.shape[1]

    # auto-trim template to non-zero region (with 1-sample margin)
    row_energy = np.sum(template**2, axis=1)
    nonzero = np.where(row_energy > 0)[0]
    if len(nonzero) == 0:
        raise ValueError("Template is entirely zero.")

    t0, t1 = max(nonzero[0] - 1, 0), min(nonzero[-1] + 2, T_t_full)
    template = template[t0:t1]
    T_t = template.shape[0]

    # select good fft size
    fft_len = next_fast_len(T + T_t - 1)

    # batched rfft
    # n=fft_len zero-pads automatically
    D = rfft(data, n=fft_len, axis=0)  # (fft_len//2+1, F)
    P = rfft(template, n=fft_len, axis=0)  # (fft_len//2+1, F)

    # correlation via conjugate product and irfft
    R_full = np.real(irfft(D * P.conj(), n=fft_len, axis=0)).sum(axis=1)

    # crop to "same" size (template centre at each data sample)
    center = T_t // 2
    R = R_full[center : center + T]

    # overlap normalization to handle edge clipping
    row_e = np.sum(template**2, axis=1)
    cum = np.concatenate([[0.0], np.cumsum(row_e)])
    i = np.arange(T)
    ls = np.maximum(0, center - i)
    le = np.minimum(T_t, T - i + center)
    overlap_energy = cum[le] - cum[ls]

    R_norm = R / np.where(overlap_energy > 0, overlap_energy, 1.0)

    return R_norm


def main():
    args = parse_args()
    check_args(args)
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
    if args.twodcorr:
        impulse = get_impulse_template(
            robj, fch1, foff, shape=(nchan, gulp), idx_middle=gulp // 2, window=False
        )

    # build static frequency zap mask
    if args.mask is not None:
        print(f"Using static frequency zap mask: {args.mask}")
        fchans = fch1 + foff * np.arange(nchan)

        # do not zap by default
        freqzapmask = np.zeros(nchan, dtype=bool)

        for item in args.mask:
            _mask = (fchans >= item[0]) & (fchans <= item[1])
            freqzapmask[_mask] = True

        print(f"Static frequency zap mask: {freqzapmask}")

    # determine these from the first gulp
    data_offset = None
    data_scale = None

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

        # apply static frequency zap mask
        if args.mask is not None:
            data[freqzapmask, :] = 0.0

        if args.plot:
            plot_data_2d(data, "frequency masked")

        # normalise
        data, _, _ = normalise(data)

        if args.plot:
            plot_data_2d(data, "normalised")

        if args.twodcorr:
            # find the 6-min analogue re-pointing jumps by 2d correlating
            # with a matched boxcar filter in the frequency - time domain
            _chunk = data - np.mean(data, axis=0)
            stdev = np.std(_chunk, axis=0)
            mask = stdev > 0
            _chunk[mask, :] /= stdev[mask]

            normcorr = time_shift_detection(_chunk, impulse)
            print("Maximum cross-correlation: {0:.4f}".format(np.max(normcorr)))
            idx_middle = np.argmax(normcorr)

            if args.plot:
                plot_data_2d(data, "data")
                plot_data_2d(impulse, "template")
                plot_data_1d(normcorr, "correlation")

            if np.max(normcorr) > 0.1:
                # remove the time samples corresponding to the maximum sweep
                # in the lowest frequency channel
                # we remove the entire box, because neighboring samples are affected
                # even in the high-frequency part of the band
                f_hi = fch1 + (nchan - 1) * foff
                f_lo = f_hi + foff
                sweep_end = 0.5 * get_dm_smearing(f_lo, f_hi, coherent_dm)
                delta_sweep = int(np.ceil(sweep_end / robj.your_header.tsamp))

                # left side
                idx_start = idx_middle - delta_sweep
                if idx_start < 0:
                    idx_start = 0

                # right side
                idx_end = idx_middle + delta_sweep
                if idx_end > nsamp - 1:
                    idx_end = nsamp - 1

                assert idx_start < idx_end
                data[:, idx_start:idx_end] = 0.0

                if args.plot:
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

        # determine data offset and scale from the first gulp
        if data_offset is None:
            data_offset = np.median(data)
            data_max = np.quantile(data, 0.99)
            zero_offset = 2 ** (nbit - 1) - 0.5
            data_scale = (data_max - data_offset) / (zero_offset - 0.5)

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
        data_int = convert_to_integer(
            data, nbit=nbit, data_offset=data_offset, data_scale=data_scale
        )

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
