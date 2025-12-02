#
#   2025 Fabian Jankowski
#   Time align and frequency splice search mode data.
#

import argparse
import os.path
import sys

import numpy as np
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
        description="Time align and frequency splice search mode data.",
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

    # sanity check files
    for item in args.files:
        if not os.path.isfile(item):
            print(f"File does not exist: {item}")
            sys.exit(1)


def get_global_start_and_end_times(files):
    """
    Get the global start and end times from the files.

    Parameters
    ----------
    files: ~np.array
        The input data filenames.

    Returns
    -------
    gtstart, gtend: float
        The global start and end MJDs.
    """

    gtstart = None
    gtend = None

    for item in files:
        yobj = Your(item)
        print(yobj.tstart, yobj.tend)

        # tstart
        if gtstart is None:
            gtstart = yobj.tstart
        else:
            if yobj.tstart > gtstart:
                gtstart = yobj.tstart

        # tend
        if gtend is None:
            gtend = yobj.tend
        else:
            if yobj.tend < gtend:
                gtend = yobj.tend

    return gtstart, gtend


def fft_shift_2d(signal_2d, shift_samples, pad_factor=2):
    """
    Fine shift a time domain signal by an non-integer samples amount
    in the Fourier domain.

    The function is vectorized along the frequency axis (axis 1).

    Parameters
    ----------
    signal_2d: ~np.array
        The input time domain data of shape (nsamp, nfreq).
    shift_samples: float
        The shift in time samples.
    pad_factor: int
        Pad the FFT by this amount.

    Returns
    -------
    result: ~np.array
        The shifted signal of shape (nsamp, nfreq).
    """

    # time along axis 0
    # freq along axis 1
    assert len(signal_2d.shape) == 2

    N = signal_2d.shape[0]
    N_padded = N * pad_factor

    # compute FFT with zero padding
    fft_signal = np.fft.fft(signal_2d, n=N_padded, axis=0)

    # frequency axis for padded signal
    freq = np.fft.fftfreq(N_padded)

    # phase shift for non-integer delay
    phase_shift = np.exp(-1j * 2 * np.pi * freq * shift_samples)

    # apply shift and inverse FFT
    shifted_fft = fft_signal * phase_shift[:, None]
    shifted_signal = np.fft.ifft(shifted_fft, axis=0)

    # return real part and trim back to original length
    result = np.real(shifted_signal[:N, :])

    # cast back to input dtype
    result = result.astype(signal_2d.dtype)

    assert result.shape == signal_2d.shape
    assert result.dtype == signal_2d.dtype

    return result


#
# MAIN
#


def main():
    # handle command line arguments
    args = parse_args()
    check_args(args)

    files = args.files
    files = np.array(files)
    print(files)

    nchans = []
    fch1s = []

    for item in files:
        yobj = Your(item)
        assert yobj.foff < 0
        nchans.append(yobj.nchans)
        fch1s.append(yobj.fch1)

    gnchan = int(np.sum(nchans))
    print(f"Global number of channels: {gnchan}")

    # sort by decreasing frequency
    idx_sort = np.argsort(fch1s)
    idx_sort = idx_sort[::-1]
    files = files[idx_sort]

    print(files)

    # global start and end times
    gtstart, gtend = get_global_start_and_end_times(files)

    print(f"Global start and end times: {gtstart}, {gtend}")

    # global number of samples
    gnsamp = (gtend - gtstart) * 60 * 60 * 24.0 / yobj.tsamp
    gnsamp = int(np.floor(gnsamp))

    print(f"Total number of samples: {gnsamp}")

    gulp = 16 * 1024

    # prepare writer
    yobj = Your(files[0])

    outname = "merged.fil"
    if os.path.isfile(outname):
        print(f"Output file exists: {outname}")
        sys.exit(1)

    sigproc_object = make_sigproc_object(
        rawdatafile=outname,
        source_name=yobj.source_name.decode("utf-8"),
        nchans=gnchan,
        foff=yobj.foff,  # MHz
        fch1=yobj.fch1,  # MHz
        tsamp=yobj.tsamp,  # seconds
        tstart=gtstart,  # MJD
        src_raj=yobj.src_raj,  # HHMMSS.SS
        src_dej=yobj.src_dej,  # DDMMSS.SS
        machine_id=yobj.machine_id,
        nbeams=yobj.nbeams,
        ibeam=yobj.ibeam,
        nbits=yobj.nbits,
        nifs=1,
        barycentric=0,
        pulsarcentric=0,
        telescope_id=yobj.telescope_id,
        data_type=yobj.data_type,
        az_start=yobj.az_start,
        za_start=yobj.za_start,
    )

    sigproc_object.write_header(outname)

    pbar = tqdm(total=gnsamp)
    samples_left = gnsamp
    sstart = 0

    while samples_left > 0:
        if samples_left < gulp:
            gulp = samples_left

        total = None

        for item in files:
            yobj = Your(item)

            # starting sample offset
            soffset = (gtstart - yobj.tstart) * 60 * 60 * 24.0 / yobj.tsamp
            soffset_int = int(np.floor(soffset))
            soffset_frac = soffset - soffset_int
            assert soffset >= 0
            assert soffset_int >= 0
            assert soffset_frac >= 0

            # print(soffset, soffset_int, soffset_frac)

            # get the signal offset by the integer sample amount
            _data = yobj.get_data(sstart + soffset_int, gulp, pol=0, npoln=1)
            assert _data.shape[0] == gulp

            # shift the data by the fractional part
            if soffset_frac > 0:
                _data = fft_shift_2d(_data, soffset_frac)
                assert _data.shape[0] == gulp

            if total is None:
                total = _data.copy()
            else:
                total = np.hstack([total, _data])

        sigproc_object.append_spectra(total, outname)
        sstart += gulp
        samples_left -= gulp
        pbar.update(sstart - pbar.n)

    pbar.update(0)
    pbar.close()

    print("All done.")


if __name__ == "__main__":
    main()
