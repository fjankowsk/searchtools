#
#   2023 - 2026 Fabian Jankowski
#   Plot a DM trial - time diagram from PRESTO's prepsubband output.
#

import argparse
import glob
import os.path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from tqdm import tqdm
from presto import infodata


def parse_args():
    """
    Parse the commandline arguments.

    Returns
    -------
    args: populated namespace
        The commandline arguments.
    """

    parser = argparse.ArgumentParser(
        description="Plot a DM trial - time diagram from PRESTO's prepsubband output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        action="store_true",
        default=False,
        help="Output plots to file rather than to screen.",
    )

    parser.add_argument(
        "--tscrunch",
        type=int,
        dest="tscrunch",
        default=None,
        help="Decimate the data in time by this factor.",
    )

    parser.add_argument(
        "--window",
        type=float,
        dest="window",
        default=30.0,
        help="Time window of the running median in seconds.",
    )

    args = parser.parse_args()

    return args


def get_iqr(data, axis):
    """
    Compute the interquartile range.

    Parameters
    ----------
    data: np.array
        The input data as numpy float array.

    Returns
    -------
    iqr: float
        The interquartile range.
    """

    q25, q75 = np.nanpercentile(data, [25.0, 75.0], axis=axis)

    iqr = q75 - q25

    return iqr


def get_robust_std(data, axis, mode="Gaussian"):
    """
    Compute the robust standard deviation using the IQR.

    Parameters
    ----------
    data: np.array
        The input data.
    mode: str
        The underlying distribution to assume when computing the robust standard deviation.

    Returns
    -------
    sr: float
        The robust standard deviation.
    """

    iqr = get_iqr(data, axis=axis)

    if mode == "Gaussian":
        fact = 0.7413
    elif mode == "Rice":
        fact = 0.9183
    elif mode == "Exponential":
        fact = 0.9102
    else:
        raise RuntimeError("Distribution not known: {0}".format(mode))

    sr = fact * iqr

    return sr


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


def main():
    args = parse_args()

    # exclude the zerodm files
    all_files = glob.glob("prepsub_DM*.dat")
    zerodm_files = glob.glob("prepsub_DM0.*.dat")

    files = [item for item in all_files if item not in zerodm_files]
    files = sorted(files)
    print("Number of files: {0}".format(len(files)))

    _data = np.fromfile(files[0], dtype=np.float32)
    print("Number of time samples: {0}".format(_data.shape[0]))

    if args.tscrunch is not None:
        assert _data.shape[0] % args.tscrunch == 0

    total = np.zeros(shape=(len(files), _data.shape[0]), dtype=np.float32)
    print("Total array shape: {0}".format(total.shape))

    # info file
    info_file = os.path.splitext(files[0])[0] + ".inf"
    info = infodata.infodata(info_file)
    tsamp = info.dt
    print("Sampling time: {0:.2f} ms".format(tsamp * 1e3))

    # running median window length
    print(f"Running median window length: {args.window} s")
    rparams = {"window": int(np.ceil(args.window / tsamp))}

    dms = []

    for i, file in enumerate(tqdm(files)):
        _temp = np.fromfile(file, dtype=np.float32)

        # subtract a running median from each time series
        _temp -= running_median(_temp, params=rparams)

        total[i, :] = _temp

        # save dm
        info_file = os.path.splitext(file)[0] + ".inf"
        info = infodata.infodata(info_file)
        dms.append(info.DM)

    # sort by increasing dm
    dms = np.array(dms)
    _idxs = np.argsort(dms)
    dms = dms[_idxs]
    total = total[_idxs, :]

    if args.tscrunch is not None:
        # average time samples
        tsamp *= args.tscrunch
        total = total.reshape(total.shape[0], -1, args.tscrunch)
        total = np.mean(total, axis=-1)
        print("Time integrated to: {0:.2f} ms".format(tsamp * 1e3))

    figsize = (10.4, 7.4)

    fig, axs = plt.subplots(
        figsize=figsize,
        nrows=5,
        ncols=1,
        sharex=True,
        sharey=False,
        gridspec_kw={
            "height_ratios": [0.2, 0.2, 0.2, 0.2, 1],
            "hspace": 0,
        },
    )

    plot_range = np.arange(total.shape[1]) * tsamp

    axs[0].plot(plot_range, total[total.shape[0] - 1, :], color="black", zorder=3)
    axs[0].set_ylabel("{0:.1f} DM".format(dms[-1]))

    axs[1].plot(plot_range, total[total.shape[0] // 2, :], color="black", zorder=3)
    axs[1].set_ylabel("{0:.1f} DM".format(dms[len(dms) // 2]))

    axs[2].plot(plot_range, total[0, :], color="black", zorder=3)
    axs[2].set_ylabel("{0:.1f} DM".format(dms[0]))

    ampmax = np.max(total, axis=0)

    axs[3].plot(plot_range, ampmax, color="C0", zorder=3)

    bin_width_x = tsamp
    bin_width_y = abs(dms[1] - dms[0])

    # check that the dm step is always the same
    assert np.all(np.isclose(np.abs(np.diff(dms)), bin_width_y))

    extent = (
        0 - 0.5 * bin_width_x,
        (total.shape[1] - 1) * bin_width_x + 0.5 * bin_width_x,
        dms[0] - 0.5 * bin_width_y,
        dms[-1] + 0.5 * bin_width_y,
    )

    # colour scale selection
    vmin = np.min(total)
    vmax = 0.75 * np.max(total)
    print(vmin, vmax)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax, clip=True)

    axs[4].imshow(
        total, cmap="Greys", aspect="auto", extent=extent, origin="lower", norm=norm
    )

    axs[4].set_xlabel("Time (s)")
    axs[4].set_ylabel("Trial DM")

    fig.align_ylabels()

    fig.tight_layout()

    fig.savefig("candidate_fdmt_timeline.png", bbox_inches="tight", dpi=300)

    if not args.output:
        plt.show()

    print("All done.")


if __name__ == "__main__":
    main()
