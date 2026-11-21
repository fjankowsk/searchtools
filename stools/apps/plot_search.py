#
#   2023 - 2026 Fabian Jankowski
#   Plot single pulse search output.
#

import argparse
import glob

import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    """
    Parse the commandline arguments.

    Returns
    -------
    args: populated namespace
        The commandline arguments.
    """

    parser = argparse.ArgumentParser(
        description="Plot single pulse search output.",
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

    args = parser.parse_args()

    return args


def main():
    # handle command line arguments
    args = parse_args()

    files = glob.glob("*.singlepulse")

    cands = []

    for fname in files:
        with open(fname, "r") as fd:
            lines = fd.readlines()
            for line in lines:
                if line.startswith("#"):
                    pass
                else:
                    cands.append(line)

    with open("sp_cands.csv", "w") as fd:
        fd.write("".join(cands))

    names = ["dm", "snr", "time", "sample", "downsamp"]
    df = pd.read_csv("sp_cands.csv", names=names, delim_whitespace=True)
    df.info()

    fig = plt.figure()
    ax = fig.add_subplot()

    size = 100.0 * df["snr"] / df["snr"].max()
    size[size < 16.0] = 16.0

    sc = ax.scatter(
        df["time"],
        df["dm"],
        c=df["snr"],
        s=size,
        edgecolor="black",
        lw=0.1,
        marker=".",
        zorder=3,
    )

    plt.colorbar(sc, label="S/N")

    ax.grid()
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("DM")

    fig.tight_layout()

    fig.savefig("candidate_dm_timeline.png", bbox_inches="tight", dpi=300)

    fig = plt.figure()
    ax = fig.add_subplot()

    mask = df["dm"] > 0.5
    ax.hist(
        df.loc[mask, "dm"], bins="auto", color="black", lw=2, histtype="step", zorder=3
    )

    ax.grid()
    ax.set_xlabel("DM")
    ax.set_ylabel("Number")

    fig.tight_layout()

    if not args.output:
        plt.show()

    print("All done.")


if __name__ == "__main__":
    main()
