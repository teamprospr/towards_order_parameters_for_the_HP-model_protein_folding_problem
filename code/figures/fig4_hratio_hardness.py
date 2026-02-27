#! /usr/bin/env python3
"""
File:           fig4_hratio_hardness.py
Description:    Create a heatmap of the H-ratio against hardness percentile with
                the number of recursions as a black line on a 2nd Y-axis.
"""
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scienceplots  # noqa: F401
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import PercentFormatter

os.chdir(Path(__file__).parent)
sys.path.insert(0, "..")
from common import (  # noqa: E402
    DATASET_NAME_HRATIO,
    DATASET_NAME_RANDOM,
    ROOT_DIR,
    load_results,
)


def plot_results_length(
    df: pd.DataFrame,
    length: int,
    ax: plt.Axes,
    color_idx: int = 0,
    annotate_y_axis: str = "both",
) -> None:
    ax_h_ratio = ax
    ax_hardness = ax_h_ratio.twinx()
    max_recusions = df["placed"].max()
    df = df[df["length"] == length]
    print(f"Plotting {len(df)} instances of lenght {length}")
    df = df.sort_values(by=["placed"])
    x = [x / (len(df) - 1) for x in range(len(df))]
    y = df["h_ratio"]
    c = plt.rcParams["axes.prop_cycle"].by_key()["color"][color_idx]
    r = np.corrcoef(x, y)[0][1]
    plt.text(
        0.07,
        0.925,
        f"$r = {r:.2f}$",
        transform=plt.gca().transAxes,
        verticalalignment="top",
        horizontalalignment="left",
        color=c,
        bbox=dict(
            boxstyle="round,pad=0.3", edgecolor="black", facecolor="white",
            alpha=0.75
        ),
    )
    cmap = LinearSegmentedColormap.from_list(f"h-ratio-l{length}", ["white", c])
    cmap = LinearSegmentedColormap.from_list(f"h-ratio-l{length}", ["white", c])
    ax_h_ratio.hist2d(x, y, bins=25, cmap=cmap)
    ax_h_ratio.set_ylabel(r"\texttt{H}-ratio")
    #ax_h_ratio.yaxis.set_major_formatter(PercentFormatter(1))
    ax_h_ratio.set_ylim(0, 1)
    ax_h_ratio.xaxis.set_major_formatter(PercentFormatter(1))
    ax_h_ratio.set_xlim(-0.05, 1.05)
    # ax_h_ratio.set_xlabel("Hardness Percentile")
    if annotate_y_axis not in ["left", "both"]:
        ax_h_ratio.set_ylabel("")
        ax_h_ratio.set_yticks([])

    y = df["placed"]
    ax_hardness.plot(x, y, label=f"Length {length}", color="k")
    ax_hardness.set_ylabel("Recursions")
    ax_hardness.set_yscale("log")
    ax_hardness.set_ylim(1, max_recusions)
    if annotate_y_axis not in ["right", "both"]:
        ax_hardness.set_ylabel("")
        ax_hardness.set_yticks([])


def main():
    dataset_names = [DATASET_NAME_HRATIO, DATASET_NAME_RANDOM]
    for dataset_name in dataset_names:
        # Load results and remove any points for lengths larger than 30.
        df = load_results(dataset_name)
        df = df[df["length"] <= 30]
        print(f"Plotting {len(df)} instances")

        plt.style.use(["science"])
        # Style following:
        #   https://www.elsevier.com/about/policies-and-standards/author/artwork-and-media-instructions/artwork-sizing
        plt.rcParams.update({"font.size": 7})
        width = 7.48
        dpi = 1000
        lengths = sorted(df["length"].unique())
        ncols = len(lengths)  # 6
        nrows = len(lengths) // ncols
        if nrows * ncols < len(lengths):
            nrows += 1
        fig, axs = plt.subplots(nrows, ncols)
        for i, (ax, length) in enumerate(zip(axs.flatten(), lengths)):
            ax.set_title(f"Length {length}")
            annotate_y_axis = "none"
            if i == 0:
                annotate_y_axis = "left"
            if i == len(lengths) - 1:
                annotate_y_axis = "right"
            plot_results_length(df, length, ax, i, annotate_y_axis)
        fig.supxlabel("Hardness Percentile")
        fig.set_size_inches(width, nrows * 1.25 + 0.5, forward=True)
        fig.subplots_adjust(wspace=0.35, hspace=0.25)
        fig.tight_layout()

        # Save .eps for LaTeX document.
        base_output_path = ROOT_DIR / "figures"
        output_path_eps = base_output_path / \
            f"fig4_hratio_hardness_{dataset_name}.eps"
        print(f"Saving figure to {output_path_eps}")
        output_path_eps.parent.mkdir(exist_ok=True, parents=True)
        fig.savefig(output_path_eps, dpi=dpi)

        # Save .png for local viewing.
        output_path_png = output_path_eps.with_suffix(".png")
        print(f"Saving figure to {output_path_png}")
        output_path_png.parent.mkdir(exist_ok=True, parents=True)
        fig.savefig(output_path_png, dpi=dpi)


if __name__ == "__main__":
    main()
