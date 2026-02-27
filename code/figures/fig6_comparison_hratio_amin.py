#! /usr/bin/env python3
"""
File:           fig6_comparison_hratio_amin.py
Description:    Plot a side-by-side comparison of using H-ratio and Amin as
                primary order parameters.
"""
# System packages.
import os
import sys
from pathlib import Path

# Scientific packages.
import pandas as pd
import numpy as np

# Plotting packages.
import scienceplots  # noqa: F401
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.axes import Axes
from matplotlib.legend_handler import HandlerTuple

# Project packages.
os.chdir(Path(__file__).parent)
sys.path.insert(0, "..")
from common import (  # noqa: E402
    DATASET_NAME_HRATIO,
    ROOT_DIR,
    load_results,
)


def save_fig(fig, figname: str):
    """Generic figure saving function."""
    # Save .eps for LaTeX document.
    dpi=1000
    base_output_path = ROOT_DIR / "figures"
    output_path_eps = base_output_path / f"{figname}.eps"
    output_path_eps.parent.mkdir(exist_ok=True, parents=True)
    print(f"Saving figure to {output_path_eps}")
    fig.savefig(output_path_eps, dpi=dpi)

    # Save .png for local viewing.
    output_path_png = output_path_eps.with_suffix(".png")
    print(f"Saving figure to {output_path_png}")
    output_path_png.parent.mkdir(exist_ok=True, parents=True)
    fig.savefig(output_path_png, dpi=dpi)


def lighten_color(color, amount=0.5):
    """
    Lightens the given color by mixing it with white.

    Parameters:
        color: Matplotlib color string, hex, or RGB tuple
        amount: 0 = original color, 1 = white

    Returns:
        Lightened RGB tuple
    """
    try:
        c = mcolors.cnames[color]
    except KeyError:
        c = color
    rgb = mcolors.to_rgb(c)
    lightened = [1 - (1 - x) * (1 - amount) for x in rgb]
    return lightened


def plot_max_trace(df: pd.DataFrame, order_param: str, xlabel: str, ax: Axes, plot_legend: str):
    """Plot a lineplot of the maximum values."""
    # Remove any points for lengths larger than 30.
    lengths = sorted(
        [length for length in df["length"].unique() if length <= 30],
    )
    handles = []
    labels = []
    previous_max = None

    for i, length in enumerate(lengths):
        df_local = df[df["length"] == length]
        color = f"C{i}"

        # Line plot trace of maximum values.
        df_max = df_local.groupby(order_param, as_index=False)["placed"].max()
        df_max.sort_values(order_param, inplace=True)
        (lines,) = ax.plot(
            df_max[order_param],
            df_max["placed"],
            c=color,
            label=f"{length}",
            zorder=2
        )

        # Scatter plot the data in a lighter color.
        y_scatter = df_local["placed"]
        if previous_max is not None:
            # Align previous max by order_param
            prev_max_interp = np.interp(df_local[order_param], previous_max[order_param], previous_max['placed'])
            y_scatter = np.where(y_scatter > prev_max_interp, y_scatter, np.nan)

        points = ax.scatter(
            df_local[order_param],
            y_scatter,
            color=lighten_color(color),
            s=0.5,
            alpha=0.3,
        )
        # Set alpha to 1.0 for better visibility in legend.
        points.set_alpha(1.0)

        # Keep track of handles and labels.
        labels.append(str(int(length)))
        handles += [(lines, points)]

        # Update previous max for clipping scatter.
        previous_max = df_max

    # Plot legend if specified.
    if plot_legend == "left":
        box = ax.get_position()
        ax.set_position((box.x0, box.y0, box.width * 0.8, box.height))

        # Put a legend to the right of axis.
        ax.legend(
            handles,
            labels,
            handler_map={tuple: HandlerTuple(ndivide=None)},
            ncols=1,
            loc="center right",
            frameon=True,
            bbox_to_anchor=(0.0, 0.5),
        )
    elif plot_legend == "right":
        box = ax.get_position()
        ax.set_position((box.x0, box.y0, box.width * 0.8, box.height))

        # Reverse legend order to be in numerical order.
        ax.legend(
            handles[::-1],
            labels[::-1],
            handler_map={tuple: HandlerTuple(ndivide=None)},
            ncols=1,
            loc="center left",
            frameon=True,
            bbox_to_anchor=(1.0, 0.5),
        )
    else:
        # Legend is plotted with the right axis, so the left plots Ylabel.
        ax.set_ylabel("Recursions")

    # Style axis.
    ax.set_yscale("log")
    ax.set_ylim(bottom=1, top=1e12)
    ax.set_xlabel(xlabel)


def plot_comparison(df: pd.DataFrame):
    """Plot a comparison of the Hratio with Amin order parameters."""
    # Create subplots for both order parameters.
    # Setup figure details, following:
    #   https://www.elsevier.com/about/policies-and-standards/author/artwork-and-media-instructions/artwork-sizing
    plt.style.use(["science"])
    plt.rcParams.update({"font.size": 8})
    nrows, ncols = 1, 2
    sns.set_style("whitegrid", {"font.family": "serif"})
    fig, axs = plt.subplots(nrows, ncols, sharey="row")

    # Plot H-ratio and random figures.
    legend_pos = "right"
    plot_max_trace(df, "h_ratio", r"\texttt{H}-ratio", axs[0], "")
    plot_max_trace(df, "A_min_norm", r"Normalized $A_{\text{min}}$", axs[1], legend_pos)

    # Style overall figure and save.
    width = 7.48
    fig.set_size_inches(width, width / 2.75, forward=True)
    fig.subplots_adjust(wspace=0.25)
    fig.tight_layout()
    save_fig(fig, f"fig6_comparison_hratio_amin_legend_{legend_pos}")


def main():
    # Set styling of figures.
    plt.style.use(["science"])
    plt.rcParams.update({
        "legend.frameon": True,
        "legend.framealpha": 1,
        "legend.facecolor": "white",
        "legend.edgecolor": "black",
        "font.size": 7,
    })

    # Load results data and add Amin order parameter values.
    df = load_results(DATASET_NAME_HRATIO)
    df["h_count"] = df["sequence"].str.count("H")
    df["num_p_singlets"] = df["sequence"]\
        .map(lambda s: sum([1 if s[i] == "H" and s[i+1] == "P" and s[i+2] == "H"
                            else 0 for i in range(len(s)-2) ]))
    df["A_min"] = df["h_count"] + df["num_p_singlets"]
    df["A_min_norm"] = df["A_min"] / df["sequence"].str.len()

    # Create and save the figure.
    plot_comparison(df)


if __name__ == "__main__":
    main()
