#! /usr/bin/env python3
"""
File:           fig7_contact_map.py
Description:    Plot the average contact maps for the 0.1% hardest and median
                difficult instances. Compute the average contact order and show
                as well.
"""
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
import scienceplots  # noqa: F401
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.axes import Axes

from prospr import Protein

os.chdir(Path(__file__).parent)
sys.path.insert(0, "../..")
from common import (  # noqa: E402
    DATASET_NAME_HRATIO,
    DATASET_NAME_RANDOM,
    ROOT_DIR,
    load_results_hardest,
    load_results_median,
)


def _create_heatmap(df: pd.DataFrame, length: int):
    heatmap = np.zeros((length, length))
    for _, row in df.iterrows():
        sequence = row["sequence"]
        protein = Protein(sequence)
        moves = [int(s) for s in row["hash"][1:-1].split(",")]
        for move in moves:
            protein.place_amino(move)
        bond_pairs = protein.get_bonds()
        for i, j in bond_pairs:
            heatmap[i, j] += 1

    # Index pair probability as fraction.
    heatmap /= sum(map(sum, heatmap))
    return heatmap


def _apply_triangle_mask(df_upper: pd.DataFrame, df_lower: pd.DataFrame):
    """Apply triangle mask to two DataFrames and merge."""
    # Mask upper triangle, then overwrite df_upper with df_lower where mask is
    # False.
    upper_mask = np.tri(df_upper.shape[0], k=1)
    upper_triangle = np.where(upper_mask, df_upper, np.nan)
    lower_triangle = np.where(upper_mask, np.nan, df_lower)
    return upper_triangle, lower_triangle


def gradient_image(ax, direction=0.3, cmap_range=(0, 1), **kwargs):
    """
    Draw a gradient image based on a colormap.
    Copied from the official Matplotlib gallary example:
        https://matplotlib.org/stable/gallery/lines_bars_and_markers/gradient_bar.html

    Parameters
    ----------
    ax : Axes
        The Axes to draw on.
    direction : float
        The direction of the gradient. This is a number in
        range 0 (=vertical) to 1 (=horizontal).
    cmap_range : float, float
        The fraction (cmin, cmax) of the colormap that should be
        used for the gradient, where the complete colormap is (0, 1).
    **kwargs
        Other parameters are passed on to `.Axes.imshow()`.
        In particular, *cmap*, *extent*, and *transform* may be useful.
    """
    phi = direction * np.pi / 2
    v = np.array([np.cos(phi), np.sin(phi)])
    X = np.array([[v @ [1, 0], v @ [1, 1]],
                  [v @ [0, 0], v @ [0, 1]]])
    a, b = cmap_range
    X = a + (b - a) / X.max() * X
    im = ax.imshow(X, interpolation='bicubic', clim=(0, 1),
                   aspect='auto', **kwargs)
    return im


def plot_colormap(ax):
    """Plot generic colormap scaling from 0 till 100%. """
    # Add 5 subplots inside this axes, one for each color.
    cmap_axs = []
    x0 = 0.0
    num_lengths = 5
    for _ in range(num_lengths):
        cmap_axs.append(ax.inset_axes((x0, 0.0, 1/num_lengths, 1)))
        x0 += 1/num_lengths

    # Add color gradients to each subplot.
    for i in range(num_lengths):
        cmap = LinearSegmentedColormap.from_list(
            f"white_to_C{i}",
            [(1, 1, 1), plt.get_cmap("tab10")(i)]
        )
        gradient_image(
            cmap_axs[i], direction=0, extent=(0, 1, 0, 1),
            transform=ax.transAxes, cmap=cmap
        )
        # cmap_axs[i].xaxis.set_ticklabels([])
        cmap_axs[i].xaxis.set_visible(False)
        cmap_axs[i].yaxis.set_ticklabels([])

    # Style axes.
    ax.set_title(r"\,", fontsize=18)
    ax.xaxis.set_ticks([])
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0))
    ax.yaxis.tick_right()
    ax.set_ylabel("Fraction of Contacts")
    ax.yaxis.set_label_position("left")


def plot_contact_map(
    df_hardest: pd.DataFrame,
    df_median: pd.DataFrame,
    length: int,
    ax: Axes,
    color_idx: int = 0,
):
    """Plot contact map in grid format."""
    # If length is -1, plot colormap.
    if length == -1:
        plot_colormap(ax)
        return

    # Slice DataFrame on length and report number of instances.
    df_hardest = df_hardest[df_hardest["length"] == length]
    df_median = df_median[df_median["length"] == length]
    print(
        f"Plotting {len(df_hardest)} hard and {len(df_median)} median instances"
        + f" of lenght {length}"
    )

    # Create heatmap for hardest instances.
    heatmap_hardest = _create_heatmap(df_hardest, length)
    heatmap_median = _create_heatmap(df_median, length)
    heatmaps = _apply_triangle_mask(heatmap_hardest, heatmap_median)

    for heatmap in heatmaps:
        # Setup and plot heatmap.
        base_color = plt.get_cmap("tab10")(color_idx)
        cmap = LinearSegmentedColormap.from_list(
            f"white_to_C{color_idx}", [(1, 1, 1), base_color]
        )
        ax.imshow(heatmap, cmap=cmap, origin="lower", aspect="equal")

    # Add diagonal line to separate hardest from median.
    h, w = heatmaps[0].shape
    ax.plot([-0.5, w - 0.5], [-0.5, h - 0.5], color="black", linewidth=1)

    # Format axes.
    ax.grid(False)
    ax.minorticks_off()
    step = 5
    if length == 10:
        step = 2
    elif length == 15:
        step = 3
    ticks = list(range(0, length, step))[1:]
    tick_labels = [str(i) for i in ticks]
    ax.set_xticks(ticks, tick_labels)
    ax.set_yticks(ticks, tick_labels)


def _get_contact_order(row) -> float:
    """Compute the contact order for the given DataFrame row."""
    # Compute and return the number of unqiue bondable amino acids.
    p = Protein(row["sequence"])
    p_hash = [int(x.strip()) for x in row["hash"].strip()[1:-1].split(",")]
    p.set_hash(p_hash)
    return p.get_contact_order()


def compute_contact_order(df: pd.DataFrame, length: int) -> float:
    """Compute the average contact of a specific length from the DataFrame."""
    df_length = df[df["length"] == length].copy()
    df_length["contact_order"] = df_length.apply(_get_contact_order, axis=1)
    return df_length["contact_order"].mean()


def save_fig(fig, figname: str):
    """Generic figure save function."""
    # Save .eps for LaTeX document.
    dpi=1000
    base_output_path = ROOT_DIR / "figures/bigger"
    output_path_eps = base_output_path / f"{figname}.eps"
    output_path_eps.parent.mkdir(exist_ok=True, parents=True)
    print(f"Saving figure to {output_path_eps}")
    fig.savefig(output_path_eps, dpi=dpi)

    # Save .png for local viewing.
    output_path_png = output_path_eps.with_suffix(".png")
    print(f"Saving figure to {output_path_png}")
    output_path_png.parent.mkdir(exist_ok=True, parents=True)
    fig.savefig(output_path_png, dpi=dpi)


def main():
    """Entrypoint for code."""
        # Create figure for both the H-ratio and random results.
    dataset_names = [DATASET_NAME_HRATIO, DATASET_NAME_RANDOM]
    for dataset_name in dataset_names:
        # Get results and only save 0.1% hardest and 0.1% median of each length.
        perc = 0.1
        df_hardest = load_results_hardest(dataset_name, fraction=perc / 100)
        df_median = load_results_median(dataset_name, fraction=perc / 100)

        # Setup plotting style.
        plt.style.use(["science"])
        plt.rcParams.update({"font.size": 7})
        width = 7.48
        lengths = sorted([l for l in df_hardest["length"].unique() if l <= 30])
        ncols = len(lengths) + 1
        nrows = 1

        # Add a bogus length for the colormap.
        ax_width=1/ncols - 0.0075
        width_ratios = [ax_width for _ in range(len(lengths))]
        width_ratios.append(0.0075 * len(lengths))
        lengths.append(-1)

        # Setup figure and plot data.
        fig, axs = plt.subplots(nrows, ncols, width_ratios=width_ratios)
        for i, (ax, length) in enumerate(zip(axs.flatten(), lengths)):
            # Plot contact map with length in title.
            ax.set_title(f"Length {length}", y= 1.08)
            plot_contact_map(df_hardest, df_median, length, ax, i)

            # Compute avg. CO and add annotation.
            if length != -1:
                co_hardest = compute_contact_order(df_hardest, length)
                co_median = compute_contact_order(df_median, length)
                print(f"\tCO_H: {co_hardest}\n\tCO_M: {co_median}")

                # Construct text to plot.
                text = r"$\text{CO}_{\text{H}}$="
                text += f"{round(co_hardest, 2):.2f}"
                text += r"\,\textbar\,$\text{CO}_{\text{M}}$="
                text += f"{round(co_median, 2):.2f}"

                ax.text(
                    0.5,
                    1.025,
                    text,
                    transform=ax.transAxes,
                    verticalalignment="bottom",
                    horizontalalignment="center",
                    fontsize=6,
                    color="k",
                )

        # Style figure.
        fig.supxlabel("Amino Acid Index", y=0.05)
        fig.supylabel("Amino Acid Index")
        fig.set_size_inches(width, nrows * 1.25 + 0.5, forward=True)
        # fig.subplots_adjust(wspace=0.35, hspace=0.25)
        fig.tight_layout()

        # Save figures.
        save_fig(fig, f"fig7_contact_map_{dataset_name}")


if __name__ == "__main__":
    main()
