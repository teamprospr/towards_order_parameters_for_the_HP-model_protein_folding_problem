#! /usr/bin/env python3
"""
File:           fig8_hardest_each_length.py
Description:    Visualize the hardest instances for each length from the H-ratio
                dataset.
"""
# System packages.
import os
import sys
from pathlib import Path

# Scientific packages.
import pandas as pd
from prospr import Protein
from prospr.visualize import plot_protein

# Plotting packages.
import scienceplots  # noqa: F401
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

# Project packages.
os.chdir(Path(__file__).parent)
sys.path.insert(0, "..")
from common import (  # noqa: E402
    DATASET_NAME_HRATIO,
    ROOT_DIR,
    load_results,
)


def save_fig(fig, figname: str):
    """Generic function for saving figures."""
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


def plot_instance(hardest_protein: pd.Series, ax: Axes, plot_legend: bool=False):
    """Plot the given proteins in the given order."""
    # Plot protein.
    protein = Protein(hardest_protein["sequence"])
    protein.set_hash(list(map(int, hardest_protein["hash"][1:-1].split(","))))

    # Add legend only if specified.
    if not plot_legend:
        plot_protein(
            protein,
            style="paper",
            markersize=35,
            linewidth=1.25,
            ax=ax,
            show=False,
            legend=False,
            annotate_first=True,
        )
    else:
        plot_protein(
            protein,
            style="paper",
            markersize=35,
            linewidth=1.25,
            ax=ax,
            show=False,
            annotate_first=True,
            legend=True,
            legend_style="outer",
            fontsize=6,
        )

    # Style figure.
    ax.set_title(
        f"Length {hardest_protein["length"]}", y=1.075
    )
    recursions = f"{hardest_protein["placed"]:,}".replace(",", r"\,")
    ax.text(
        0.5,
        1.075,
        f"{recursions} recursions",
        horizontalalignment="center",
        verticalalignment="center",
        transform = ax.transAxes,
        fontsize="small"
    )

    # Morph figure to make all segments of equal length.
    ax.axis('equal')

def main():
    """Entrypoint."""
    # Load the data.
    df = load_results(DATASET_NAME_HRATIO)
    df.sort_values("placed", inplace=True)
    lengths = df[df["length"] <= 30]["length"].unique()

    # Setup figure with a window per length.
    # Style following:
    #   https://www.elsevier.com/about/policies-and-standards/author/artwork-and-media-instructions/artwork-sizing
    plt.style.use(["science"])
    plt.rcParams.update({
        "legend.frameon": True,
        "legend.framealpha": 1,
        "legend.facecolor": "white",
        "legend.edgecolor": "black",
        "font.size": 7,
    })
    ncols = len(lengths)
    nrows = 1

    # Create subplots and fill in with hardest instances.
    plot_legend = False
    fig, axs = plt.subplots(nrows, ncols)
    for i, (ax, length) in enumerate(zip(axs.flatten(), lengths)):
        df_hardest = df[df["length"] == length].tail(1).iloc[0]
        if i == len(lengths) - 1:
            plot_legend = True
        plot_instance(df_hardest, ax, plot_legend=plot_legend)

    width, height = 7.48, 2.5
    fig.set_size_inches(width, height, forward=True)  # Height was 1.75
    fig.subplots_adjust(wspace=0.0)
    fig.tight_layout()
    save_fig(fig, f"fig8_hardest_each_length")


if __name__ == "__main__":
    main()
