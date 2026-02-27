#! /usr/bin/env python3
"""
File:           fig2_dataset_hratios.py
Description:    Compute the H-ratio distributions for each dataset and create a
                comparison figure.
"""
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scienceplots  # noqa: F401
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.axes import Axes

os.chdir(Path(__file__).parent)
sys.path.insert(0, "..")
from common import (  # noqa: E402
    DATASET_NAME_HRATIO,
    DATASET_NAME_RANDOM,
    ROOT_DIR,
    load_results
)

# Define text sizes.
FONTSIZE_LABELS = 10
FONTSIZE_LEGEND = 9
FONTSIZE_TICKS = 8
FONTSIZE_ANNOTATION = 6


def apply_hratios_bins(df: pd.DataFrame) -> pd.DataFrame:
    """Add a column 'h_bin' to the dataframe with the H-ratio bin."""
    hratio_bins = np.arange(0, 1.01, 0.1).tolist()
    df["h_bin"] = pd.cut(
        df["h_ratio"],
        bins=hratio_bins,
        include_lowest=True,
    )
    result = (
        df.groupby(["length", "h_bin"], observed=False)
            .size()
            .reset_index(name="count")
    )
    return result


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


def _add_text_to_plot(ax: Axes, xs, width, sep):
    """Add text annotations to the plot."""
    total_proteins = [1024, 9152, 15422, 20652, 25932]

    # Loop over bin locations.
    for i, [x, num_proteins] in enumerate(zip(xs, total_proteins)):
        # Add R and H to annotate the bars with the datasets.
        ax.text(
            x - width / 2 - 2.5 * sep,
            -4,
            "R",
            fontsize=FONTSIZE_ANNOTATION,
            fontweight="black",
        )
        ax.text(
            x + width / 2 - 0.5 * sep,
            -4,
            "H",
            fontsize=FONTSIZE_ANNOTATION,
            fontweight="black",
        )

        # Add number of proteins above the bars.
        num_prot_text = f"{num_proteins:,}".replace(",", r"\,")
        ax.text(
            x - width / 2 - 0.019 * len(num_prot_text),
            101,
            num_prot_text,
            # fontsize="x-small",
            fontsize=FONTSIZE_TICKS,
            fontweight="black",
        )


def plot_barchart(df_hratio, df_random):
    """Plot a double bar-chart per length."""
    # Setup figure.
    plt.style.use(["science"])
    plt.rcParams.update({
        "legend.frameon": True,
        "legend.framealpha": 1,
        "legend.facecolor": "white",
        "legend.edgecolor": "black",
        "font.size": 7,
    })
    width = 3.54330709  # 90mm
    aspect_ratio = 3 / 4
    fig = plt.figure(figsize=(width, aspect_ratio * width), dpi=1000)
    sns.set_style("whitegrid", {"font.family": "serif"})
    ax = fig.gca()

    # Scalar Mappable for the H-ratio colormaps.
    norm = Normalize(vmin=0.1, vmax=1.0)
    hratio_sm = plt.cm.ScalarMappable(cmap="rocket_r", norm=norm)

    # Setup barchart variables.
    lengths = df_hratio.index
    x = np.arange(len(lengths))
    width = 0.275
    sep = 0.03

    # Plot random dataset.
    bottom1 = np.zeros(len(lengths))
    for col in df_random.columns:
        ax.bar(
            x - sep - width/2,
            df_random[col],
            width,
            bottom=bottom1,
            label=col.right,
            color=hratio_sm.to_rgba(col.right)
        )
        bottom1 += df_random[col].values

    # Plot hratio dataset.
    bottom2 = np.zeros(len(lengths))
    for col in df_hratio.columns:
        ax.bar(
            x + sep + width/2,
            df_hratio[col],
            width,
            bottom=bottom2,
            color=hratio_sm.to_rgba(col.right)
        )
        bottom2 += df_hratio[col].values

    # Annotate plot with text.
    _add_text_to_plot(ax, x, width, sep)

    # Add legend to the right side of the plot.
    box = ax.get_position()
    ax.set_position((box.x0, box.y0, box.width * 0.8, box.height))
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles[::-1],
        labels[::-1],
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        title="H-ratio Bin",
        fontsize=FONTSIZE_LEGEND,
    )

    # Create second x-axis label.
    ax2 = ax.twiny()
    ax2.set_xticks([])
    ax2.set_xticklabels([])
    ax2.set_xlabel("Number of Proteins", labelpad=10, fontsize=FONTSIZE_LABELS)

    # Setup primary axis.
    ax.set_xticks(x)
    ax.xaxis.set_tick_params(labelsize=FONTSIZE_TICKS)
    ax.yaxis.set_tick_params(labelsize=FONTSIZE_TICKS)
    ax.set_xticklabels(lengths, fontsize=FONTSIZE_TICKS)
    ax.set_xlabel("Protein Length", labelpad=3, fontsize=FONTSIZE_LABELS)
    ax.set_ylabel(r"Percentage of Dataset", fontsize=FONTSIZE_LABELS)
    ax.set_ylim(0, 100)

    # Remove tick params.
    ax.tick_params(
        axis="x",
        which="both",
        bottom=False,
        top=False,
    )
    ax2.tick_params(
        axis="x",
        which="both",
        bottom=False,
        top=False,
    )

    plt.tight_layout()
    save_fig(fig, "fig2_dataset_hratios")


def main():
    """Entrypoint for code."""
    # Load the proteins from both datasets limit till length 30.
    df_hratio = load_results(DATASET_NAME_HRATIO)
    df_hratio = df_hratio[df_hratio.length <= 30]
    df_random = load_results(DATASET_NAME_RANDOM)
    df_random = df_random[df_random.length <= 30]

    # Bin hratio values per 0.1 window.
    df_hratio_binned = apply_hratios_bins(df_hratio)
    df_random_binned = apply_hratios_bins(df_random)

    # Pivot dataframes to be plotted by iterating over columns.
    df_hratio_pivot = df_hratio_binned.pivot(
        index="length", columns="h_bin", values="count"
    ).fillna(0)
    df_random_pivot = df_random_binned.pivot(
        index="length", columns="h_bin", values="count"
    ).fillna(0)

    # Ensure identical bin order.
    df_random_pivot = df_random_pivot.reindex(
        columns=df_hratio_pivot.columns, fill_value=0
    )

    # Convert to percentages per length.
    df_hratio_pct = df_hratio_pivot.div(df_hratio_pivot.sum(axis=1), axis=0) \
        * 100
    df_random_pct = df_random_pivot.div(df_random_pivot.sum(axis=1), axis=0) \
        * 100

    # Plot the figure.
    plot_barchart(df_hratio_pct, df_random_pct)


if __name__ == "__main__":
    main()