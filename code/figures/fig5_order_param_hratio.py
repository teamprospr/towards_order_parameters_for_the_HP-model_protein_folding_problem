#! /usr/bin/env python3
"""
File:           fig5_order_param_hratio.py
Description:    Apply the H-ratio as primary order parameter and fit the maxima
                per H-ratio bin with the GEV distribution. Add a Cheeseman inset
                for comparison.
"""
# System packages.
import os
import sys
from pathlib import Path

# Scientific packages.
import numpy as np
import pandas as pd
from scipy.stats import genextreme

# Plotting packages.
import scienceplots  # noqa: F401
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors

os.chdir(Path(__file__).parent)
sys.path.insert(0, "..")
from common import (  # noqa: E402
    DATASET_NAME_HRATIO,
    DATASET_NAME_RANDOM,
    ROOT_DIR,
    load_results,
)


def save_fig(fig, figname: str):
    """Generic save function for figures."""
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


def r2_score(y, y_pred):
    """
    Compute the R^2 metric.
    :param  [float]    y:  Observations
    :param  [float]    y_pred:  Predictions.
    """
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - (ss_res / ss_tot)
    return r2


def fit_hardest(df: pd.DataFrame):
    """Fit the upper envelope with the given scipy.stats function."""
    # Get upper envelope of max values per H-ratio.
    df_max = df.groupby("h_ratio", as_index=False)["placed"].max()
    df_max.sort_values("h_ratio", inplace=True)

    # Transform maximums into histogram data.
    sample_size = 1000
    placed_sum = df_max["placed"].sum()
    hist_data = []
    for row in df_max.iterrows():
        hist_data.extend([
            row[1]["h_ratio"]
            for _ in range(round(sample_size * (row[1]["placed"] / placed_sum)))
        ])

    # Fit histogram data to the GEV.
    popt = genextreme.fit(hist_data)
    print(f"popt:\t{popt}")

    # Sample PDF and scale to original data.
    xs = np.linspace(0.0, 1.0, 50)
    ys = genextreme.pdf(xs, *popt)
    ys *= df_max["placed"].max() / max(ys)

    # Compute the R^2 goodness-of-fit score.
    ys_pred = genextreme.pdf(df_max["h_ratio"], *popt)
    ys_pred *= df_max["placed"].max() / max(ys_pred)
    r2 = r2_score(df_max["placed"], ys_pred)
    print(f"R^2: {r2:.2f}")

    return xs, ys, r2


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


def plot_cheeseman(ax: Axes):
    """Plot Cheeseman results on provided axis."""
    # Setup subplot axis to be empty with Title above and load image.
    img = mpimg.imread(f"./images/cheeseman_4-col_no_title.png")
    ax.set_title(f"Cheeseman's 4-COL results")
    ax.axis("off")

    # Add new axis which is slightly bigger to show Cheeseman's results.
    left, bottom, width, height = ax.get_position().bounds
    ax2 = plt.axes((left+0.0485, bottom+0.005, width-0.005, height-0.005))
    ax2.axis("off")
    ax2.imshow(img)


def plot_results_length(
    df: pd.DataFrame,
    length: int,
    ax: Axes,
    color_idx: int = 0,
    annotate_y_axis: str = "both",
) -> None:
    # If length is -1, plot Cheeseman figure.
    if length == -1:
        plot_cheeseman(ax)
        return

    # Plot only the given length.
    df = df[df["length"] == length]
    print(f"Plotting {len(df)} instances of lenght {length}")
    color = f"C{color_idx}"

    # Scatter plot H-ratio and number of recursions.
    x = df["h_ratio"]
    y = df["placed"]
    ax.scatter(x, y, c=color, s=0.5, rasterized=True)

    # Overlap points with hardest ones.
    df_max = df.groupby("h_ratio", as_index=False)["placed"].max()
    df_max.sort_values("h_ratio", inplace=True)
    x = df_max["h_ratio"]
    y = df_max["placed"]
    ax.scatter(
        x,
        y,
        s=7,
        marker="*",
        edgecolors="k",
        facecolors=color,
        linewidths=0.25,
        rasterized=True
    )

    # Fit the envelope of max values and plot the line.
    xs_fit, ys_fit, r2 = fit_hardest(df)
    ax.plot(xs_fit, ys_fit, "--", lw=0.75, c=lighten_color(color), zorder=-1)

    # Add R^2 score as text.
    x_coord = 0.93 if color_idx > 1 else 0.07
    ha = "right" if color_idx > 1 else "left"
    ax.text(
        x_coord,
        0.925,
        f"$R^2 = {r2:.2f}$",
        transform=ax.transAxes,
        verticalalignment="top",
        horizontalalignment=ha,
        fontsize=7,
        color=color,
        bbox=dict(
            boxstyle="round,pad=0.3",
            edgecolor="black",
            facecolor="white",
            alpha=0.75,
        ),
    )

    # Style figures.
    ax.set_ylabel("Normalized Recursions")
    ax.set_yticks([0, y.max() // 2, y.max()])
    ax.set_yticklabels(["0", "0.5", "1"])

    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%g'))
    ax.set_xlim(-0.05, 1.05)
    if annotate_y_axis not in ["left"]:
        ax.set_ylabel("")
        ax.set_yticklabels([])


def main():
    dataset_names = [DATASET_NAME_HRATIO, DATASET_NAME_RANDOM]
    for dataset_name in dataset_names:
        # Load results and remove any points for lengths larger than 30.
        df = load_results(dataset_name)
        df = df[df["length"] <= 30]
        print(f"Plotting {len(df)} instances")

        # Setup plotting style.
        plt.style.use(["science"])
        plt.rcParams.update({
            "legend.frameon": True,
            "legend.framealpha": 1,
            "legend.facecolor": "white",
            "legend.edgecolor": "black",
        })
        plt.rcParams.update({"font.size": 7})
        # width = 7.48
        width = 6.48
        lengths = sorted(df["length"].unique())
        ncols = 3
        nrows = 2

        # Append -1 to signify plotting Cheeseman.
        lengths.append(-1)

        # Setup figure, plot data, and style axes.
        fig, axs = plt.subplots(nrows, ncols)
        for i, (ax, length) in enumerate(zip(axs.flatten(), lengths)):
            ax.set_title(f"Length {length}")
            annotate_y_axis = "none"
            if i % ncols == 0:
                annotate_y_axis = "left"
            if i == len(lengths) - 1:
                annotate_y_axis = "right"
            plot_results_length(df, length, ax, i, annotate_y_axis)
        fig.supxlabel(r"\texttt{H}-ratio", y=0.05)
        fig.set_size_inches(width, nrows * 1.25 + 0.5, forward=True)
        fig.subplots_adjust(wspace=0.35, hspace=0.25)
        fig.tight_layout()

        # Save figure.
        save_fig(fig, f"fig5_order_param_hratio_{dataset_name}")


if __name__ == "__main__":
    main()
