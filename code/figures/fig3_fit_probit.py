#! /usr/bin/env python3
# System packages.
import math
import os
import sys
from pathlib import Path
from typing import List, Tuple

# Scientific packages.
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.special import erfinv
from scipy.stats import normaltest

# Plotting packages.
import scienceplots  # noqa: F401
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import matplotlib.ticker as mtick
import matplotlib.colors as mcolors
from matplotlib.legend_handler import HandlerTuple
import seaborn as sns

# Project packages.
os.chdir(Path(__file__).parent)
sys.path.insert(0, "..")
from common import (  # noqa: E402
    ROOT_DIR,
)


def get_data(dataset: str):
    """Fetch data."""
    def get_csvs(path: Path) -> List[Path]:
        assert path.is_dir()
        csv_paths = []
        for child_path in path.iterdir():
            if child_path.is_file() and child_path.name.lower().endswith(".csv"):
                csv_paths.append(child_path)
            if child_path.is_dir():
                csv_paths += get_csvs(child_path)
        return csv_paths

    data = {}

    for p in Path("jobs").iterdir():
        if not p.is_dir() or not p.name.startswith(dataset + "_"):
            continue
        try:
            # Folder name must follow "<prefix>_l<lenght>_..."
            length = int(p.name.split("_")[1][1:])

            # Don't parse lengths larger than 30.
            if length > 30:
                continue
        except Exception:
            continue
        filenames = [
            p2
            for p2 in get_csvs(p)
            if (
                p2.name.startswith("HP_")
                and
                # Filter out intermediate results with run number suffix (e.g. ..._r0.csv)
                not p2.stem.split("_")[-1].startswith("r")
            )
        ]

        source_dfs = [pd.read_csv(f).dropna() for f in filenames]
        source_df = pd.concat(
            [df for df in source_dfs if len(df) > 0], ignore_index=True
        )
        data[length] = source_df
    return data


def probit_func(x, alpha, beta):
    """
    Probit function with scaling parameters to fit data to.
    Note:   Log-scale y-values when fitting this function.
    :param  int      x:      Rank of the point in the figure (x-axis).
    :param  float    alpha:  Scaling factor (sqrt(2) in original probit)
    :param  float    beta:   Y-axis shift.
    """
    return alpha * erfinv(2 * x - 1) + beta


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


def fit_probit_to_data(data: dict) -> Tuple[List, pd.DataFrame]:
    """Fit a probit function to all lengths individually."""
    # Setup variables for plotting and fitting.
    lengths = sorted(data.keys())
    lengths_ys = []
    popts = pd.DataFrame(columns=[
        "length", "n_samples", "alpha", "beta", "r2", "pval_normaltest"
    ])

    # Fit the data for each length.
    for length in lengths:
        source_df = data[length]
        ys = sorted(source_df["placed"].values, reverse=False)
        xs = np.array(range(1, len(ys) + 1))

        # Fit data to the scaled probit function.
        # Scale x-data to [0,1] and log-scale y-data for fit.
        machine_precision = np.finfo(float).eps
        x = np.linspace(machine_precision, 1 - machine_precision, len(ys))
        lengths_ys.append(ys.copy())
        log_ys = np.log10(ys)
        p0 = [math.sqrt(2), log_ys[xs[len(xs) // 2]]]
        popt = curve_fit(probit_func, x, log_ys, p0=p0)[0]
        r2 = r2_score(log_ys, probit_func(x, *popt))
        pval_normaltest = normaltest(log_ys)

        # Print statistics on the loaded data.
        print(f"\tOptimal params length {length}:  {popt}")
        popts.loc[length] = {
            "length": length,
            "n_samples": len(ys),
            "alpha": popt[0],
            "beta": popt[1],
            "r2": r2,
            "pval_normaltest": pval_normaltest.pvalue,
        }
    print(f"R^2 scores:  {popts['r2'].values}")
    print(f"Normaltest p-values:  {popts['pval_normaltest'].values}\n")

    # Return y-values and optimal fit parameters per length.
    return lengths_ys, popts


def plot_data(data: dict, title: str, ax_line: Axes, ax_r2: Axes, annotate_y: bool, plot_legend: str):
    """Plot the data on the given two axes. Only plot legend if specified."""
    # Fit probits to each length.
    lengths_ys, popts = fit_probit_to_data(data)

    # Plot the data and fitted probit functions.
    handles = []
    labels = []
    machine_precision = np.finfo(float).eps
    for i, (_, row) in list(enumerate(popts.iterrows())):
        color = f"C{i}"
        light_color = lighten_color(color, amount=0.5)
        ys = lengths_ys[i]
        xs = [i / len(ys) for i in range(len(ys))]
        (length, n_samples, alpha, beta) = row.values[:4]
        labels.append(str(int(length)))
        (instances,) = ax_line.plot(xs, ys, color=light_color)
        x = np.linspace(machine_precision, 1 - machine_precision, int(n_samples))
        # Un-transform and clamp (probit was fitted on log transformed data)
        y = [max(1, v) for v in 10 ** probit_func(x, alpha, beta)]
        (fit,) = ax_line.plot(x, y, ":", color=color)
        handles += [(instances, fit)]

    # Plot legend if specified.
    if plot_legend == "left":
        box = ax_line.get_position()
        ax_line.set_position((box.x0, box.y0, box.width * 0.8, box.height))

        # Put a legend to the right of axis.
        ax_line.legend(
            handles,
            labels,
            handler_map={tuple: HandlerTuple(ndivide=None)},
            ncols=1,
            loc="center right",
            frameon=True,
            bbox_to_anchor=(0.0, 0.5),
        )
    elif plot_legend == "right":
        box = ax_line.get_position()
        ax_line.set_position((box.x0, box.y0, box.width * 0.8, box.height))

        # Put a legend to the right of axis.
        ax_line.legend(
            handles,
            labels,
            handler_map={tuple: HandlerTuple(ndivide=None)},
            ncols=1,
            loc="center left",
            frameon=True,
            bbox_to_anchor=(1.0, 0.5),
        )

    # Plot R^2 values
    lengths = sorted(data.keys())
    r2_values = popts['r2'].values
    colors = [f"C{i}" for i in range(len(lengths))]
    ax_r2.scatter(lengths, r2_values, fc=colors, ec="white", s=15)

    # Style axes.
    ax_line.set_yscale("log")
    ax_line.set_yscale("log")
    ax_line.set_ylim(bottom=1, top=1e12 / 1.5)
    ax_line.set_xlabel("Hardness Percentile")
    ax_r2.set_xticklabels([])
    ax_r2.set_ylim(0.75, 1.05)
    ax_r2.set_title(title)
    if annotate_y:
        ax_line.set_ylabel("Recursions")
        ax_r2.set_ylabel("$R^2$")


def save_fig(fig, figname: str):
    """Plot given figure with name."""
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


def plot_probits():
    """
    Fit a Probit to each length, store the best parameters and create a plot
    showing the development of the params across lengths.
    """
    # Load data for both figures.
    random_data = get_data("random")
    hratio_data = get_data("hratio")

    # Create subplots for lines and R^2 values.
    # Setup figure details, following:
    #   https://www.elsevier.com/about/policies-and-standards/author/artwork-and-media-instructions/artwork-sizing
    plt.style.use(["science"])
    plt.rcParams.update({"font.size": 8})
    nrows, ncols = 2, 2
    sns.set_style("whitegrid", {"font.family": "serif"})
    fig, axs = plt.subplots(nrows, ncols, height_ratios=[0.1, 0.9], sharey="row")

    # Plot H-ratio and random figures.
    legend_pos = "right"
    plot_data(
        random_data,
        r"$\mathcal{D}_{\text{random}}$",
        axs[1,0],
        axs[0,0],
        annotate_y=True,
        plot_legend="",
    )
    plot_data(
        hratio_data,
        r"$\mathcal{D}_{\text{\texttt{H}-ratio}}$",
        axs[1,1],
        axs[0,1],
        annotate_y=False,
        plot_legend=legend_pos,
    )

    # Style overall figure and save.
    width = 7.48
    fig.set_size_inches(width, nrows * 1.25 + 0.5, forward=True)
    fig.subplots_adjust(wspace=0.35, hspace=0.1)
    fig.tight_layout()
    save_fig(fig, f"fig3_fit_probit_legend_{legend_pos}")


def main():
    plt.style.use(["science"])
    plt.rcParams.update({
        "legend.frameon": True,
        "legend.framealpha": 1,
        "legend.facecolor": "white",
        "legend.edgecolor": "black",
    })
    root = Path(__file__).parent.parent.parent
    os.chdir(root)
    plot_probits()


if __name__ == "__main__":
    main()
