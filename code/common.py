#! /usr/bin/env python3
"""
File:           common.py
Description:    Set of common functions used for loading the experiment results.
"""
import re
from pathlib import Path

import pandas as pd

ROOT_DIR = Path(__file__).resolve().parent.parent
ROOT_DATASTE_DIR = ROOT_DIR / "data"
DATASET_NAME_HRATIO = "vanEck_hratio"
DATASET_NAME_RANDOM = "vanEck_random"


def _load_results_length_csv(
    csv_path: Path,
    path_dataset: Path,
    df: pd.DataFrame,
    length: int,
    h_count: int | None = None,
) -> pd.DataFrame:
    # Only load the file containing the final results
    if re.search(r"_r\d+(?:_tmp)?\.csv$", csv_path.name.lower()):
        return df
    print(f"Loading {csv_path}")
    assert length == int(csv_path.name.split("_")[1])
    if h_count is not None:
        assert h_count == int(csv_path.name.split("_")[2][1:])
    df_h_ratio_results = pd.read_csv(csv_path)
    if df_h_ratio_results.empty:
        return df
    df_h_ratio_dataset = pd.read_csv(path_dataset)
    # assert df_h_ratio_results["protein_id"].is_unique, csv_path
    assert df_h_ratio_dataset["id"].is_unique
    df_h_ratio_results = df_h_ratio_results.merge(
        df_h_ratio_dataset[["id", "sequence"]],
        left_on="protein_id",
        right_on="id",
        how="left",
    )
    df_h_ratio_results = df_h_ratio_results[
        ["protein_id", "placed", "sequence", "hash"]
    ].copy()
    df_h_ratio_results["length"] = length
    df_h_ratio_results["h_ratio"] = df_h_ratio_results["sequence"].apply(
        lambda s: len([c for c in s if c == "H"]) / length
    )
    df_h_ratio_results["path"] = csv_path
    df = pd.concat([df, df_h_ratio_results], ignore_index=True)
    return df


def load_results_length(
    base_path_results: Path,
    base_path_dataset: Path,
    length: int,
    df: pd.DataFrame | None = None,
) -> None:
    for path_results in base_path_results.iterdir():
        if path_results.is_file():
            # Random dataset:
            path_dataset = base_path_dataset.with_suffix(".csv")
            if path_results.name.endswith(".csv"):
                df = _load_results_length_csv(path_results, path_dataset, df, length)
        else:
            # H-ratio dataset:
            if path_results.is_dir() and not path_results.name.startswith("H"):
                continue
            h_count = int(path_results.name[1:])
            path_dataset = base_path_dataset / f"{length}_H{h_count}.csv"
            for csv_path in path_results.glob("*.csv"):
                df = _load_results_length_csv(
                    csv_path, path_dataset, df, length, h_count
                )
    return df


# TODO load "flat results" from reversed sequence experiments
def load_results(dataset_name: str, match_suffix: str = "_results") -> pd.DataFrame:
    df: pd.DataFrame | None = None
    match_prefix = dataset_name.split("_")[-1] + "_"
    for path_results in (ROOT_DIR / Path("jobs")).iterdir():
        if (
            not path_results.is_dir()
            or not path_results.name.startswith(match_prefix)
            or not path_results.name.endswith(match_suffix)
        ):
            continue
        length = -1
        name_parts = path_results.name.split("_")
        # assert name_parts[0] == dataset_name, f"{path_results.name} != {dataset_name}"
        assert name_parts[-1] == "results"
        length = int(name_parts[-2][1:])
        path_dataset = ROOT_DATASTE_DIR / dataset_name / str(length)
        df = load_results_length(path_results, path_dataset, length, df)
    df["H_count"] = df["sequence"].map(lambda s: sum(1 if c == "H" else 0 for c in s))
    return df


def load_results_hardest(
    dataset_name: str,
    fraction=0.05,
    match_suffix: str = "_results",
) -> pd.DataFrame:
    """Load the results and only return the hardest <fraction>."""
    df = load_results(dataset_name, match_suffix=match_suffix)
    df.sort_values("placed", inplace=True)
    result_data = []
    for _, group in df.groupby("length"):
        n = max(1, int(len(group) * fraction))
        top_rows = group.nlargest(n, "placed")
        result_data.append(top_rows)
    return pd.concat(result_data, ignore_index=True)


def load_results_median(
    dataset_name: str,
    fraction=0.05,
    match_suffix: str = "_results",
) -> pd.DataFrame:
    """
    Load the results and only return the <fraction> number of rows around the
    median.
    """
    # Load job results and create a grouped DataFrame per length.
    df = load_results(dataset_name, match_suffix=match_suffix)
    df.sort_values("placed", inplace=True)
    result_rows = []

    for _, group in df.groupby("length"):
        # Sort by hardness.
        g = group.sort_values("placed").reset_index(drop=True)
        n = len(g)

        # Compute number of rows to return.
        window = max(1, int(n * fraction))

        # Centered window around median of values.
        median_idx = n // 2
        half = window // 2
        start = max(0, median_idx - half)
        end = min(n, start + window)

        # If clipped at end, adjust start. Then add rows.
        start = max(0, end - window)
        result_rows.append(g.iloc[start:end])

    return pd.concat(result_rows, ignore_index=True)
