#! /usr/bin/env python3
"""
File:           check_uniqueness.py
Description:    Verify that all datasets only have unique proteins.
"""

from pathlib import Path
import pandas as pd


def check_uniqueness_recursive(path):
    """
    Check if a given dataset contains duplicate proteins.
    Return the number of duplicates.
    """
    # Iterate over subfolders.
    if path.is_dir():
        return sum(check_uniqueness_recursive(p) for p in path.iterdir())
    # Filter non-CSVs.
    if path.suffix.lower() != ".csv":
        print("No CSV: ", path)
        return 0

    # If is file and CSV, count duplicates and return.
    duplicates = 0
    try:
        df = pd.read_csv(path)
        duplicate_series = df["sequence"].duplicated(keep=False).value_counts()
        if True in duplicate_series.index:
            duplicates = duplicate_series[True]
    except:  # noqa: E722
        print(f"Error: Skipping file {path}")
    return duplicates


if __name__ == "__main__":
    datasets_dir = Path(__file__).parent.parent.parent.joinpath("data")
    for p in datasets_dir.iterdir():
        if not p.is_dir():
            continue
        print(f"{p.name}\t({check_uniqueness_recursive(p)} duplicates)")
