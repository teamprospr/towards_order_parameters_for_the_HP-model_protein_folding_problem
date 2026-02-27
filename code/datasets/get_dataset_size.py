#! /usr/bin/env python3
"""
File:           get_dataset_size.py
Description:    Get the size of the generated datasets.
                Set VERBOSE for more info.
"""

from pathlib import Path
import pandas as pd

VERBOSE = True


def count_csv_rows_recursive(path):
    if path.is_dir():
        return sum(count_csv_rows_recursive(p) for p in path.iterdir())
    if path.suffix.lower() != ".csv":
        print("foo", path)
        return 0
    rows = 0
    try:
        df = pd.read_csv(path)
        rows = len(df)
        if VERBOSE:
            print(f"\t{path.name}: {rows}")
    except:  # noqa: E722
        print(f"Error: Skipping file {path}")
    return rows


if __name__ == "__main__":
    datasets_dir = Path(__file__).parent.parent.parent.joinpath("data")
    for p in datasets_dir.iterdir():
        if not p.is_dir():
            continue
        if VERBOSE:
            print(f"{p.name}:")
            total_proteins = count_csv_rows_recursive(p)
            print(f"{p.name}:\t({total_proteins} proteins)\n")
        else:
            print(f"{p.name}\t({count_csv_rows_recursive(p)} proteins)")
