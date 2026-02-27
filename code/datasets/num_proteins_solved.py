#! /usr/bin/env python3
"""
File:           num_proteins_solved.py
Description:    Create a LaTeX table with the number of proteins solved per
                length.
"""

import os
import sys
from pathlib import Path

os.chdir(Path(__file__).parent)
sys.path.insert(0, "..")
from common import (  # noqa: E402
    DATASET_NAME_HRATIO,
    DATASET_NAME_RANDOM,
    load_results,
)

MAX_LENGTH = 30


def main():
    df_hratio = load_results(DATASET_NAME_HRATIO)
    df_random = load_results(DATASET_NAME_RANDOM)
    for df in [df_hratio, df_random]:
        df.drop(df[df["length"] > MAX_LENGTH].index, inplace=True)
    print()
    print("Solved per length (LaTeX table):")
    print("% Columns: Length & D_hratio & D_random")
    for length in list(sorted(df_hratio["length"].unique())) + ["All"]:
        if length == "All":
            n_hratio = len(df_hratio)
            n_random = len(df_random)
            print("\\hdashline")
        else:
            n_hratio = len(df_hratio[df_hratio["length"] == length])
            n_random = len(df_random[df_random["length"] == length])
        print(f"{length} & \\num{{{n_hratio}}} & \\num{{{n_random}}} \\\\")
    print()


if __name__ == "__main__":
    main()
