#! /usr/bin/env python3
"""
File:           generate_hratio_dataset.py
Description:    This file generates data for testing different H-ratios.
"""

import random
from math import factorial
from pathlib import Path

import pandas as pd


def generate_hratio(p_len=10, size=1000):
    """
    Generates a refind version of the H-ratio dataset given the provided
    arguments. The script creates a bin per possible H-ratio, instead of using
    the intervals from the previous paper. For example, a length 10 protein will
    have the bins: [0/10, 1/10, 2/10, .., 10/10]. This way of generating the
    proteins increases the diversity of the dataset. Per bin we generate up to
    1000 proteins to represent the different H-ratio values equally.
    """
    aminos = ["H", "P"]

    # Abort if there is already data on the given length.
    len_path = f"{ds_path}/{p_len}"
    if Path(len_path).exists():
        print(f"Dataset for length {p_len} already exists.")
        return
    Path(len_path).mkdir(exist_ok=True)

    # Generate proteins per H-ratio value.
    for H_count in range(p_len + 1):
        # Compute how many unique proteins exist and limit to 1000.
        P_count = p_len - H_count
        unique_proteins = int(
            factorial(p_len) / (factorial(H_count) * factorial(P_count))
        )
        num_proteins = min(1000, unique_proteins)
        print(f"H-ratio {H_count}/{p_len}:  {num_proteins}")

        # Generate the number of proteins by permutating the protein sequence.
        proteins = set()

        while len(proteins) != num_proteins:
            proteins.add(
                "".join(
                    random.sample(
                        aminos,
                        counts=[H_count, p_len - H_count],
                        k=p_len,
                    )
                )
            )

        # Write proteins to a file.
        res_file = Path(f"{len_path}/{p_len}_H{H_count}.csv")
        with open(res_file, "w") as fp:
            fp.write("id,sequence\n")
        pd.DataFrame(proteins).to_csv(res_file, header=False, mode="a")


if __name__ == "__main__":
    # Fetch arguments from user required for generating data.
    print("Requesting needed information.\nLeave blank for the default value.")
    print("==================================", end="\n\n")

    p_len_str = input("Protein length (default=10): ").strip()
    p_len = int(p_len_str) if p_len_str else 10

    size_str = input("Maximum #proteins per H-ratio (default=1000): ").strip()
    size = int(size_str) if size_str else 1000

    # Create dataset folder if does not exists already.
    ds_path = Path(__file__).parent.parent.parent.joinpath("data/vanEck_hratio")
    Path(ds_path).mkdir(exist_ok=True)

    # Set the random generator's seed to a fixed date for reproducibility.
    generation_epoch = 1735755521.6005082
    random.seed(generation_epoch)

    # Generate the data using the provided arguments.
    generate_hratio(p_len, size)
