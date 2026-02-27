#! /usr/bin/env python3
"""
File:           generate_random_dataset.py
Description:    This file generates data for testing using random sampling.
"""

import random
from math import factorial
from pathlib import Path

import pandas as pd


def generate_random(p_len=10, size=1000):
    """
    Generates a new random dataset where the number of proteins is equal to the
    number of proteins that was generated in the refined H-ratio dataset.
    However, contrary to the H-ratio dataset, there is no filtering on
    uniqueness.
    """
    aminos = ["H", "P"]

    # Compute the number of proteins to generate from the H-ratios.
    num_proteins = 0

    for H_count in range(p_len + 1):
        # Compute how many unique proteins exist and limit to 1000.
        P_count = p_len - H_count
        unique_proteins = int(
            factorial(p_len) / (factorial(H_count) * factorial(P_count))
        )
        num_proteins += min(1000, unique_proteins)
    print(f"Generating {num_proteins} proteins for length {p_len}..")

    # Generate the number of proteins by permutating the protein sequence.
    new_proteins = []
    proteins = set()

    while len(proteins) < num_proteins:
        new_proteins = ["".join(random.choices(aminos, k=p_len)) for _ in range(500)]
        proteins.update(set(new_proteins))

    # Remove access proteins if set contains too many.
    while len(proteins) > num_proteins:
        proteins.pop()

    # Write proteins to a file.
    res_file = Path(f"{ds_path}/{p_len}.csv")
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
    ds_path = Path(__file__).parent.parent.parent.joinpath("data/vanEck_random")
    Path(ds_path).mkdir(exist_ok=True)

    # Set the random generator's seed to a fixed date for reproducibility.
    generation_epoch = 1735755521.6005082
    random.seed(generation_epoch)

    # Generate the data using the provided arguments.
    generate_random(p_len, size)
