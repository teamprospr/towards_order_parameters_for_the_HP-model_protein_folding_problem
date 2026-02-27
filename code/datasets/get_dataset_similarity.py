#! /usr/bin/env python3
"""
File:           get_dataset_similarity.py
Description:    Check for similarity between the two generated datasets.
"""

from pathlib import Path
import pandas as pd


def get_dataset_sequences(dataset_path: Path, max_length: int=100):
    """Get dataset sizes with a limiter on length."""
    proteins = set()
    for p in dataset_path.iterdir():
        if p.is_file() and p.suffix == ".csv" \
        and int(p.stem.split("_")[0]) <= max_length:
            # Get list of protein sequences.
            df = pd.read_csv(p)
            proteins.update(df["sequence"].unique())
        if p.is_dir() and int(p.name) <= max_length:
            # Process dir recursively.
            proteins.update(get_dataset_sequences(p))
    return proteins


def main():
    """Entrypoint."""
    # Setup paths to the datasets.
    datasets_dir = Path(__file__).parent.parent.parent.joinpath("data")
    hratio_path = Path(datasets_dir, "vanEck_hratio")
    random_path = Path(datasets_dir, "vanEck_random")

    # Limit count to length 30.
    max_lenght = 30

    # Store protein sequences as set of strings.
    proteins_hratio = get_dataset_sequences(hratio_path, max_lenght)
    proteins_random = get_dataset_sequences(random_path, max_lenght)

    # Compute stats.
    len_hratio = len(proteins_hratio)
    len_random = len(proteins_random)
    total_size = len_hratio + len_random
    num_unique = len(proteins_hratio.union(proteins_random))
    num_double = len(proteins_hratio.intersection(proteins_random))

    # Report stats.
    print("Length H-ratio: ", len_hratio)
    print("Length random:  ", len_random)
    print("============================")
    print("Total proteins:     ", len_hratio + len_random)
    print("Recurring proteins: ", num_double)
    print("Unique proteins:    ", num_unique)
    print("============================")
    print(f"Percentage overlap: {num_double / total_size * 100:.2f} %")
    print(f"Percentage unique:  {num_unique / total_size * 100:.2f} %")


if __name__ == "__main__":
    main()
