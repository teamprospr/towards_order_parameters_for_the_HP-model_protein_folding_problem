#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


def get_dataset_sequences(dataset_path: Path, max_length: int=100) -> set:
    """Get dataset sizes with a limiter on length."""
    proteins = set()
    for p in dataset_path.iterdir():
        if p.is_file() and p.suffix == ".csv" \
        and int(p.stem.split("_")[0]) <= max_length:
            # Get list of protein sequences.
            df = pd.read_csv(p)
            proteins.update(df["sequence"])
        if p.is_dir() and int(p.name) <= max_length:
            # Process dir recursively.
            proteins.update(get_dataset_sequences(p))
    return proteins


def report_double_proteins(proteins_double: set):
    """Report the sequence statistics of the double proteins."""
    # Transform the set into a DataFrame and compute H-ratio.
    df = pd.DataFrame(proteins_double, columns=["sequence"])
    df["length"] = df["sequence"].str.len()
    df["hratio"] = df["sequence"].apply(
        lambda s: len([c for c in s if c == "H"])
    ) / df["length"]

    # Compute statistics over set of unique proteins.
    print("Location of double proteins:")
    lengths = sorted(df["length"].unique())

    for length in lengths:
        print(f"    Length {length}:")
        df_length = df[df.length == length]
        total_proteins = 0
        hratios = sorted(df_length["hratio"].unique())

        for hratio in hratios:
            num_proteins = df_length[df_length.hratio == hratio]["sequence"].count()
            total_proteins += num_proteins
            print(f"        H-ratio {hratio}:\t{num_proteins}")

        print(f"        Total:\t{total_proteins}")

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
    proteins_double = proteins_hratio.intersection(proteins_random)
    num_double = len(proteins_double)

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
    print("============================")

    # Report where the double proteins reside.
    report_double_proteins(proteins_double)


if __name__ == "__main__":
    main()
