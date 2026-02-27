# ------------------------------------------------------------------------------
# Script for constructing a CSV with missing proteins.
# Can be used as input for the solve_remaining experiment.
# ------------------------------------------------------------------------------

import os
import sys
import csv


def parse_input_dir(input_dir: str) -> list[str]:
    """Parse input dir into set of folders to potentially merge."""
    # Verify directory exists.
    if not os.path.isdir(input_dir):
        print("Given results directory does not exist.")
        exit(1)

    # Parse random results by simply pasing only the input directory.
    if input_dir.startswith("random_l"):
        return [input_dir]
    # Otherwise, return the subfolders if hratio results.
    elif input_dir.startswith("hratio_l"):
        res_dirs = []
        for root, dirs, _ in os.walk(input_dir):
            if dirs:
                res_dirs.extend([f"{root}/{dir}" for dir in dirs])
        return res_dirs
    # Error on any other directory name pattern.
    else:
        print("Provided results directory is not random or hratio.")
        exit(1)


def get_length_hratio_from_resdir(resdir: str):
    """Get length and hratio from the resdir."""
    # Get path to input data file.
    length = resdir.split("l")[1].split("_")[0]
    if resdir.startswith("random_l"):
        return length, ""
    elif resdir.startswith("hratio_l"):
        hratio = resdir.split("/")[1]
        return length, hratio
    else:
        print("Result directory does not start with 'random' or 'hratio'.")
        exit(1)


def get_data_path(length: str, hratio: str):
    """Construct dataset path based on given length and hratio."""
    if hratio:
        return f"../data/vanEck_hratio/{length}/{length}_{hratio}.csv"
    else:
        return f"../data/vanEck_random/{length}.csv"


def get_num_proteins(resdir: str):
    """Get number of proteins in input file based on results directory."""
    data_path = get_data_path(*get_length_hratio_from_resdir(resdir))

    # Read last line with highest protein ID.
    with open(data_path, "rb") as file:
        # Move cursor to before last line break.
        file.seek(-2, os.SEEK_END)
        # Move cursor till after the line break before.
        while file.read(1) != b"\n":
            file.seek(-2, os.SEEK_CUR)
        # Read the last line.
        last_protein = file.readline().decode()

    # Get protein ID from last line.
    return int(last_protein.split(",")[0]) + 1


def get_sequences_from_ids(length: str, hratio: str, missing_ids: list):
    """Construct lookup dict of protein sequences from list of protein ids."""
    data_path = get_data_path(length, hratio)
    wanted = set(missing_ids)
    lookup = {}

    with open(data_path) as fp:
        reader = csv.reader(fp)
        for id_str, seq in reader:
            if id_str in wanted:
                lookup[id_str] = seq

    return lookup


def get_proteins_from_file(resfile: str, hratio: str, protein_id_hash: dict):
    """
    Expand or construct a dictionary mapping {(H-ratio, protein_id): hash}.
    """
    with open(resfile, "r") as fp:
        # Setup CSV reader and skip header.
        reader = csv.reader(fp)
        next(reader)

        # Loop over all IDs and keep track of duplicates and missing IDs.
        for row in reader:
            protein_id = row[0]
            current_hash = row[-1]
            cur_key = (hratio, protein_id)

            # Store protein and report if two records for the same ID differ.
            if cur_key not in protein_id_hash.keys():
                protein_id_hash[cur_key] = current_hash
            elif protein_id_hash[cur_key] != current_hash:
                print(f"{resfile}:\tDifferent hashes for {cur_key}!")
                print(f"\tCurrent:  {current_hash}")
                print(f"\tStored:   {protein_id_hash[cur_key]}")
    
    return protein_id_hash


def get_missing_proteins(resdir: str):
    """Check if all protein IDs are present, and merge results if so."""
    # Get merged file name.
    length, hratio = get_length_hratio_from_resdir(resdir)
    if hratio:
        merged_file = f"{resdir}/HP_{length}_{hratio}_dfs_bnb.csv"
    else:
        merged_file = f"{resdir}/HP_{length}_dfs_bnb.csv"

    # Get missing proteins from either the merged file or partial results.
    if os.path.isfile(merged_file):
        protein_id_hash = get_proteins_from_file(merged_file, hratio, {})
    else:
        # Loop over all partial results and construct one dict of IDs.
        protein_id_hash = {}
        for filename in os.listdir(resdir):
            if filename.startswith("HP_") and filename.endswith(".csv"):
                filepath = os.path.join(resdir, filename)
                protein_id_hash = get_proteins_from_file(
                    filepath, hratio, protein_id_hash
                )

    # Get missing protein IDs for this results directory.
    num_proteins = get_num_proteins(resdir)
    protein_ids = sorted([key[1] for key in protein_id_hash.keys()], key=int)
    missing_ids = [
        str(i) for i in range(1, num_proteins) if str(i) not in protein_ids
    ]

    # Report if protein IDs are missing and return rows to be written to file.
    if missing_ids:
        print(f"{resdir}:\tMissing IDs:\t{missing_ids}")
        sequences = get_sequences_from_ids(length, hratio, missing_ids)
        return [[length, hratio, id, sequences[id]] for id in missing_ids]
    else:
        return []


def main():
    """Entrypoint for execution to scope variables locally."""
    # Get results directories from user provided main results folder.
    if len(sys.argv) < 2:
        print("No results directory provided.")
        print("Usage:\tpython merge_results.py <results_directory>")
    input_dir = sys.argv[1]
    res_dirs = parse_input_dir(input_dir)

    # Merge result subdirectories, if all IDs are present.
    missing_rows = []
    for resdir in res_dirs:
        missing_rows.extend(get_missing_proteins(resdir))

    # Get missing filename from user input or set default.
    length = input_dir.split("l")[1].split("_")[0]
    hratio = True if input_dir.startswith("hratio") else False
    if len(sys.argv) >= 3:
        missing_file = sys.argv[2]
    elif hratio:
        missing_file = f"missing_hratio_l{length}.csv"
    else:
        missing_file = f"missing_random_l{length}.csv"
    
    # Write header and missing IDs to output file.
    if missing_rows:
        with open(missing_file, "w", newline="") as fp:
            csvwriter = csv.writer(fp, lineterminator="\n")
            csvwriter.writerow(["length", "hratio", "protein_id", "sequence"])
            csvwriter.writerows(missing_rows)
    else:
        print(f"{input_dir}:\tNo missing IDs found.")


if __name__ == "__main__":
    main()
