# ------------------------------------------------------------------------------
# Script for merging results that are fully finished.
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


def get_num_proteins(resdir: str):
    """Get number of proteins in input file based on results directory."""
    length, hratio = get_length_hratio_from_resdir(resdir)
    if hratio:
        data_path = f"../data/vanEck_hratio/{length}/{length}_{hratio}.csv"
    else:
        data_path = f"../data/vanEck_random/{length}.csv"

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


def check_merged_results(resdir: str, merged_file: str):
    """Check if all IDs are present in a merged CSV file."""
    # Get maximum ID for verification.
    num_proteins = get_num_proteins(resdir)

    # Initialize a dictionary to store {protein_id: hash} pairs.
    protein_id_hash = {}
    all_data_lines = []
    duplicate_ids = []

    with open(merged_file, "r") as fp:
        # Setup CSV reader and skip header.
        reader = csv.reader(fp)
        next(reader)

        # Loop over all IDs and keep track of duplicates and missing IDs.
        for row in reader:
            protein_id = row[0]
            current_hash = row[-1]

            # Check for duplicates.
            if protein_id in protein_id_hash.keys():
                duplicate_ids.append(protein_id)
                if protein_id_hash[protein_id] != current_hash:
                    print(
                        f"\t\tWARNING: Different hashes found for "
                        + f"{protein_id}!"
                    )
            else:
                # Store protein for lookups and CSV writing.
                protein_id_hash[protein_id] = current_hash
                all_data_lines.append(row)

    # Report duplicate IDs.
    if duplicate_ids:
        print(f"\tDuplicate IDs:\n\t\t{duplicate_ids}")

    # Get all protein IDs and check for missing ones.
    protein_ids = sorted(protein_id_hash.keys(), key=int)
    missing_ids = [
        str(i) for i in range(1, num_proteins) if str(i) not in protein_ids
    ]

    # Report if no problems found, or print empty line for readability.
    if not duplicate_ids and not missing_ids:
        print(f"\tAlready merged.")
    
    # Report if protein IDs are missing, otherwise create merged file.
    if missing_ids:
        print(f"\tMissing IDs:\n\t\t{missing_ids}")
    elif duplicate_ids:
        print("\tRemoved duplicates from merged results.")

        # Sort all_data_lines by protein_id.
        all_data_lines.sort(key=lambda x: int(x[0]))

        # Write merged results to file.
        with open(merged_file, mode="w", newline="") as fp:
            writer = csv.writer(fp)
            writer.writerow(
                [
                    "protein_id",
                    "algorithm",
                    "time",
                    "score",
                    "checked",
                    "placed",
                    "hash",
                ]
            )
            writer.writerows(all_data_lines)

    # Print writeline for keeping output tidy.
    print("")


def gather_partial_results(resdir: str, outfile: str):
    """
    Gather partial results from rank files and check for duplicates and missing
    proteins.
    """
    # Get maximum ID for verification.
    num_proteins = get_num_proteins(resdir)

    # Initialize a dictionary to store {protein_id: hash} pairs.
    protein_id_hash = {}
    all_data_lines = []
    duplicate_ids = []

    # Loop over all files in the results folder
    for filename in os.listdir(resdir):
        if filename.startswith("HP_") and filename.endswith(".csv"):
            filepath = os.path.join(resdir, filename)
            with open(filepath, mode="r") as csvfile:
                # Setup CSV reader and skip header.
                reader = csv.reader(csvfile)
                next(reader)

                for row in reader:
                    protein_id = row[0]
                    current_hash = row[-1]

                    # Check for duplicates
                    if protein_id in protein_id_hash.keys():
                        duplicate_ids.append(protein_id)
                        if protein_id_hash[protein_id] != current_hash:
                            print(
                                f"\t\tWARNING: Different hashes found for "
                                + f"{protein_id}!"
                            )
                    else:
                        # Store protein for lookups and CSV writing.
                        protein_id_hash[protein_id] = current_hash
                        all_data_lines.append(row)

    # Report duplicate IDs.
    if duplicate_ids:
        print(f"\tDuplicate IDs:\n\t\t{duplicate_ids}")

    # Get all protein IDs and check for missing ones.
    protein_ids = sorted(protein_id_hash.keys(), key=int)
    missing_ids = [
        str(i) for i in range(1, num_proteins) if str(i) not in protein_ids
    ]

    # Report if protein IDs are missing, otherwise create merged file.
    if missing_ids:
        print(f"\tMissing IDs:\n\t\t{missing_ids}")
    else:
        print("\tAll IDs present, writing merged results.")

        # Sort all_data_lines by protein_id.
        all_data_lines.sort(key=lambda x: int(x[0]))

        # Write merged results to file.
        with open(outfile, mode="w", newline="") as fp:
            writer = csv.writer(fp)
            writer.writerow(
                [
                    "protein_id",
                    "algorithm",
                    "time",
                    "score",
                    "checked",
                    "placed",
                    "hash",
                ]
            )
            writer.writerows(all_data_lines)

    # Print writeline for keeping output tidy.
    print("")


def merge_results(resdir: str):
    """Check if all protein IDs are present, and merge results if so."""
    # Get merged file name.
    length, hratio = get_length_hratio_from_resdir(resdir)
    if hratio:
        merged_file = f"{resdir}/HP_{length}_{hratio}_dfs_bnb.csv"
    else:
        merged_file = f"{resdir}/HP_{length}_dfs_bnb.csv"

    # For this results dir, check its merged file or combine partial results.
    if os.path.isfile(merged_file):
        print(f"{resdir} (merged):")
        check_merged_results(resdir, merged_file)
    else:
        print(f"{resdir}:")
        gather_partial_results(resdir, merged_file)


def main():
    """Entrypoint for execution to scope variables locally."""
    # Get results directories from user provided main results folder.
    if len(sys.argv) < 2:
        print("No results directory provided.")
        print("Usage:\tpython merge_results.py <results_directory>")
    input_dir = sys.argv[1]
    res_dirs = parse_input_dir(input_dir)

    # Merge result subdirectories, if all IDs are present.
    for resdir in res_dirs:
        merge_results(resdir)


if __name__ == "__main__":
    main()
