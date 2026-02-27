#!/bin/bash

[ -z "$VERBOSE" ] && VERBOSE=0

# Check if directory is given
if [ -z "$1" ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

# Store root of results.
rdir="$1"

# Check for any *bnb.csv in the result directory
shopt -s nullglob
files=("$rdir"/*bnb.csv)
shopt -u nullglob

# Report if final file is not there.
if [ ${#files[@]} -eq 0 ]; then
    echo -e "$rdir:\tUnfinished\n"
else
    # Setup paths to the files.
    res_file=${files[0]##*/}
    length="${res_file//[^0-9]/}"
    input_fpath="../data/vanEck_random/${length}.csv"

    # Extract IDs from files (first column).
    expected_ids=$(cut -d',' -f1 ${input_fpath} | tail -n +2 | sort)
    result_ids=$(cut -d',' -f1 ${files[0]} | tail -n +2 | sort)
    
    # Compute the missing and duplicate IDs.
    missing=$(comm -23 <(echo "$expected_ids") <(echo "$result_ids") | sort -n | tr '\n' ' ')
    duplicates=$(echo "$result_ids" | uniq -d | sort -n | tr '\n' ' ')
    
    # Print if there are missing or duplicates.
    if [[ -n "$missing" ]] || [[ -n "$duplicates" ]]; then
        if [ -n "$missing" ]; then
            echo -e "$rdir:\tMissing: $missing"
            really_missing=""

            # Search for missing IDs in temprorary files.
            for id in $missing; do
                matches=$(grep -l "^$id," "$rdir"/*_r*.csv)
                if [[ -n "$matches" ]]; then
                    if [[ $VERBOSE -eq 1 ]]; then
                        echo -e "\t$id FOUND in: $matches"
                    fi
                else
                    # Construct new list of missing IDs if not found.
                    really_missing="$really_missing $id"
                fi
            done
            
            # Report if all missing were found, otherwise new list.
            if [ -z "$really_missing" ]; then
                echo -e "$rdir:\tAll missing found."
            else
                echo -e "$rdir:\tReally missing: $really_missing"
            fi
        fi
        
        [[ -n "$duplicates" ]] && echo -e "$rdir: Duplicates: $duplicates"
        echo ""
    else
        echo -e "$rdir:\tFinished succesfully!\n"
    fi
fi

