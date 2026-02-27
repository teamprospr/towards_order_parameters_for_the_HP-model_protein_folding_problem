# Get JobIDs in a list.
jobids=($(find . -maxdepth 1 -type f \( -name "*.out" -o -name "*.err" \) \
    | sed -E 's/.*-([0-9]+)\.(out|err)/\1/' \
    | sort -u))

# Process per 10 to not get dumped cores.
for ((i=0; i<${#jobids[@]}; i+=5)); do
    batch=$(IFS=,; echo "${jobids[*]:i:5}")
    echo "=== Batch $((i/5+1)) ==="
    echo -e "\nJobs: $batch"
    eacct -j "$batch" | awk '$1 ~ /-sb$/ { sum += $13 } END { print sum }'
done

