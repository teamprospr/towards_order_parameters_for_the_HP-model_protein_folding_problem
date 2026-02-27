# Get JobIDs in a list.
jobids=($(find . -maxdepth 1 -type f \( -name "*.out" -o -name "*.err" \) \
    | sed -E 's/.*-([0-9]+)\.(out|err)/\1/' \
    | sort -u))

# Setup CSV file.
echo "jobid,joules" >> job_joules.csv

# Process per 10 to not get dumped cores.
for ((i=0; i<${#jobids[@]}; i+=5)); do
    batch=$(IFS=,; echo "${jobids[*]:i:5}")
    echo "=== Batch $((i/5+1)) ==="
    echo -e "\nJobs: $batch"
    eacct -j "$batch" | awk '$1 ~ /-sb$/ {
            split($1, a, "-");
            printf "%s,%s\n", a[1], $13
        }' >> job_joules.csv
done

