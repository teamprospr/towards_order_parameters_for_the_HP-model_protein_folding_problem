jobids=$(./get_jobids_sb.sh)

sacct -j "$jobids" --format=Elapsed,AllocCPUS --noheader --parsable2 | \
    awk -F'|' '
    function elapsed_to_hours(s) {
        # Check if there is a day prefix
        if (s ~ /-/) {
            split(s, dpart, "-");
            days = dpart[1];
            t = dpart[2];
        } else {
            days = 0;
            t = s;
        }
        split(t, a, ":");
        # HH:MM:SS -> hours
        if (length(a) == 3) return days*24 + a[1] + a[2]/60 + a[3]/3600;
        # MM:SS -> hours
        else if (length(a) == 2) return days*24 + a[1]/60 + a[2]/3600;
        else return days*24 + a[1];  # fallback
    }
    {
        hrs = elapsed_to_hours($1);
        total_hours += hrs;
        cpu_hours += hrs * $2;
    }
    END {
        printf "Total runtime (hours): %.2f\n", total_hours;
        printf "Total CPU core-hours: %.2f\n", cpu_hours;
    }'

