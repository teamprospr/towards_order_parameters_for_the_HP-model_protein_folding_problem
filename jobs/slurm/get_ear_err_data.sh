awk '
BEGIN {
    print "jobid,wall_time_s,dc_power_w,energy_kwh"
}

/job.step:/ {
    if (match($0, /job\.step:[[:space:]]*([0-9]+)/, m))
        jobid = m[1]
}

/Wall time:/ {
    if (match($0, /Wall time:[[:space:]]*([0-9.]+)/, m))
        wall = m[1]
}

/DC\/DRAM\/PCK power:/ {
    if (match($0, /DC\/DRAM\/PCK power:[[:space:]]*([0-9.]+)/, m))
        dc = m[1]
}

jobid && wall && dc {
    kwh = (wall * dc) / 3.6e6
    printf "%s,%.3f,%.3f,%.6f\n", jobid, wall, dc, kwh
    jobid = wall = dc = ""
}
' *.err > energy.csv

