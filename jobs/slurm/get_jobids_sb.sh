jobids=$(find . -maxdepth 1 -type f \( -name "*.out" -o -name "*.out" \) \
  | sed -E 's/.*-([0-9]+)\.(out|err)/\1-sd/' \
  | sort -u \
  | paste -sd, -)
echo "$jobids"
