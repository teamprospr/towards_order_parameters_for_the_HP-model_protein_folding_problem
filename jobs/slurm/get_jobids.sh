jobids=$(find . -maxdepth 1 -type f \( -name "*.out" -o -name "*.err" \) \
  | sed -E 's/.*-([0-9]+)\.(out|err)/\1/' \
  | sort -u \
  | paste -sd, -)
echo "$jobids"
