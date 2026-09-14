#!/usr/bin/env bash
#
###############################################################################
# Defaults with override behavior
###############################################################################
: "${efield1_x:=1.000000}"
: "${efield1_y:=0.000000}"
: "${efield1_z:=0.000000}"
: "${rm_time_rev:=yes}"

export efield1_x efield1_y efield1_z rm_time_rev

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <input_file> [output_file]" >&2
  exit 1
fi

infile="$1"
outfile="${2:-${infile}.modified}"

awk -v ex="$efield1_x" \
    -v ey="$efield1_y" \
    -v ez="$efield1_z" \
    -v rmtr="$rm_time_rev" '
{
  ###########################################################################
  # 1) Efield1 block:
  #    % Efield1
  #      1.000000 | 1.000000 | 0.000000 |
  ###########################################################################
  if (($1 == "%" && $2 == "Efield1") || $1 == "%Efield1") {
    print
    if (getline > 0) {
      printf " %s | %s | %s |        # First external Electric Field\n", ex, ey, ez
    }
    next
  }

  ###########################################################################
  # 2) Uncomment RmTimeRev if requested:
  #    #RmTimeRev  ->  RmTimeRev
  ###########################################################################
  if ($0 ~ /^[[:space:]]*#[[:space:]]*RmTimeRev/ && rmtr == "yes") {
    sub(/^[[:space:]]*#[[:space:]]*/, "")
    print
    next
  }

  # default: pass through unchanged
  print
}
' "$infile" > "$outfile"

###############################################################################
# Safety checks
###############################################################################
grep -q '1\.000000[[:space:]]*\|[[:space:]]*1\.000000[[:space:]]*\|[[:space:]]*0\.000000[[:space:]]*\|[[:space:]]*#[[:space:]]*First external Electric Field' "$outfile"

if [[ "$rm_time_rev" == "yes" ]]; then
  grep -q '^[[:space:]]*RmTimeRev' "$outfile"
fi

echo "Patched input written to $outfile"
echo "  efield1_x   = $efield1_x"
echo "  efield1_y   = $efield1_y"
echo "  efield1_z   = $efield1_z"
echo "  rm_time_rev = $rm_time_rev"