#!/usr/bin/env bash
#
###############################################################################
# Defaults with override behavior
###############################################################################
: "${en_range_min:=0.00000}"
: "${en_range_max:=10.00000}"
: "${damp_mode:=LORENTZIAN}"
: "${damp_factor:=0.100000}"

export en_range_min en_range_max damp_mode damp_factor

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <input_file> [output_file]" >&2
  exit 1
fi

infile="$1"
outfile="${2:-${infile}.modified}"

awk -v emin="$en_range_min" \
    -v emax="$en_range_max" \
    -v dmode="$damp_mode" \
    -v dfac="$damp_factor" '
{
  ###########################################################################
  # 1) EnRngeRt block:
  #    % EnRngeRt
  #      0.00000 | 20.00000 |
  ###########################################################################
  if (($1 == "%" && $2 == "EnRngeRt") || $1 == "%EnRngeRt") {
    print
    if (getline > 0) {
      printf "  %s | %s |         eV    # Energy range\n", emin, emax
    }
    next
  }

  ###########################################################################
  # 2) DampMode line:
  #    DampMode= "NONE"
  ###########################################################################
  if ($1 == "DampMode=") {
    printf "DampMode= \"%s\"                 # Damping type ( NONE | LORENTZIAN | GAUSSIAN )\n", dmode
    next
  }

  ###########################################################################
  # 3) DampFactor line:
  #    DampFactor= 0.000000       eV
  ###########################################################################
  if ($1 == "DampFactor=") {
    printf "DampFactor= %s       eV    # Damping parameter\n", dfac
    next
  }

  # default: pass through unchanged
  print
}
' "$infile" > "$outfile"

echo "Patched input written to $outfile"
echo "  en_range_min = $en_range_min"
echo "  en_range_max = $en_range_max"
echo "  damp_mode    = $damp_mode"
echo "  damp_factor  = $damp_factor"