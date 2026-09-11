#!/usr/bin/env bash
#
###############################################################################
# Defaults with override behavior
###############################################################################
: "${Polar_band:=40}"
: "${Block_size:=6000}"
: "${k_final:=1}"
: "${i_corband:=3}"
: "${f_corband:=8}"

export Polar_band Block_size k_final i_corband f_corband

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <input_file> [output_file]" >&2
  exit 1
fi

infile="$1"
outfile="${2:-${infile}.modified}"

awk -v pb="$Polar_band" \
    -v bs="$Block_size" \
    -v kf="$k_final" \
    -v ic="$i_corband" \
    -v fc="$f_corband" '
{
  ###########################################################################
  # 1) Polarization function bands block
  #    % BndsRnXp
  #      1 |   80 |
  ###########################################################################
  if ($1 == "%" && $2 == "BndsRnXp") {
    print                     # the "% BndsRnXp" line
    if (getline > 0) {        # consume and replace the next line
      printf "   1 |  %s |                         # [Xp] Polarization function bands\n", pb
    }
    next
  }

  ###########################################################################
  # 2) Response block size:
  #    NGsBlkXp=  <Block_size>  mHa
  ###########################################################################
  if ($1 == "NGsBlkXp=") {
    printf "NGsBlkXp= %s               mHa    # [Xp] Response block size\n", bs
    next
  }

  ###########################################################################
  # 3) QPkrange block:
  #    %QPkrange
  #      1 |  k_final |  i_corband |  f_corband |
  ###########################################################################
  if ($1 == "%QPkrange" || ($1 == "%" && $2 == "QPkrange")) {
    print                     # the "%QPkrange" line
    if (getline > 0) {        # read and discard the old numeric line
      printf "  1 | %s | %s | %s |   \n", kf, ic, fc
    }
    next
  }

  # default: just pass through
  print
}
END {
  print ""                      # blank line before extras
  print "SE_ROLEs= \"q qp b\""
  print "SE_CPU= \"1 4 4\"                     # Parallelism over q points only"
}
' "$infile" > "$outfile"

echo "Patched input written to $outfile"
echo "  Polar_band = $Polar_band"
echo "  Block_size = $Block_size"
echo "  k_final    = $k_final"
echo "  i_corband  = $i_corband"
echo "  f_corband  = $f_corband"