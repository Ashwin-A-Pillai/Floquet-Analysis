#!/usr/bin/env bash
#
# patch_ypp.sh — modify a YPP bands input:
#   • INTERP_mode -> BOLTZ
#   • % BANDS_bands -> [1 | $Num_bands]
#   • BANDS_steps  -> $Num_steps
#   • GfnQPdb      -> by $YPP_MODE (KS -> "none", QP -> "E < GW0/ndb.QP")
#   • append a %BANDS_kpts block
#
# Usage: ./patch_ypp.sh <input_file> [output_file]
#

set -euo pipefail

# User-tunable (can be overridden from environment)
export Num_bands=${Num_bands:-8}
export Num_steps=${Num_steps:-40}
export YPP_MODE=${YPP_MODE:-KS}     # KS | QP

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <input_file> [output_file]" >&2
  exit 1
fi

infile="$1"
outfile="${2:-${infile}.modified}"

: "${Num_bands:?}" "${Num_steps:?}" "${YPP_MODE:?}"

awk -v nb="$Num_bands" -v ns="$Num_steps" '
BEGIN{
  mode = ENVIRON["YPP_MODE"];   # KS or QP
}
{
  # Keep comments as-is
  if ($0 ~ /^[[:space:]]*#/) { print; next }

  # 0) Drop any existing %BANDS_kpts … % block
  if ($1=="%BANDS_kpts") {
    while (getline) {
      if ($1=="%") break
    }
    next
  }

  # 1) INTERP_mode -> BOLTZ
  if ($0 ~ /^[[:space:]]*INTERP_mode[[:space:]]*=/) {
    print "INTERP_mode= \"BOLTZ\"                # Interpolation mode (NN/BOLTZ)"
    next
  }

  # 2) % BANDS_bands
  if ($1=="%" && $2=="BANDS_bands") {
    print                 # header line
    getline               # old "1 | N |"
    printf "   1 |  %s |                         # Number of bands\n", nb
    next
  }

  # 3) BANDS_steps
  if ($0 ~ /^[[:space:]]*BANDS_steps[[:space:]]*=/) {
    printf "BANDS_steps= %s                  # Number of divisions\n", ns
    next
  }

  # 4) GfnQPdb depends on mode
  if ($0 ~ /^[[:space:]]*GfnQPdb[[:space:]]*=/) {
    if (mode=="KS") {
      print "GfnQPdb= \"none\"                  # [EXTQP G] Database action"
    } else {
      print "GfnQPdb= \"E < GW0/ndb.QP\"        # [EXTQP G] Database action"
    }
    next
  }

  # default: passthrough
  print
}
END{
  print ""
  print "%BANDS_kpts                      # hBN path: Gamma-M-K-Gamma"
  print " 0.0       | 0.0       | 0.0      |   # Gamma"
  print " 0.5       | 0.0       | 0.0      |   # M"
  print " 0.333333  | 0.333333  | 0.0      |   # K"
  print " 0.0       | 0.0       | 0.0      |   # Gamma"
  print "%"
}
' "$infile" > "$outfile"

echo "Patched input written to $outfile"