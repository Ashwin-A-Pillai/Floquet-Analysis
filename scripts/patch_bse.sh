#!/usr/bin/env bash
#
# patch_optics.sh — tweak an optics BSE/Yambo input:
#   • BSENGBlk ? –1
#   • KfnQPdb ? “E <GW0/ndb.QP”
#   • % BEnRange ? $iE | $fE
#   • BEnSteps ? $stepE
#   • % BndsRnXp ? $Polar_band
#   • NGsBlkXp ? $Block_size
#
# Usage: ./patch_optics.sh <input_file> [output_file]
# Requires:
  export iE=5.0
  export fE=20.0
  export i_corband=2
  export f_corband=7
  export stepE=400
  export Polar_band=80
  export Block_size=3000

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <input_file> [output_file]" >&2
  exit 1
fi

infile="$1"
outfile="${2:-${infile}.patched}"

: "${iE:?Please set \$iE}" \
  "${fE:?Please set \$fE}" \
  "${i_corband:?Please set \$i_corband}" \
  "${f_corband:?Please set \$f_corband}"\
  "${stepE:?Please set \$stepE}" \
  "${Polar_band:?Please set \$Polar_band}" \
  "${Block_size:?Please set \$Block_size}"

awk -v ie="$iE" \
    -v fe="$fE" \
    -v ic="$i_corband" \
    -v fc="$f_corband" \
    -v se="$stepE" \
    -v pb="$Polar_band" \
    -v bs="$Block_size" '
{

  # 2) KfnQPdb tweak
  if ($1 ~ /^KfnQPdb=/) {
    print "KfnQPdb= \"E < GW0/ndb.QP\"         # [EXTQP BSK BSS] Database action"
    next
  }

  # 3) % BEnRange block
  if ($1=="%" && $2=="BEnRange") {
    print                             # “% BEnRange” line
    getline                           # old range line
    printf "  %s | %s |         eV    # [BSS] Energy range\n", ie, fe
    next
  }

  # 4) BEnSteps
  if ($1=="BEnSteps=") {
    printf "BEnSteps= %s                    # [BSS] Energy steps\n", se
    next
  }

  # 5) % BndsRnXp block
  if ($1=="%" && $2=="BndsRnXp") {
    print                             # “% BndsRnXp” line
    getline                           # old bands line
    printf "   1 |  %s |                         # [Xp] Polarization function bands\n", pb
    next
  }
  
  # 6) % BSEQptR block: force it to 1 | 1
  if ($1=="%" && $2=="BSEQptR") {
    print
    getline
    printf "  1 | 1 |                             # [BSK] Transferred momenta range\n"
    getline
    print
    next
  }
  
  # 5) % BndsRnXp block
  if ($1=="%" && $2=="BSEBands") {
    print                             # “% BndsRnXp” line
    getline                           # old bands line
    printf "  %s  |  %s |", ic, fc
    next
  }

  # default: passthrough
  print
}
END {
  print ""                      # blank line before extras
  print "SE_ROLEs= \"k eh t\""
  print "SE_CPU= \"1 8 8\"                     # Parallelism over q points only"
  print "BS_nCPU_invert = 0"            
  print "BS_nCPU_diago  = 0"            
}
' "$infile" > "$outfile"

echo "Patched input written to $outfile"
