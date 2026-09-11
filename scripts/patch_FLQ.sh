#!/usr/bin/env bash
#
###############################################################################
# Defaults with override behavior
###############################################################################
: "${nl_band_i:=2}"
: "${nl_band_f:=7}"
: "${fl_order:=-4}"
: "${field1_freq:=4.100000}"
: "${field1_kind:=SIN}"
: "${efield1_x:=1.000000}"
: "${efield1_y:=0.000000}"
: "${efield1_z:=0.000000}"

export nl_band_i nl_band_f fl_order \
       field1_freq field1_kind \
       efield1_x efield1_y efield1_z

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <input_file> [output_file]" >&2
  exit 1
fi

infile="$1"
outfile="${2:-${infile}.modified}"

awk -v bi="$nl_band_i" \
    -v bf="$nl_band_f" \
    -v flo="$fl_order" \
    -v ffq="$field1_freq" \
    -v fkind="$field1_kind" \
    -v dx="$efield1_x" \
    -v dy="$efield1_y" \
    -v dz="$efield1_z" '
{
  ###########################################################################
  # 1) NLBands block:
  #    % NLBands
  #      4 | 5 |
  ###########################################################################
  if (($1 == "%" && $2 == "NLBands") || $1 == "%NLBands") {
    print
    if (getline > 0) {
      printf "   %s |  %s |                         # [NL] Bands range\n", bi, bf
    }
    next
  }

  ###########################################################################
  # 2) FLOrder line:
  #    FLOrder=-4
  ###########################################################################
  if ($1 ~ /^FLOrder=/) {
    printf "FLOrder=%s                       # [NL] Fourier Order of Floquet Analysis\n", flo
    next
  }

  ###########################################################################
  # 3) Field1_Freq line:
  #    Field1_Freq= 4.100000      eV
  ###########################################################################
  if ($1 == "Field1_Freq=") {
    printf "Field1_Freq= %s      eV    # [RT Field1] Frequency\n", ffq
    next
  }

  ###########################################################################
  # 4) Field1_kind line:
  #    Field1_kind= "DELTA"
  ###########################################################################
  if ($1 == "Field1_kind=") {
    printf "Field1_kind= \"%s\"             # [RT Field1] Kind(SIN|COS|RES|ANTIRES|GAUSS|DELTA|QSSIN)\n", fkind
    next
  }

  ###########################################################################
  # 5) Field1_Dir block:
  #    % Field1_Dir
  #      1.000000 | 1.000000 | 0.000000 |
  ###########################################################################
  if (($1 == "%" && $2 == "Field1_Dir") || $1 == "%Field1_Dir") {
    print
    if (getline > 0) {
      printf " %s | %s | %s |        # [RT Field1] Versor\n", dx, dy, dz
    }
    next
  }

  # default: pass through unchanged
  print
}
' "$infile" > "$outfile"

echo "Patched input written to $outfile"
echo "  nl_band_i   = $nl_band_i"
echo "  nl_band_f   = $nl_band_f"
echo "  fl_order    = $fl_order"
echo "  field1_freq = $field1_freq"
echo "  field1_kind = $field1_kind"
echo "  field1_dir  = $efield1_x $efield1_y $efield1_z"