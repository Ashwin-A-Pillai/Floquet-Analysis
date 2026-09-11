#!/usr/bin/env bash
#
###############################################################################
# Defaults with override behavior
###############################################################################

: "${iq:=1}"
: "${fq:=1}"

: "${nl_band_i:=2}"
: "${nl_band_f:=7}"

: "${nl_en_min:=2.000000}"
: "${nl_en_max:=8.000000}"
: "${nl_esteps:=200}"

: "${efield1_x:=1.000000}"
: "${efield1_y:=0.000000}"
: "${efield1_z:=0.000000}"

export iq fq \
       nl_band_i nl_band_f \
       nl_en_min nl_en_max nl_esteps \
       efield1_x efield1_y efield1_z

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <input_file> [output_file]" >&2
  exit 1
fi

infile="$1"
outfile="${2:-${infile}.modified}"

if [[ ! -f "$infile" ]]; then
  echo "[ERROR] Input file not found: $infile" >&2
  exit 2
fi

###############################################################################
# Patch only the requested entries.
#
# Everything else in the original Yambo/Lumen input is preserved unchanged.
###############################################################################

awk -v iqv="$iq" \
    -v fqv="$fq" \
    -v bi="$nl_band_i" \
    -v bf="$nl_band_f" \
    -v emin="$nl_en_min" \
    -v emax="$nl_en_max" \
    -v esteps="$nl_esteps" \
    -v ex="$efield1_x" \
    -v ey="$efield1_y" \
    -v ez="$efield1_z" '
{
  ###########################################################################
  # 1) QpntsRXd block
  #
  #    % QpntsRXd
  #      1 | 1 |
  ###########################################################################
  if (($1 == "%" && $2 == "QpntsRXd") || $1 == "%QpntsRXd") {
    print

    if (getline > 0) {
      printf "   %s | %s |                         # [Xd] Transferred momenta\n", \
             iqv, fqv
    }

    next
  }

  ###########################################################################
  # 2) BndsRnXd block
  #
  #    % BndsRnXd
  #      1 | 50 |
  ###########################################################################
  if (($1 == "%" && $2 == "BndsRnXd") || $1 == "%BndsRnXd") {
    print

    if (getline > 0) {
      printf "  %s | %s |                           # [Xd] Polarization function bands\n", \
             bi, bf
    }

    next
  }

  ###########################################################################
  # 3) EnRngeXd block
  #
  #    % EnRngeXd
  #      0.000000 | 10.000000 | eV
  ###########################################################################
  if (($1 == "%" && $2 == "EnRngeXd") || $1 == "%EnRngeXd") {
    print

    if (getline > 0) {
      printf "  %s | %s |         eV    # [Xd] Energy range\n", \
             emin, emax
    }

    next
  }

  ###########################################################################
  # 4) ETStpsXd line
  #
  #    ETStpsXd= 200
  ###########################################################################
  if ($1 == "ETStpsXd=" || $1 ~ /^ETStpsXd=/) {
    printf "ETStpsXd= %s                    # [Xd] Total Energy steps\n", \
           esteps
    next
  }

  ###########################################################################
  # 5) LongDrXd block
  #
  #    % LongDrXd
  #      1.000000 | 0.000000 | 0.000000 |
  ###########################################################################
  if (($1 == "%" && $2 == "LongDrXd") || $1 == "%LongDrXd") {
    print

    if (getline > 0) {
      printf " %s | %s | %s |        # [Xd] [cc] Electric Field\n", \
             ex, ey, ez
    }

    next
  }

  ###########################################################################
  # Default: preserve everything exactly as generated
  ###########################################################################
  print
}
' "$infile" > "$outfile"

###############################################################################
# Summary
###############################################################################

echo "Patched input written to $outfile"
echo "  q-points      = $iq $fq"
echo "  nl_band_i     = $nl_band_i"
echo "  nl_band_f     = $nl_band_f"
echo "  nl_en_min     = $nl_en_min"
echo "  nl_en_max     = $nl_en_max"
echo "  nl_esteps     = $nl_esteps"
echo "  efield1_dir   = $efield1_x $efield1_y $efield1_z"