#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Defaults (override by exporting env vars before running)
###############################################################################
: "${k_final:=29}"
: "${iE:=5.0}"
: "${fE:=20.0}"
: "${i_corband:=2}"
: "${f_corband:=7}"
: "${stepE:=400}"

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <input_file> [output_file]" >&2
  exit 1
fi

infile="$1"
outfile="${2:-${infile}.modified}"

awk -v kf="$k_final" \
    -v ie="$iE" \
    -v fe="$fE" \
    -v ic="$i_corband" \
    -v fc="$f_corband" \
    -v se="$stepE" '
function comment_part(line,   c) {
  c = ""
  if (match(line, /#.*/)) c = substr(line, RSTART)
  return c
}

# --- % QpntsRXd block: replace the next (numeric) line ---
/^[[:space:]]*%[[:space:]]*QpntsRXd([[:space:]]|$)/ {
  print
  if (getline line) {
    c = comment_part(line)
    printf "   1 |  %s |", kf
    if (c != "") printf "                         %s", c
    printf "\n"
  }
  next
}

# --- % BndsRnXd block: replace the next (numeric) line ---
/^[[:space:]]*%[[:space:]]*BndsRnXd([[:space:]]|$)/ {
  print
  if (getline line) {
    c = comment_part(line)
    printf "  %s | %s |", ic, fc
    if (c != "") printf "                           %s", c
    printf "\n"
  }
  next
}

# --- % EnRngeXd block: replace the next (numeric) line ---
/^[[:space:]]*%[[:space:]]*EnRngeXd([[:space:]]|$)/ {
  print
  if (getline line) {
    c = comment_part(line)
    # Keep unit as eV (as in your example)
    printf "  %s | %s |         eV", ie, fe
    if (c != "") printf "    %s", c
    printf "\n"
  }
  next
}

# --- ETStpsXd= line: replace the value ---
/^[[:space:]]*ETStpsXd[[:space:]]*=/ {
  c = comment_part($0)
  printf "ETStpsXd= %s", se
  if (c != "") printf "                    %s", c
  printf "\n"
  next
}

# default passthrough
{ print }
' "$infile" > "$outfile"

echo "Patched input written to: $outfile"
echo "  k_final   = $k_final"
echo "  i_corband = $i_corband"
echo "  f_corband = $f_corband"
echo "  iE        = $iE"
echo "  fE        = $fE"
echo "  stepE     = $stepE"
