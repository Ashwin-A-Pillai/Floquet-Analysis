#!/usr/bin/env bash
#
###############################################################################
# Defaults with override behavior
###############################################################################
: "${xorder:=4}"
: "${time_min:=-40.000000}"
: "${time_max:=-1.000000}"

: "${energy_steps:=200}"
: "${energy_min:=0.000000}"
: "${energy_max:=10.000000}"

: "${damp_mode:=NONE}"
: "${damp_factor:=0.10000}"

: "${pump_path:=none}"

export xorder time_min time_max \
       energy_steps energy_min energy_max \
       damp_mode damp_factor pump_path

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
# Write through a temporary file so that input and output can safely be the
# same file.
###############################################################################
outdir="$(dirname -- "$outfile")"
outbase="$(basename -- "$outfile")"

mkdir -p -- "$outdir"

tmpfile="$(mktemp "${outdir}/.${outbase}.tmp.XXXXXX")"

cleanup() {
  rm -f -- "$tmpfile"
}

trap cleanup EXIT

awk -v xo="$xorder" \
    -v tmin="$time_min" \
    -v tmax="$time_max" \
    -v esteps="$energy_steps" \
    -v emin="$energy_min" \
    -v emax="$energy_max" \
    -v dmode="$damp_mode" \
    -v dfactor="$damp_factor" \
    -v ppath="$pump_path" '
{
  ###########################################################################
  # 1) Xorder line:
  #
  #    Xorder= 1
  ###########################################################################
  if ($1 == "Xorder=") {
    printf "Xorder= %s                        # Max order of the response/exc functions\n", xo
    next
  }

  ###########################################################################
  # 2) TimeRange block:
  #
  #    % TimeRange
  #    -1.000000 |-1.000000 | fs
  ###########################################################################
  if (($1 == "%" && $2 == "TimeRange") || $1 == "%TimeRange") {
    print

    if (getline > 0) {
      printf "  %s | %s |         fs    # Time-window where processing is done\n", \
             tmin, tmax
    }

    next
  }

  ###########################################################################
  # 3) ETStpsRt line:
  #
  #    ETStpsRt= 200
  ###########################################################################
  if ($1 == "ETStpsRt=") {
    printf "ETStpsRt= %s                    # Total Energy steps\n", esteps
    next
  }

  ###########################################################################
  # 4) EnRngeRt block:
  #
  #    % EnRngeRt
  #      0.000000 | 10.000000 | eV
  ###########################################################################
  if (($1 == "%" && $2 == "EnRngeRt") || $1 == "%EnRngeRt") {
    print

    if (getline > 0) {
      printf "  %s | %s |         eV    # Energy range\n", emin, emax
    }

    next
  }

  ###########################################################################
  # 5) DampMode line:
  #
  #    DampMode= "NONE"
  ###########################################################################
  if ($1 == "DampMode=") {
    printf "DampMode= \"%s\"                 # Damping type ( NONE | LORENTZIAN | GAUSSIAN )\n", \
           dmode
    next
  }

  ###########################################################################
  # 6) DampFactor line:
  #
  #    DampFactor= 0.000000 eV
  ###########################################################################
  if ($1 == "DampFactor=") {
    printf "DampFactor= %s       eV    # Damping parameter\n", dfactor
    next
  }

  ###########################################################################
  # 7) PumpPATH line:
  #
  #    PumpPATH= "none"
  ###########################################################################
  if ($1 == "PumpPATH=") {
    printf "PumpPATH= \"%s\"                 # Path of the simulation with the Pump only\n", \
           ppath
    next
  }

  ###########################################################################
  # Default: pass through unchanged
  ###########################################################################
  print
}
' "$infile" > "$tmpfile"

mv -f -- "$tmpfile" "$outfile"
trap - EXIT

###############################################################################
# Print the applied settings
###############################################################################
echo "Patched input written to $outfile"
echo "  xorder       = $xorder"
echo "  time_min     = $time_min"
echo "  time_max     = $time_max"
echo "  energy_steps = $energy_steps"
echo "  energy_min   = $energy_min"
echo "  energy_max   = $energy_max"
echo "  damp_mode    = $damp_mode"
echo "  damp_factor  = $damp_factor"
echo "  pump_path    = $pump_path"