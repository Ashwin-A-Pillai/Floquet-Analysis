#!/usr/bin/env bash
#
###############################################################################
# Defaults with override behavior
###############################################################################

# Random integration / Coulomb cutoff
: "${RQpts:=5000000}"
: "${RGvec:=29}"
: "${cgeo:=slab z}"

# BSE kernel
: "${BSE_Ex:=20000}"

# QP database
: "${QPdir:=E < GW0/ndb.QP}"

# BSE q-point range
: "${iq:=1}"
: "${fq:=1}"

# BSE bands
: "${i_band:=1}"
: "${f_band:=15}"

# BSE spectrum
: "${iE:=0.000000}"
: "${fE:=10.000000}"
: "${stepE:=100}"

# Polarization bands for W
: "${xp_band_i:=1}"
: "${xp_band_f:=15}"
: "${ngs_blk_xp:=3000}"

# Electric-field direction
: "${efield1_x:=1.000000}"
: "${efield1_y:=0.000000}"
: "${efield1_z:=0.000000}"

export RQpts RGvec cgeo \
       BSE_Ex QPdir \
       iq fq i_band f_band \
       iE fE stepE \
       xp_band_i xp_band_f ngs_blk_xp \
       efield1_x efield1_y efield1_z

set -euo pipefail

###############################################################################
# Arguments
###############################################################################

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
# Everything else in the original Lumen/Yambo input is preserved unchanged.
###############################################################################

awk -v rq="$RQpts" \
    -v rg="$RGvec" \
    -v cg="$cgeo" \
    -v bex="$BSE_Ex" \
    -v qpdir="$QPdir" \
    -v iqv="$iq" \
    -v fqv="$fq" \
    -v ib="$i_band" \
    -v fb="$f_band" \
    -v emin="$iE" \
    -v emax="$fE" \
    -v esteps="$stepE" \
    -v xpbi="$xp_band_i" \
    -v xpbf="$xp_band_f" \
    -v ngs="$ngs_blk_xp" \
    -v ex="$efield1_x" \
    -v ey="$efield1_y" \
    -v ez="$efield1_z" '
{
  ###########################################################################
  # 1) RandQpts
  ###########################################################################
  if ($1 == "RandQpts=" || $1 ~ /^RandQpts=/) {
    printf "RandQpts= %s                        # [RIM] Number of random q-points in the BZ\n", rq
    next
  }

  ###########################################################################
  # 2) RandGvec
  ###########################################################################
  if ($1 == "RandGvec=" || $1 ~ /^RandGvec=/) {
    printf "RandGvec= %s                RL    # [RIM] Coulomb interaction RS components\n", rg
    next
  }

  ###########################################################################
  # 3) CUTGeo
  ###########################################################################
  if ($1 == "CUTGeo=" || $1 ~ /^CUTGeo=/) {
    printf "CUTGeo= \"%s\"                    # [CUT] Coulomb Cutoff geometry: box/cylinder/sphere/ws/slab X/Y/Z/XY..\n", cg
    next
  }

  ###########################################################################
  # 4) BSENGexx
  ###########################################################################
  if ($1 == "BSENGexx=" || $1 ~ /^BSENGexx=/) {
    printf "BSENGexx= %s              mHa    # [BSK] Exchange components\n", bex
    next
  }

  ###########################################################################
  # 5) BSENGBlk
  #
  #    Always use all G-vectors available in the screening database.
  ###########################################################################
  if ($1 == "BSENGBlk=" || $1 ~ /^BSENGBlk=/) {
    printf "BSENGBlk= -1              RL    # [BSK] Screened interaction block size [if -1 uses all the G-vectors of W_{G,Gp}]\n"
    next
  }

  ###########################################################################
  # 6) KfnQPdb
  ###########################################################################
  if ($1 == "KfnQPdb=" || $1 ~ /^KfnQPdb=/) {
    printf "KfnQPdb= \"%s\"                  # [EXTQP BSK BSS] Database action\n", qpdir
    next
  }

  ###########################################################################
  # 7) BSEQptR
  ###########################################################################
  if (($1 == "%" && $2 == "BSEQptR") || $1 == "%BSEQptR") {
    print

    if (getline > 0) {
      printf " %s | %s |                             # [BSK] Transferred momenta range\n", \
             iqv, fqv
    }

    next
  }

  ###########################################################################
  # 8) BSEBands
  ###########################################################################
  if (($1 == "%" && $2 == "BSEBands") || $1 == "%BSEBands") {
    print

    if (getline > 0) {
      printf "  %s | %s |                           # [BSK] Bands range\n", \
             ib, fb
    }

    next
  }

  ###########################################################################
  # 9) BEnRange
  ###########################################################################
  if (($1 == "%" && $2 == "BEnRange") || $1 == "%BEnRange") {
    print

    if (getline > 0) {
      printf "  %s | %s |         eV    # [BSS] Energy range\n", \
             emin, emax
    }

    next
  }

  ###########################################################################
  # 10) BEnSteps
  ###########################################################################
  if ($1 == "BEnSteps=" || $1 ~ /^BEnSteps=/) {
    printf "BEnSteps= %s                    # [BSS] Energy steps\n", esteps
    next
  }

  ###########################################################################
  # 11) BLongDir
  ###########################################################################
  if (($1 == "%" && $2 == "BLongDir") || $1 == "%BLongDir") {
    print

    if (getline > 0) {
      printf " %s | %s | %s |       # [BSS] [cc] Electric Field versor\n", \
             ex, ey, ez
    }

    next
  }

  ###########################################################################
  # 12) BndsRnXp
  ###########################################################################
  if (($1 == "%" && $2 == "BndsRnXp") || $1 == "%BndsRnXp") {
    print

    if (getline > 0) {
      printf "  %s | %s |                           # [Xp] Polarization function bands\n", \
             xpbi, xpbf
    }

    next
  }

  ###########################################################################
  # 13) NGsBlkXp
  ###########################################################################
  if ($1 == "NGsBlkXp=" || $1 ~ /^NGsBlkXp=/) {
    printf "NGsBlkXp= %s        mHa   # [Xp] Response block size\n", ngs
    next
  }

  ###########################################################################
  # 14) LongDrXp
  ###########################################################################
  if (($1 == "%" && $2 == "LongDrXp") || $1 == "%LongDrXp") {
    print

    if (getline > 0) {
      printf " 0.000000 | 0.100000E-4 | 0.000000 |       # [Xp] [cc] Electric Field\n"
    }

    next
  }

  ###########################################################################
  # Default: preserve everything else exactly as generated
  ###########################################################################
  print
}
' "$infile" > "$outfile"

###############################################################################
# Summary
###############################################################################

echo "Patched input written to $outfile"
echo "  RQpts        = $RQpts"
echo "  RGvec        = $RGvec"
echo "  cgeo         = $cgeo"
echo "  BSE_Ex       = $BSE_Ex mHa"
echo "  BSENGBlk     = -1 RL"
echo "  QPdir        = $QPdir"
echo "  BSEQptR      = $iq $fq"
echo "  BSEBands     = $i_band $f_band"
echo "  BEnRange     = $iE $fE eV"
echo "  BEnSteps     = $stepE"
echo "  XpBands      = $xp_band_i $xp_band_f"
echo "  NGsBlkXp     = $ngs_blk_xp mHa"
echo "  field_dir    = $efield1_x $efield1_y $efield1_z"