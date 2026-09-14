#!/usr/bin/env bash
#
###############################################################################
# Defaults with override behavior
###############################################################################

# Random integration / Coulomb cutoff
: "${RQpts:=5000000}"
: "${RGvec:=29}"
: "${cgeo:=slab z}"

# Static screening
: "${xp_band_i:=1}"
: "${xp_band_f:=15}"
: "${ngs_blk_xp:=3000}"

# Electric-field direction
: "${efield1_x:=1.000000}"
: "${efield1_y:=0.000000}"
: "${efield1_z:=0.000000}"

# Collision bands
: "${i_cb:=4}"
: "${f_cb:=5}"

# Hartree / exchange / correlation cutoffs
: "${har_rlvcs:=3000}"
: "${exx_rlvcs:=3000}"
: "${cor_rlvcs:=3000}"

export RQpts RGvec cgeo \
       xp_band_i xp_band_f ngs_blk_xp \
       efield1_x efield1_y efield1_z \
       i_cb f_cb \
       har_rlvcs exx_rlvcs cor_rlvcs

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
    -v xbi="$xp_band_i" \
    -v xbf="$xp_band_f" \
    -v ngs="$ngs_blk_xp" \
    -v ex="$efield1_x" \
    -v ey="$efield1_y" \
    -v ez="$efield1_z" \
    -v icb="$i_cb" \
    -v fcb="$f_cb" \
    -v har="$har_rlvcs" \
    -v exx="$exx_rlvcs" \
    -v corr="$cor_rlvcs" '
{
  ###########################################################################
  # 1) RandQpts
  #
  #    RandQpts= 5000000
  ###########################################################################
  if ($1 == "RandQpts=" || $1 ~ /^RandQpts=/) {
    printf "RandQpts= %s                         # [RIM] Number of random q-points in the BZ\n", rq
    next
  }

  ###########################################################################
  # 2) RandGvec
  #
  #    RandGvec= 29 RL
  ###########################################################################
  if ($1 == "RandGvec=" || $1 ~ /^RandGvec=/) {
    printf "RandGvec= %s                 RL    # [RIM] Coulomb interaction RS components\n", rg
    next
  }

  ###########################################################################
  # 3) CUTGeo
  #
  #    CUTGeo= "slab z"
  ###########################################################################
  if ($1 == "CUTGeo=" || $1 ~ /^CUTGeo=/) {
    printf "CUTGeo= \"%s\"                    # [CUT] Coulomb Cutoff geometry: box/cylinder/sphere/ws/slab X/Y/Z/XY..\n", cg
    next
  }

  ###########################################################################
  # 4) BndsRnXs block
  #
  #    % BndsRnXs
  #      1 | 15 |
  ###########################################################################
  if (($1 == "%" && $2 == "BndsRnXs") || $1 == "%BndsRnXs") {
    print

    if (getline > 0) {
      printf "  %s | %s |                           # [Xs] Polarization function bands\n", \
             xbi, xbf
    }

    next
  }

  ###########################################################################
  # 5) NGsBlkXs
  #
  #    NGsBlkXs= 3000 mHa
  ###########################################################################
  if ($1 == "NGsBlkXs=" || $1 ~ /^NGsBlkXs=/) {
    printf "NGsBlkXs= %s          mHa    # [Xs] Response block size\n", ngs
    next
  }

  ###########################################################################
  # 6) LongDrXs block
  #
  #    % LongDrXs
  #      1.000000 | 0.000000 | 1.000000 |
  ###########################################################################
  if (($1 == "%" && $2 == "LongDrXs") || $1 == "%LongDrXs") {
    print

    if (getline > 0) {
      printf " %s | %s | %s |         # [Xs] [cc] Electric Field\n", \
             ex, ey, ez
    }

    next
  }

  ###########################################################################
  # 7) COLLBands block
  #
  #    % COLLBands
  #      4 | 5 |
  ###########################################################################
  if (($1 == "%" && $2 == "COLLBands") || $1 == "%COLLBands") {
    print

    if (getline > 0) {
      printf "  %s | %s |                           # [COLL] Bands for the collisions\n", \
             icb, fcb
    }

    next
  }

  ###########################################################################
  # 8) HARRLvcs
  #
  #    HARRLvcs= 3000 mHa
  ###########################################################################
  if ($1 == "HARRLvcs=" || $1 ~ /^HARRLvcs=/) {
    printf "HARRLvcs= %s            mHa    # [HA] Hartree     RL components\n", har
    next
  }

  ###########################################################################
  # 9) EXXRLvcs
  #
  #    EXXRLvcs= 3000 mHa
  ###########################################################################
  if ($1 == "EXXRLvcs=" || $1 ~ /^EXXRLvcs=/) {
    printf "EXXRLvcs= %s            mHa    # [XX] Exchange    RL components\n", exx
    next
  }

  ###########################################################################
  # 10) CORRLvcs
  #
  #     CORRLvcs= 3000 mHa
  ###########################################################################
  if ($1 == "CORRLvcs=" || $1 ~ /^CORRLvcs=/) {
    printf "CORRLvcs= %s            mHa    # [GW] Correlation RL components [if -1 uses all the G-vectors of W_{G,Gp}]\n", corr
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
echo "  XsBands      = $xp_band_i $xp_band_f"
echo "  NGsBlkXs     = $ngs_blk_xp mHa"
echo "  field_dir    = $efield1_x $efield1_y $efield1_z"
echo "  COLLBands    = $i_cb $f_cb"
echo "  HARRLvcs     = $har_rlvcs mHa"
echo "  EXXRLvcs     = $exx_rlvcs mHa"
echo "  CORRLvcs     = $cor_rlvcs mHa"