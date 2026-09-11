#!/usr/bin/env bash
#
###############################################################################
# Defaults with override behavior
###############################################################################
: "${rand_qpts:=5000000}"
: "${cut_geo:=slab z}"
: "${ngs_blk_xs:=3000}"
: "${harr_rlvcs:=3000}"
: "${exx_rlvcs:=3000}"
: "${corr_rlvcs:=3000}"

export rand_qpts cut_geo ngs_blk_xs \
       harr_rlvcs exx_rlvcs corr_rlvcs

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <input_file> [output_file]" >&2
  exit 1
fi

infile="$1"
outfile="${2:-${infile}.modified}"

awk -v rq="$rand_qpts" \
    -v cg="$cut_geo" \
    -v nb="$ngs_blk_xs" \
    -v ha="$harr_rlvcs" \
    -v xx="$exx_rlvcs" \
    -v co="$corr_rlvcs" '
{
  ###########################################################################
  # 1) RandQpts line:
  #    RandQpts=0
  ###########################################################################
  if ($1 ~ /^RandQpts=/) {
    printf "RandQpts=%s                       # [RIM] Number of random q-points in the BZ\n", rq
    next
  }

  ###########################################################################
  # 2) CUTGeo line:
  #    CUTGeo= "none"
  ###########################################################################
  if ($1 == "CUTGeo=") {
    printf "CUTGeo= \"%s\"                   # [CUT] Coulomb Cutoff geometry: box/cylinder/sphere/ws/slab X/Y/Z/XY..\n", cg
    next
  }

  ###########################################################################
  # 3) NGsBlkXs line:
  #    NGsBlkXs= 1                RL
  ###########################################################################
  if ($1 == "NGsBlkXs=") {
    printf "NGsBlkXs= %s mHa                # [Xs] Response block size\n", nb
    next
  }

  ###########################################################################
  # 4) HARRLvcs line:
  #    HARRLvcs=  6735            RL
  ###########################################################################
  if ($1 == "HARRLvcs=") {
    printf "HARRLvcs=  %s            mHa    # [HA] Hartree     RL components\n", ha
    next
  }

  ###########################################################################
  # 5) EXXRLvcs line:
  #    EXXRLvcs=  6735            RL
  ###########################################################################
  if ($1 == "EXXRLvcs=") {
    printf "EXXRLvcs=  %s            mHa    # [XX] Exchange    RL components\n", xx
    next
  }

  ###########################################################################
  # 6) CORRLvcs line:
  #    CORRLvcs=  6735            RL
  ###########################################################################
  if ($1 == "CORRLvcs=") {
    printf "CORRLvcs=  %s            mHa    # [GW] Correlation RL components\n", co
    next
  }

  # default: pass through unchanged
  print
}
' "$infile" > "$outfile"

echo "Patched input written to $outfile"
echo "  rand_qpts  = $rand_qpts"
echo "  cut_geo    = $cut_geo"
echo "  ngs_blk_xs = $ngs_blk_xs"
echo "  harr_rlvcs = $harr_rlvcs"
echo "  exx_rlvcs  = $exx_rlvcs"
echo "  corr_rlvcs = $corr_rlvcs"