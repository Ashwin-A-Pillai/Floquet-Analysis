#!/usr/bin/env bash
#
###############################################################################
# Defaults with override behavior
###############################################################################

: "${RQpts:=5000000}"
: "${RGvec:=29}"
: "${cgeo:=slab z}"

: "${exx_rlvcs:=3000}"
: "${vxc_rlvcs:=3000}"
: "${ngs_blk_xp:=3000}"

: "${xp_band_i:=1}"
: "${xp_band_f:=15}"

: "${efield1_x:=1.000000}"
: "${efield1_y:=0.000000}"
: "${efield1_z:=0.000000}"

: "${Gbnd:=40}"

: "${iq:=1}"
: "${fq:=1}"
: "${i_band:=1}"
: "${f_band:=1}"

export RQpts RGvec cgeo \
       exx_rlvcs vxc_rlvcs ngs_blk_xp \
       xp_band_i xp_band_f \
       efield1_x efield1_y efield1_z \
       Gbnd iq fq i_band f_band

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
    -v exx="$exx_rlvcs" \
    -v vxc="$vxc_rlvcs" \
    -v ngs="$ngs_blk_xp" \
    -v xbi="$xp_band_i" \
    -v xbf="$xp_band_f" \
    -v ex="$efield1_x" \
    -v ey="$efield1_y" \
    -v ez="$efield1_z" \
    -v gb="$Gbnd" \
    -v iqv="$iq" \
    -v fqv="$fq" \
    -v ib="$i_band" \
    -v fb="$f_band" '
{
  ###########################################################################
  # 1) RandQpts
  #
  #    RandQpts= 5000000
  ###########################################################################
  if ($1 == "RandQpts=" || $1 ~ /^RandQpts=/) {
    printf "RandQpts= %s                       # [RIM] Number of random q-points in the BZ\n", rq
    next
  }

  ###########################################################################
  # 2) RandGvec
  #
  #    RandGvec= 29 RL
  ###########################################################################
  if ($1 == "RandGvec=" || $1 ~ /^RandGvec=/) {
    printf "RandGvec= %s                RL    # [RIM] Coulomb interaction RS components\n", rg
    next
  }

  ###########################################################################
  # 3) CUTGeo
  #
  #    CUTGeo= "slab z"
  ###########################################################################
  if ($1 == "CUTGeo=" || $1 ~ /^CUTGeo=/) {
    printf "CUTGeo= \"%s\"                   # [CUT] Coulomb Cutoff geometry: box/cylinder/sphere/ws/slab X/Y/Z/XY..\n", cg
    next
  }

  ###########################################################################
  # 4) EXXRLvcs
  ###########################################################################
  if ($1 == "EXXRLvcs=" || $1 ~ /^EXXRLvcs=/) {
    printf "EXXRLvcs= %s                mHa    # [XX] Exchange    RL components\n", exx
    next
  }

  ###########################################################################
  # 5) VXCRLvcs
  ###########################################################################
  if ($1 == "VXCRLvcs=" || $1 ~ /^VXCRLvcs=/) {
    printf "VXCRLvcs= %s                mHa    # [XC] XCpotential RL components\n", vxc
    next
  }

  ###########################################################################
  # 6) BndsRnXp block
  #
  #    % BndsRnXp
  #      1 | 15 |
  ###########################################################################
  if (($1 == "%" && $2 == "BndsRnXp") || $1 == "%BndsRnXp") {
    print

    if (getline > 0) {
      printf "  %s | %s |                           # [Xp] Polarization function bands\n", \
             xbi, xbf
    }

    next
  }

  ###########################################################################
  # 7) NGsBlkXp
  #
  #    NGsBlkXp= 3000 mHa
  ###########################################################################
  if ($1 == "NGsBlkXp=" || $1 ~ /^NGsBlkXp=/) {
    printf "NGsBlkXp= %s                mHa    # [Xp] Response block size\n", ngs
    next
  }

  ###########################################################################
  # 8) LongDrXp block
  #
  #    % LongDrXp
  #      1.000000 | 0.000000 | 1.000000 |
  ###########################################################################
  if (($1 == "%" && $2 == "LongDrXp") || $1 == "%LongDrXp") {
    print

    if (getline > 0) {
      printf " %s | %s | %s |        # [Xp] [cc] Electric Field\n", \
             ex, ey, ez
    }

    next
  }

  ###########################################################################
  # 9) GbndRnge block
  #
  #    % GbndRnge
  #      1 | 40 |
  ###########################################################################
  if (($1 == "%" && $2 == "GbndRnge") || $1 == "%GbndRnge") {
    print

    if (getline > 0) {
      printf "  1 | %s |                           # [GW] G[W] bands range\n", gb
    }

    next
  }

  ###########################################################################
  # 10) QPkrange block
  #
  #     % QPkrange
  #       iq | fq | i_band | f_band |
  #
  # Also handles:
  #
  #     %QPkrange
  ###########################################################################
  if (($1 == "%" && $2 == "QPkrange") || $1 == "%QPkrange") {
    print

    if (getline > 0) {
      printf " %s | %s | %s | %s |\n", \
             iqv, fqv, ib, fb
    }

    next
  }

  ###########################################################################
  # Default: preserve everything else unchanged
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
echo "  exx_rlvcs    = $exx_rlvcs"
echo "  vxc_rlvcs    = $vxc_rlvcs"
echo "  ngs_blk_xp   = $ngs_blk_xp"
echo "  xp_bands     = $xp_band_i $xp_band_f"
echo "  efield1_dir  = $efield1_x $efield1_y $efield1_z"
echo "  Gbnd         = $Gbnd"
echo "  QPkrange     = $iq $fq $i_band $f_band"