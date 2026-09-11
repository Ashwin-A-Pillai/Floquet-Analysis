#!/usr/bin/env bash
#
# patch_bse.sh — tweak a Yambo BSE input (readable style, preserves structure)
#   • RandQpts        -> ${Rand_qpts} (default 1000000)
#   • RandGvec        -> ${Rand_gvec} (default 100)
#   • CUT block       -> CUTGeo="slab z" (collapse block)
#   • BSEmod          -> "retarded"
#   • BSENGexx (Ry)   -> ${BSE_Ex} (default 40)
#   • BSENGBlk (Ry)   -> ${BSE_Blk} (default 8)
#   • Lkind           -> "full"
#   • #ALLGexx        -> ALLGexx (uncomment)
#   • % BSEQptR       -> 1 | ${f_qk} |
#   • KfnQPdb → "E < GW0/ndb.QP"
#   • % BEnRange → $iE | $fE
#   • BEnSteps → $stepE
#   • % BndsRnXp → $Polar_band
#   • NGsBlkXp → ${Block_size} (Ry) if present
#   • #WRbsWF → WRbsWF (uncomment)
#   • % QpntsRXp → 1 | ${Xqk} |
#   • KEEP: em1d and PPA; REMOVE: em1s
#   • DELETE: BSKIOmode="2D_standard"
#   • REMOVE Xs (static screening) blocks/lines:
#       % QpntsRXs, % BndsRnXs, % DmRngeXs, % EhEngyXs, % LongDrXs
#       NGsBlkXs, GrFnTpXs, CGrdSpXs, DrudeWXs
#
# Usage: ./patch_bse.sh <input_file> [output_file]
#
# Defaults (override with exports before running):
  export Rand_qpts=${Rand_qpts:-1000000}
  export Rand_gvec=${Rand_gvec:-100}
  export BSE_Ex=${BSE_Ex:-40}
  export BSE_Blk=${BSE_Blk:-8}
  export f_qk=${f_qk:-195}
  export Block_size=${Block_size:-3000}
  export iE=${iE:-5.0}
  export fE=${fE:-20.0}
  export i_corband=${i_corband:-2}
  export f_corband=${f_corband:-7}
  export stepE=${stepE:-400}
  export Polar_band=${Polar_band:-80}
  export Xqk=${Xqk:-1}

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <input_file> [output_file]" >&2
  exit 1
fi

infile="$1"
outfile="${2:-${infile}.patched}"

# Require the variables to be set (defaults were exported above, but allow user override)
: "${Rand_qpts:?}" \
  "${Rand_gvec:?}" \
  "${BSE_Ex:?}" \
  "${BSE_Blk:?}" \
  "${f_qk:?}" \
  "${Block_size:?}" \
  "${iE:?Please set \$iE}" \
  "${fE:?Please set \$fE}" \
  "${i_corband:?Please set \$i_corband}" \
  "${f_corband:?Please set \$f_corband}" \
  "${stepE:?Please set \$stepE}" \
  "${Polar_band:?Please set \$Polar_band}" \
  "${Xqk:?}"

awk -v ie="$iE" \
    -v fe="$fE" \
    -v ic="$i_corband" \
    -v fc="$f_corband" \
    -v se="$stepE" \
    -v pb="$Polar_band" \
    -v rq="$Rand_qpts" \
    -v rg="$Rand_gvec" \
    -v bex="$BSE_Ex" \
    -v bblk="$BSE_Blk" \
    -v fqk="$f_qk" \
    -v nbx="$Block_size" \
    -v xqk="$Xqk" '
# ---------- Block A: optics-style tweaks (preserved) ----------
{
  # 1) KfnQPdb tweak
  if ($1 ~ /^KfnQPdb=/) {
    print "KfnQPdb= \"E < GW0/ndb.QP\"         # [EXTQP BSK BSS] Database action"
    next
  }

  # 2) % BEnRange block
  if ($1=="%" && $2=="BEnRange") {
    print                           # “% BEnRange” line
    getline                         # old range line
    printf "  %s | %s |         eV    # [BSS] Energy range\n", ie, fe
    next
  }

  # 3) BEnSteps
  if ($1=="BEnSteps=") {
    printf "BEnSteps= %s                    # [BSS] Energy steps\n", se
    next
  }

  # 4) % BndsRnXp block (we keep Xp; may still adjust other params elsewhere)
  if ($1=="%" && $2=="BndsRnXp") {
    print                           # “% BndsRnXp” line
    getline                         # old bands line
    printf "   1 |  %s |                         # [Xp] Polarization function bands\n", pb
    next
  }

  # 5) % BSEBands block
  if ($1=="%" && $2=="BSEBands") {
    print                           # “% BSEBands” line
    getline                         # old bands line
    printf "  %s  |  %s |", ic, fc
    next
  }

  # fall-through: do not print here; let Block B handle default print
}

# ---------- Block B: new BSE/Yambo tweaks + removals ----------
BEGIN { in_cut = 0; in_xs = 0; }
{
  # KEEP em1d & PPA; REMOVE em1s (case-insensitive, token at column 1)
  tl1 = tolower($1)
  if (tl1 == "em1s") next

  # DELETE: BSKIOmode line entirely
  if ($1 ~ /^BSKIOmode=/) next

  # RandQpts
  if ($1 ~ /^RandQpts=/) {
    print "RandQpts=" rq "                       # [RIM] Number of random q-points in the BZ"
    next
  }

  # RandGvec
  if ($1 ~ /^RandGvec=/) {
    print "RandGvec= " rg "            RL    # [RIM] Coulomb interaction RS components"
    next
  }

  # CUT block: collapse to single CUTGeo line
  if ($1 ~ /^CUTGeo=/) {
    print "CUTGeo= \"slab z\"               # [CUT] Coulomb Cutoff geometry: box/cylinder/sphere X/Y/Z/XY.."
    in_cut = 1
    next
  }
  if (in_cut) {
    if ($1 ~ /^CUTwsGvec=/) { in_cut = 0; next }
    if ($1 ~ /^%$/ || ($1=="%" && $2=="CUTBox") || $1 ~ /^CUTRadius=/ || $1 ~ /^CUTCylLen=/) { next }
    in_cut = 0
  }

  # BSEmod
  if ($1 ~ /^BSEmod=/) {
    print "BSEmod= \"retarded\"               # [BSE] resonant/retarded/coupling"
    next
  }

  # BSENGexx (Ry) and BSENGBlk (Ry)
  if ($1 ~ /^BSENGexx=/) {
    print "BSENGexx=  " bex "            Ry    # [BSK] Exchange components"
    next
  }
  if ($1 ~ /^BSENGBlk=/) {
    print "BSENGBlk=" bblk "         Ry    # [BSK] Screened interaction block size [if -1 uses all the G-vectors of W(q,G,Gp)]"
    next
  }

  # Lkind
  if ($1 ~ /^Lkind=/) {
    print "Lkind=\"full\"                  #[BSE,X] bar(default)/full/tilde"
    next
  }

  # Uncomment ALLGexx (handles lines with trailing comments too)
  if ($1 ~ /^#?ALLGexx$/) {
    print "ALLGexx                      # [BSS] Force the use use all RL vectors for the exchange part"
    next
  }

  # Uncomment WRbsWF (handles lines with trailing comments too)
  if ($1 ~ /^#?WRbsWF$/) {
    print "WRbsWF                        # [BSS] Write to disk excitonic the WFs"
    next
  }

  # % BSEQptR block
  if ($1=="%" && $2=="BSEQptR") {
    print
    getline
    printf " 1 | %s |                     # [BSK] Transferred momenta range\n", fqk
    getline
    print
    next
  }

  # % QpntsRXp block (we KEEP Xp; set it to 1 | Xqk |)
  if ($1=="%" && $2=="QpntsRXp") {
    print
    getline
    printf "   1 | %s |                         # [Xp] Transferred momenta\n", xqk
    getline
    print
    next
  }

  # NGsBlkXp (Ry), if present
  if ($1 ~ /^NGsBlkXp=/) {
    print "NGsBlkXp= " nbx "                Ry    # [Xp] Polarization function block size"
    next
  }

  # ----------------- REMOVALS: Xs (static screening) -----------------
  # Remove entire % QpntsRXs, % BndsRnXs, % DmRngeXs, % EhEngyXs, % LongDrXs blocks
  if ($1=="%" && ($2=="QpntsRXs" || $2=="BndsRnXs" || $2=="DmRngeXs" || $2=="EhEngyXs" || $2=="LongDrXs")) {
    in_xs = 1; next
  }
  if (in_xs) {
    if ($1=="%") { in_xs = 0; next } else { next }
  }

  # Remove single-line Xs settings
  if ($1 ~ /^NGsBlkXs=/) next
  if ($1 ~ /^GrFnTpXs=/) next
  if ($1 ~ /^CGrdSpXs=/) next
  if ($1 ~ /^FFTGvecs=/) next
  if ($1 ~ /^DrudeWXs=/) next

  # default passthrough
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