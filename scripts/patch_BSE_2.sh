#!/usr/bin/env bash
# convert_bse_input.sh — convert an initial Yambo BSE input to your target spec
# Usage: ./convert_bse_input.sh <input_file> [output_file]

###############################################################################
# User-tunable defaults (override by exporting before run), e.g.:
#   export BSE_Ex=14000
#   export Polar_band=100
###############################################################################
: "${BSE_Ex:=8000}"      # => bex (mHa)
: "${BSE_Blk:=6000}"      # => bblk (mHa)
: "${f_qk:=1}"
: "${Block_size:=6000}"
: "${iE:=5.0}"
: "${fE:=20.0}"
: "${i_corband:=53}"
: "${f_corband:=58}"
: "${stepE:=400}"
: "${Polar_band:=120}"
: "${Xqk:=1}"
: "${lkd:=bar}"
: "${QPdir:=E < GW0/ndb.QP}"

export BSE_Ex BSE_Blk f_qk Block_size iE fE i_corband f_corband \
       stepE Polar_band Xqk lkd QPdir

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <input_file> [output_file]" >&2
  exit 1
fi

infile="$1"
outfile="${2:-${infile%.in}_converted.in}"

# sanity
: "${BSE_Ex:?}" "${BSE_Blk:?}" "${f_qk:?}" \
  "${Block_size:?}" "${iE:?Please set \$iE}" "${fE:?Please set \$fE}" \
  "${i_corband:?Please set \$i_corband}" "${f_corband:?Please set \$f_corband}" \
  "${stepE:?Please set \$stepE}" "${Polar_band:?Please set \$Polar_band}" \
  "${Xqk:?}" "${lkd:?Please set \$lkd}" "${QPdir:?Please set \$QPdir}"

tmp1="$(mktemp)"; trap 'rm -f "$tmp1"' EXIT

awk -v OFS="" \
    -v ie="$iE" -v fe="$fE" -v ic="$i_corband" -v fc="$f_corband" -v se="$stepE" \
    -v pb="$Polar_band" -v qpdir="$QPdir" \
    -v bex="$BSE_Ex" -v bblk="$BSE_Blk" -v fqk="$f_qk" -v nbx="$Block_size" -v xqk="$Xqk" \
    -v lkd="$lkd" '
BEGIN{
  printed_em1d=0; printed_ppa=0; printed_chimod=0;
  appended_XfnQP=0; appended_parallel=0;
  wrote_allgexx=0; wrote_lkind=0;
}

function print_XfnQP_block(){
  if (appended_XfnQP) return
  print "XfnQPdb= \"none\"                  # [EXTQP Xd] Database action"
  print "XfnQP_INTERP_NN= 1               # [EXTQP Xd] Interpolation neighbours (NN mode)"
  print "XfnQP_INTERP_shells= 20.00000    # [EXTQP Xd] Interpolation shells (BOLTZ mode)"
  print "XfnQP_DbGd_INTERP_mode= \"NN\"     # [EXTQP Xd] Interpolation DbGd mode"
  print "% XfnQP_E"
  print " 0.000000 | 1.000000 | 1.000000 |        # [EXTQP Xd] E parameters  (c/v) eV|adim|adim"
  print "%"
  print "XfnQP_Z= ( 1.000000 , 0.000000 )         # [EXTQP Xd] Z factor  (c/v)"
  print "XfnQP_Wv_E= 0.000000       eV    # [EXTQP Xd] W Energy reference  (valence)"
  print "% XfnQP_Wv"
  print " 0.000000 | 0.000000 | 0.000000 |        # [EXTQP Xd] W parameters  (valence) eV| 1|eV^-1"
  print "%"
  print "XfnQP_Wv_dos= 0.000000     eV    # [EXTQP Xd] W dos pre-factor  (valence)"
  print "XfnQP_Wc_E= 0.000000       eV    # [EXTQP Xd] W Energy reference  (conduction)"
  print "% XfnQP_Wc"
  print " 0.000000 | 0.000000 | 0.000000 |        # [EXTQP Xd] W parameters  (conduction) eV| 1 |eV^-1"
  print "%"
  print "XfnQP_Wc_dos= 0.000000     eV    # [EXTQP Xd] W dos pre-factor  (conduction)"
  appended_XfnQP=1
}

function print_Xp_controls(){
  print "% BndsRnXp"
  printf "   1 |  %s |                         # [Xp] Polarization function bands\n", pb
  print "%"
  printf "NGsBlkXp=  %s mHa                 # [Xp] Response block size\n", nbx
  print "% LongDrXp"
  print " 0.100000E-4 |  0.00000    |  0.00000    # [Xp] [cc] Electric Field"
  print "%"
  print "PPAPntXp= 27.21138         eV    # [Xp] PPA imaginary energy"
}

function print_parallel_hints(){
  print "BS_ROLEs = \"eh.k\"            # [PARALLEL] MPI roles for the BSE kernel"
  print "BS_CPU   = \"4.4\"             # [PARALLEL] 4 over e-h, 6 over k = 24 MPI ranks"
  print "BS_nCPU_invert = 4            # [PARALLEL] CPUs for matrix inversion"
  print "BS_nCPU_diago  = 4            # [PARALLEL] CPUs for matrix diagonalization"
}


{
  line=$0
  if (match(line,/^[[:space:]]*#/) || match(line,/^[[:space:]]*$/)) { print line; next }
  clean=line; sub(/[[:space:]]*#.*$/,"",clean)

  # Insert em1d + ppa + Chimod after dipoles
  if (clean ~ /^dipoles[[:space:]]*$/){
    print line
    if (!printed_em1d){ print "em1d                             # [R][X] Dynamically Screened Interaction"; printed_em1d=1 }
    if (!printed_ppa){  print "ppa                              # [R][Xp] Plasmon Pole Approximation for the Screened Interaction"; printed_ppa=1 }
    if (!printed_chimod){ print "Chimod= \"HARTREE\"                # [X] IP/Hartree/ALDA/LRC/PF/BSfxc"; printed_chimod=1 }
    next
  }

  # Remove EXCTemp
  if (clean ~ /^EXCTemp[[:space:]]*=/) next

  # --- Replace BSENGexx with value in mHa (requested) ---
  if (clean ~ /^[[:space:]]*BSENGexx[[:space:]]*=/){
    printf "BSENGexx=  %s            mHa    # [BSK] Exchange components\n", bex
    next
  }

  # --- (OPTIONAL) Force ALLGexx instead of BSENGexx ---
  # NOTE: This block is commented out per your request. Uncomment to force ALLGexx.
  # if (clean ~ /^[[:space:]]*BSENGexx[[:space:]]*=/) {
  #   if (!wrote_allgexx) {
  #     print "ALLGexx                      # [BSS] Force the use use all RL vectors for the exchange part"
  #     wrote_allgexx=1
  #   }
  #   next
  # }
  # if (clean ~ /^[[:space:]]*#?[[:space:]]*ALLGexx([[:space:]]|$)/) {
  #   if (!wrote_allgexx) {
  #     print "ALLGexx                      # [BSS] Force the use use all RL vectors for the exchange part"
  #     wrote_allgexx=1
  #   }
  #   next
  # }

  # Lkind -> set to desired value
  if (clean ~ /^[[:space:]]*Lkind[[:space:]]*=/) {
    printf "Lkind=\"%s\"                  #[BSE,X] bar(default)/full/tilde\n", lkd
    wrote_lkind=1
    next
  }

  # % BLongDir -> 1 | 1 | 1
  if ($1=="%" && $2=="BLongDir") {
    print
    getline
    print " 1.000000 | 1.000000 | 1.000000 |        # [BSS] [cc] Electric Field versor"
    next
  }

  # BSENGBlk -> BSE_Blk in mHa
  if (clean ~ /^BSENGBlk[[:space:]]*=/){
    printf "BSENGBlk=  %s        mHa    # [BSK] Screened interaction block size [if -1 uses all the G-vectors of W(q,G,Gp)]\n", bblk
    next
  }

  # NGsBlkXp -> Block_size in mHa
  if (clean ~ /^NGsBlkXp[[:space:]]*=/){
    printf "NGsBlkXp=  %s mHa                 # [Xp] Response block size\n", nbx
    next
  }

  # KfnQPdb
  if (clean ~ /^[[:space:]]*KfnQPdb[[:space:]]*=/){
    printf "KfnQPdb= \"%s\"                  # [EXTQP BSK BSS] Database action\n", qpdir
    next
  }

  # % BSEBands -> i_corband | f_corband |
  if ($1=="%" && $2=="BSEBands") {
    print line; getline
    printf "   %s |  %s |                         # [BSK] Bands range\n", ic, fc
    next
  }

  # % BEnRange -> iE | fE | eV
  if ($1=="%" && $2=="BEnRange") {
    print line; getline
    printf "  %s | %s |   eV    # [BSS] Energy range\n", ie, fe
    next
  }

  # BEnSteps -> stepE
  if (clean ~ /^BEnSteps[[:space:]]*=/){
    printf "BEnSteps= %s                    # [BSS] Energy steps\n", se
    next
  }

  # % BSEQptR -> 1 | f_qk |
  if ($1=="%" && $2=="BSEQptR") {
    print line; getline
    printf " 1 | %s |                             # [BSK] Transferred momenta range\n", fqk
    getline; print "%"
    next
  }

  # % BndsRnXp -> 1 | Polar_band | (if present)
  if ($1=="%" && $2=="BndsRnXp") {
    print line; getline
    printf "   1 |  %s |                         # [Xp] Polarization function bands\n", pb
    print "%"
    next
  }

  print line
}

END{
  # Do NOT append ALLGexx automatically (requested). If you later want it, uncomment next 3 lines.
  # if (!wrote_allgexx) {
  #   print "ALLGexx                      # [BSS] Force the use use all RL vectors for the exchange part"
  # }

  if (!wrote_lkind) {
    printf "Lkind=\"%s\"                  #[BSE,X] bar(default)/full/tilde\n", lkd
  }

  print_XfnQP_block()
  print_Xp_controls()
  print_parallel_hints()
}
' "$infile" > "$tmp1"

# Squeeze extra blank lines
awk '
  BEGIN{blank=0}
  /^[[:space:]]*$/ { if (!blank) { print ""; blank=1 } next }
  { print; blank=0 }
' "$tmp1" > "$outfile"

echo "Wrote: $outfile"