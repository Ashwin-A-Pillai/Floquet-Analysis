#!/usr/bin/env bash
#
###############################################################################
# Defaults with override behavior
###############################################################################

# Static screening
: "${xp_band_i:=1}"
: "${xp_band_f:=15}"
: "${ngs_blk_xp:=3000}"

# Nonlinear-response bands
: "${i_cb:=4}"
: "${f_cb:=5}"

# Nonlinear-response settings
: "${nl_time:=-1.000000}"
: "${N_intg:=CRANKNIC}"
: "${N_corr:=SEX}"

: "${nl_en_min:=0.500000}"
: "${nl_en_max:=8.000000}"
: "${nl_esteps:=12}"
: "${nl_damp:=0.200000}"

# Hartree / exchange / correlation cutoffs
: "${har_rlvcs:=3000}"
: "${exx_rlvcs:=3000}"
: "${cor_rlvcs:=3000}"

# External-QP scissor
: "${Extqp_E:=3.000000}"

# Field 1
: "${efield1_x:=1.000000}"
: "${efield1_y:=0.000000}"
: "${efield1_z:=0.000000}"
: "${F_knd:=SIN}"

export xp_band_i xp_band_f ngs_blk_xp \
       i_cb f_cb \
       nl_time N_intg N_corr \
       nl_en_min nl_en_max nl_esteps nl_damp \
       har_rlvcs exx_rlvcs cor_rlvcs \
       Extqp_E \
       efield1_x efield1_y efield1_z F_knd

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

awk -v xbi="$xp_band_i" \
    -v xbf="$xp_band_f" \
    -v ngs="$ngs_blk_xp" \
    -v ib="$i_cb" \
    -v fb="$f_cb" \
    -v nlt="$nl_time" \
    -v nli="$N_intg" \
    -v nlc="$N_corr" \
    -v emin="$nl_en_min" \
    -v emax="$nl_en_max" \
    -v esteps="$nl_esteps" \
    -v damp="$nl_damp" \
    -v har="$har_rlvcs" \
    -v exx="$exx_rlvcs" \
    -v corr="$cor_rlvcs" \
    -v extqp="$Extqp_E" \
    -v ex="$efield1_x" \
    -v ey="$efield1_y" \
    -v ez="$efield1_z" \
    -v fkind="$F_knd" '
{
  ###########################################################################
  # 1) BndsRnXs block
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
  # 2) NGsBlkXs
  #
  #    NGsBlkXs= 3000 mHa
  ###########################################################################
  if ($1 == "NGsBlkXs=" || $1 ~ /^NGsBlkXs=/) {
    printf "NGsBlkXs= %s       mHa    # [Xs] Response block size\n", ngs
    next
  }

  ###########################################################################
  # 3) LongDrXs block
  #
  #    % LongDrXs
  #      1.000000 | 0.000000 | 0.000000 |
  ###########################################################################
  if (($1 == "%" && $2 == "LongDrXs") || $1 == "%LongDrXs") {
    print

    if (getline > 0) {
      printf " %s | %s | %s |        # [Xs] [cc] Electric Field\n", \
             ex, ey, ez
    }

    next
  }

  ###########################################################################
  # 4) NLBands block
  #
  #    % NLBands
  #      4 | 5 |
  ###########################################################################
  if (($1 == "%" && $2 == "NLBands") || $1 == "%NLBands") {
    print

    if (getline > 0) {
      printf "  %s | %s |                           # [NL] Bands range\n", \
             ib, fb
    }

    next
  }

  ###########################################################################
  # 5) NLtime
  ###########################################################################
  if ($1 == "NLtime=" || $1 ~ /^NLtime=/) {
    printf "NLtime=%s           fs    # [NL] Simulation Time\n", nlt
    next
  }

  ###########################################################################
  # 6) NLintegrator
  ###########################################################################
  if ($1 == "NLintegrator=" || $1 ~ /^NLintegrator=/) {
    printf "NLintegrator= \"%s\"           # [NL] Integrator (\"EULEREXP/RK2/RK4/RK2EXP/HEUN/INVINT/CRANKNIC\")\n", \
           nli
    next
  }

  ###########################################################################
  # 7) NLCorrelation
  ###########################################################################
  if ($1 == "NLCorrelation=" || $1 ~ /^NLCorrelation=/) {
    printf "NLCorrelation= \"%s\"             # [NL] Correlation (\"IPA/HARTREE/TDDFT/LRC/LRW/JGM/SEX/LSEX/LHF\")\n", \
           nlc
    next
  }

  ###########################################################################
  # 8) NLEnRange block
  #
  #    % NLEnRange
  #      0.500000 | 8.000000 | eV
  ###########################################################################
  if (($1 == "%" && $2 == "NLEnRange") || $1 == "%NLEnRange") {
    print

    if (getline > 0) {
      printf " %s | %s |         eV    # [NL] Energy range (for loop on frequencies NLEnSteps/=0\n", \
             emin, emax
    }

    next
  }

  ###########################################################################
  # 9) NLEnSteps
  ###########################################################################
  if ($1 == "NLEnSteps=" || $1 ~ /^NLEnSteps=/) {
    printf "NLEnSteps=%s                      # [NL] Energy steps for the loop on frequencies\n", \
           esteps
    next
  }

  ###########################################################################
  # 10) NLDamping
  ###########################################################################
  if ($1 == "NLDamping=" || $1 ~ /^NLDamping=/) {
    printf "NLDamping= %s        eV    # [NL] Damping (or dephasing)\n", damp
    next
  }

  ###########################################################################
  # 11) HARRLvcs
  ###########################################################################
  if ($1 == "HARRLvcs=" || $1 ~ /^HARRLvcs=/) {
    printf "HARRLvcs= %s            mHa    # [HA] Hartree     RL components\n", har
    next
  }

  ###########################################################################
  # 12) EXXRLvcs
  ###########################################################################
  if ($1 == "EXXRLvcs=" || $1 ~ /^EXXRLvcs=/) {
    printf "EXXRLvcs= %s            mHa    # [XX] Exchange    RL components\n", exx
    next
  }

  ###########################################################################
  # 13) CORRLvcs
  ###########################################################################
  if ($1 == "CORRLvcs=" || $1 ~ /^CORRLvcs=/) {
    printf "CORRLvcs= %s                mHa    # [GW] Correlation RL components [if -1 uses all the G-vectors of W_{G,Gp}]\n", \
           corr
    next
  }

  ###########################################################################
  # 14) GfnQP_E block
  #
  #     % GfnQP_E
  #       3.000000 | 1.000000 | 1.000000 |
  #
  # Only the first component is changed.
  ###########################################################################
  if (($1 == "%" && $2 == "GfnQP_E") || $1 == "%GfnQP_E") {
    print

    if (getline > 0) {
      printf " %s | 1.000000 | 1.000000 |        # [EXTQP G] E parameters  (c/v) eV|adim|adim\n", \
             extqp
    }

    next
  }

  ###########################################################################
  # 15) Field1_kind
  ###########################################################################
  if ($1 == "Field1_kind=" || $1 ~ /^Field1_kind=/) {
    printf "Field1_kind= \"%s\"           # [RT Field1] Kind(SIN|COS|RES|ANTIRES|GAUSS|DELTA|QSSIN)\n", \
           fkind
    next
  }

  ###########################################################################
  # 16) Field1_Dir block
  #
  #     % Field1_Dir
  #       1.000000 | 0.000000 | 0.000000 |
  ###########################################################################
  if (($1 == "%" && $2 == "Field1_Dir") || $1 == "%Field1_Dir") {
    print

    if (getline > 0) {
      printf " %s | %s | %s |        # [RT Field1] Versor\n", \
             ex, ey, ez
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
echo "  XsBands        = $xp_band_i $xp_band_f"
echo "  NGsBlkXs       = $ngs_blk_xp mHa"
echo "  NLBands        = $i_cb $f_cb"
echo "  nl_time        = $nl_time fs"
echo "  nl_integrator  = $N_intg"
echo "  nl_correlation = $N_corr"
echo "  nl_en_range    = $nl_en_min $nl_en_max eV"
echo "  nl_esteps      = $nl_esteps"
echo "  nl_damp        = $nl_damp eV"
echo "  HARRLvcs       = $har_rlvcs mHa"
echo "  EXXRLvcs       = $exx_rlvcs mHa"
echo "  CORRLvcs       = $cor_rlvcs mHa"
echo "  Extqp_E        = $Extqp_E eV"
echo "  field1_kind    = $F_knd"
echo "  field1_dir     = $efield1_x $efield1_y $efield1_z"