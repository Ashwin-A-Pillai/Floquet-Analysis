#!/usr/bin/env bash
#
# patch_NLR.sh
#
# Patch ONLY existing values in a Yambo/Lumen nonlinear-optics input.
# No variables or blocks are added if they are absent.
#
# Compatible with, e.g.:
#
#   yambo_nl -u n -V all -Input ynl_temp.in
#   yambo_nl -u p -V nl  -Input ynl_temp.in
#
# Usage:
#
#   ./patch_NLR.sh ynl_temp.in YNL.in
#
# Environment variables can be used to override all defaults below.
#

set -euo pipefail

###############################################################################
# Defaults with override behavior
###############################################################################

# -----------------------------------------------------------------------------
# Parallelization
#
# Example for:
#
#   #SBATCH --ntasks=24
#   #SBATCH --cpus-per-task=2
#
# Pump/Floquet calculation with one external-field frequency:
#
#   NL_CPU   = "1 24"
#   NL_ROLEs = "w k"
#
# Yambo can take OpenMP threading from OMP_NUM_THREADS when NL_Threads=0.
# -----------------------------------------------------------------------------

: "${NL_CPU_LAYOUT:=1 ${SLURM_NTASKS:-24}}"
: "${NL_ROLES:=w k}"
: "${NL_THREADS:=0}"

# -----------------------------------------------------------------------------
# Bands
# -----------------------------------------------------------------------------

: "${nl_band_i:=3}"
: "${nl_band_f:=6}"

# -----------------------------------------------------------------------------
# Time propagation
# -----------------------------------------------------------------------------

: "${nl_step:=0.010000}"
: "${nl_time:=-1.000000}"
: "${nl_integrator:=INVINT}"
: "${nl_correlation:=IPA}"
: "${nl_lrc_alpha:=0.000000}"

# -----------------------------------------------------------------------------
# NLStepDiv
#
# true  -> uncomment NLStepDiv
# false -> comment NLStepDiv
#
# Accepted true values:
#   true, 1, yes, on
#
# Accepted false values:
#   false, 0, no, off
# -----------------------------------------------------------------------------

: "${NTDiv:=true}"

# -----------------------------------------------------------------------------
# Floquet analysis
#
# Example:
#   -4  -> Floquet analysis up to harmonic order 4
#   -1  -> disabled/default
# -----------------------------------------------------------------------------

: "${fl_order:=-4}"

# -----------------------------------------------------------------------------
# Frequency range for ordinary nonlinear-frequency-loop calculations
# -----------------------------------------------------------------------------

: "${nl_en_min:=1.000000}"
: "${nl_en_max:=5.000000}"
: "${nl_esteps:=10}"

# -----------------------------------------------------------------------------
# Damping
# -----------------------------------------------------------------------------

: "${nl_damp:=0.150000}"

# -----------------------------------------------------------------------------
# Scalar pump-mode field frequency
#
# Used for inputs generated with:
#
#   yambo_nl -u p ...
#
# which contain:
#
#   Field1_Freq= 4.100000 eV
# -----------------------------------------------------------------------------

: "${field1_freq:=0.100000}"

# -----------------------------------------------------------------------------
# Field-frequency RANGE
#
# Retained for templates containing:
#
#   % Field1_Freq
#      ...
#   %
# -----------------------------------------------------------------------------

: "${field1_freq_i:=0.100000}"
: "${field1_freq_f:=0.100000}"
: "${field1_nfreqs:=1}"

# -----------------------------------------------------------------------------
# Field parameters
# -----------------------------------------------------------------------------

: "${field1_int:=1000.00}"
: "${field1_width:=0.000000}"
: "${field1_kind:=SOFTSIN}"
: "${field1_pol:=linear}"

# -----------------------------------------------------------------------------
# Field direction
# -----------------------------------------------------------------------------

: "${efield1_x:=1.000000}"
: "${efield1_y:=1.000000}"
: "${efield1_z:=0.000000}"

# -----------------------------------------------------------------------------
# QP database
# -----------------------------------------------------------------------------

: "${gfnqpdb:=E < SAVE/ndb.QP}"

###############################################################################
# Normalize boolean switches
###############################################################################

case "${NTDiv,,}" in

    true|1|yes|on)
        NTDiv_BOOL=1
        ;;

    false|0|no|off)
        NTDiv_BOOL=0
        ;;

    *)
        echo "[ERROR] Invalid NTDiv value: ${NTDiv}" >&2
        echo "        Allowed values:" >&2
        echo "          true / false" >&2
        echo "          1 / 0" >&2
        echo "          yes / no" >&2
        echo "          on / off" >&2
        exit 3
        ;;

esac

###############################################################################
# Export values
###############################################################################

export \
    NL_CPU_LAYOUT \
    NL_ROLES \
    NL_THREADS \
    nl_band_i \
    nl_band_f \
    nl_step \
    nl_time \
    nl_integrator \
    nl_correlation \
    nl_lrc_alpha \
    fl_order \
    nl_en_min \
    nl_en_max \
    nl_esteps \
    nl_damp \
    field1_freq \
    field1_freq_i \
    field1_freq_f \
    field1_nfreqs \
    field1_int \
    field1_width \
    field1_kind \
    field1_pol \
    efield1_x \
    efield1_y \
    efield1_z \
    gfnqpdb

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
# Temporary output
###############################################################################

tmp="$(mktemp)"

trap 'rm -f "$tmp"' EXIT

###############################################################################
# Patch existing variables only
###############################################################################

awk \
    -v nlcpu="$NL_CPU_LAYOUT" \
    -v nlroles="$NL_ROLES" \
    -v nlthreads="$NL_THREADS" \
    -v bi="$nl_band_i" \
    -v bf="$nl_band_f" \
    -v nlstep="$nl_step" \
    -v nltime="$nl_time" \
    -v nlint="$nl_integrator" \
    -v nlcorr="$nl_correlation" \
    -v nlalpha="$nl_lrc_alpha" \
    -v ntdiv="$NTDiv_BOOL" \
    -v florder="$fl_order" \
    -v emin="$nl_en_min" \
    -v emax="$nl_en_max" \
    -v esteps="$nl_esteps" \
    -v damp="$nl_damp" \
    -v ff="$field1_freq" \
    -v ffi="$field1_freq_i" \
    -v fff="$field1_freq_f" \
    -v fnfreq="$field1_nfreqs" \
    -v fint="$field1_int" \
    -v fwidth="$field1_width" \
    -v fkind="$field1_kind" \
    -v fpol="$field1_pol" \
    -v ex="$efield1_x" \
    -v ey="$efield1_y" \
    -v ez="$efield1_z" \
    -v qpdb="$gfnqpdb" '

###############################################################################
# Helper
###############################################################################

function get_comment(s,    p) {
    p = index(s, "#")

    if (p > 0)
        return substr(s, p)

    return ""
}

BEGIN {
    block = ""
}

###############################################################################
# NLBands
###############################################################################

/^[[:space:]]*%[[:space:]]*NLBands[[:space:]]*$/ {

    print
    block = "NLBands"
    next
}

block == "NLBands" {

    if ($0 ~ /^[[:space:]]*%[[:space:]]*$/) {

        print
        block = ""
        next
    }

    if ($0 !~ /^[[:space:]]*#/ &&
        $0 !~ /^[[:space:]]*$/) {

        c = get_comment($0)

        printf "  %s | %s |", bi, bf

        if (c != "")
            printf "                           %s", c

        printf "\n"
        next
    }

    print
    next
}

###############################################################################
# NLEnRange
###############################################################################

/^[[:space:]]*%[[:space:]]*NLEnRange[[:space:]]*$/ {

    print
    block = "NLEnRange"
    next
}

block == "NLEnRange" {

    if ($0 ~ /^[[:space:]]*%[[:space:]]*$/) {

        print
        block = ""
        next
    }

    if ($0 !~ /^[[:space:]]*#/ &&
        $0 !~ /^[[:space:]]*$/) {

        c = get_comment($0)

        printf " %s | %s |         eV", emin, emax

        if (c != "")
            printf "    %s", c

        printf "\n"
        next
    }

    print
    next
}

###############################################################################
# Field1_Freq BLOCK
#
# Used by some nonlinear-frequency-loop input templates:
#
# % Field1_Freq
#  1.0 | 5.0 | eV
# %
###############################################################################

/^[[:space:]]*%[[:space:]]*Field1_Freq[[:space:]]*$/ {

    print
    block = "Field1_Freq"
    next
}

block == "Field1_Freq" {

    if ($0 ~ /^[[:space:]]*%[[:space:]]*$/) {

        print
        block = ""
        next
    }

    if ($0 !~ /^[[:space:]]*#/ &&
        $0 !~ /^[[:space:]]*$/) {

        c = get_comment($0)

        printf " %s | %s |         eV", ffi, fff

        if (c != "")
            printf "    %s", c

        printf "\n"
        next
    }

    print
    next
}

###############################################################################
# Field1_Dir
###############################################################################

/^[[:space:]]*%[[:space:]]*Field1_Dir[[:space:]]*$/ {

    print
    block = "Field1_Dir"
    next
}

block == "Field1_Dir" {

    if ($0 ~ /^[[:space:]]*%[[:space:]]*$/) {

        print
        block = ""
        next
    }

    if ($0 !~ /^[[:space:]]*#/ &&
        $0 !~ /^[[:space:]]*$/) {

        c = get_comment($0)

        printf " %s | %s | %s |", ex, ey, ez

        if (c != "")
            printf "        %s", c

        printf "\n"
        next
    }

    print
    next
}

###############################################################################
# Parallelization
###############################################################################

/^[[:space:]]*NL_CPU[[:space:]]*=/ {

    c = get_comment($0)

    printf "NL_CPU= \"%s\"", nlcpu

    if (c != "")
        printf "                       %s", c

    printf "\n"
    next
}

/^[[:space:]]*NL_ROLEs[[:space:]]*=/ {

    c = get_comment($0)

    printf "NL_ROLEs= \"%s\"", nlroles

    if (c != "")
        printf "                     %s", c

    printf "\n"
    next
}

/^[[:space:]]*NL_Threads[[:space:]]*=/ {

    c = get_comment($0)

    printf "NL_Threads=%s", nlthreads

    if (c != "")
        printf "                    %s", c

    printf "\n"
    next
}

###############################################################################
# NLstep
###############################################################################

/^[[:space:]]*NLstep[[:space:]]*=/ {

    c = get_comment($0)

    printf "NLstep= %s           fs", nlstep

    if (c != "")
        printf "    %s", c

    printf "\n"
    next
}

###############################################################################
# NLStepDiv
#
# NTDiv=true:
#
#   #NLStepDiv   ->   NLStepDiv
#
# NTDiv=false:
#
#   NLStepDiv    ->   #NLStepDiv
#
# The line is modified ONLY if NLStepDiv already exists in the template.
###############################################################################

/^[[:space:]]*#?[[:space:]]*NLStepDiv([[:space:]]|$)/ {

    line = $0

    # Remove leading whitespace, optional #, keyword, and whitespace
    # immediately following the keyword.
    sub(/^[[:space:]]*#?[[:space:]]*NLStepDiv[[:space:]]*/, "", line)

    if (ntdiv == 1)
        printf "NLStepDiv"
    else
        printf "#NLStepDiv"

    if (line != "")
        printf "                     %s", line

    printf "\n"
    next
}

###############################################################################
# NLtime
###############################################################################

/^[[:space:]]*NLtime[[:space:]]*=/ {

    c = get_comment($0)

    printf "NLtime=%s           fs", nltime

    if (c != "")
        printf "    %s", c

    printf "\n"
    next
}

###############################################################################
# NLintegrator
###############################################################################

/^[[:space:]]*NLintegrator[[:space:]]*=/ {

    c = get_comment($0)

    printf "NLintegrator= \"%s\"", nlint

    if (c != "")
        printf "           %s", c

    printf "\n"
    next
}

###############################################################################
# NLCorrelation
###############################################################################

/^[[:space:]]*NLCorrelation[[:space:]]*=/ {

    c = get_comment($0)

    printf "NLCorrelation= \"%s\"", nlcorr

    if (c != "")
        printf "             %s", c

    printf "\n"
    next
}

###############################################################################
# NLLrcAlpha
###############################################################################

/^[[:space:]]*NLLrcAlpha[[:space:]]*=/ {

    c = get_comment($0)

    printf "NLLrcAlpha= %s", nlalpha

    if (c != "")
        printf "             %s", c

    printf "\n"
    next
}

###############################################################################
# FLOrder
###############################################################################

/^[[:space:]]*FLOrder[[:space:]]*=/ {

    c = get_comment($0)

    printf "FLOrder=%s", florder

    if (c != "")
        printf "                       %s", c

    printf "\n"
    next
}

###############################################################################
# NLEnSteps
###############################################################################

/^[[:space:]]*NLEnSteps[[:space:]]*=/ {

    c = get_comment($0)

    printf "NLEnSteps=%s", esteps

    if (c != "")
        printf "                     %s", c

    printf "\n"
    next
}

###############################################################################
# NLDamping
###############################################################################

/^[[:space:]]*NLDamping[[:space:]]*=/ {

    c = get_comment($0)

    printf "NLDamping= %s        eV", damp

    if (c != "")
        printf "    %s", c

    printf "\n"
    next
}

###############################################################################
# Scalar Field1_Freq
#
# Used by pump/Floquet mode:
#
#   Field1_Freq= 4.100000 eV
###############################################################################

/^[[:space:]]*Field1_Freq[[:space:]]*=/ {

    c = get_comment($0)

    printf "Field1_Freq= %s      eV", ff

    if (c != "")
        printf "    %s", c

    printf "\n"
    next
}

###############################################################################
# Field1_NFreqs
###############################################################################

/^[[:space:]]*Field1_NFreqs[[:space:]]*=/ {

    c = get_comment($0)

    printf "Field1_NFreqs= %s", fnfreq

    if (c != "")
        printf "                 %s", c

    printf "\n"
    next
}

###############################################################################
# Field1_Int
###############################################################################

/^[[:space:]]*Field1_Int[[:space:]]*=/ {

    c = get_comment($0)

    printf "Field1_Int= %s       kWLm2", fint

    if (c != "")
        printf "    %s", c

    printf "\n"
    next
}

###############################################################################
# Field1_Width
###############################################################################

/^[[:space:]]*Field1_Width[[:space:]]*=/ {

    c = get_comment($0)

    printf "Field1_Width= %s     fs", fwidth

    if (c != "")
        printf "    %s", c

    printf "\n"
    next
}

###############################################################################
# Field1_kind
###############################################################################

/^[[:space:]]*Field1_kind[[:space:]]*=/ {

    c = get_comment($0)

    printf "Field1_kind= \"%s\"", fkind

    if (c != "")
        printf "           %s", c

    printf "\n"
    next
}

###############################################################################
# Field1_pol
###############################################################################

/^[[:space:]]*Field1_pol[[:space:]]*=/ {

    c = get_comment($0)

    printf "Field1_pol= \"%s\"", fpol

    if (c != "")
        printf "             %s", c

    printf "\n"
    next
}

###############################################################################
# QP database
#
# Patch ONLY if GfnQPdb is actually present in the generated template.
###############################################################################

/^[[:space:]]*GfnQPdb[[:space:]]*=/ {

    c = get_comment($0)

    printf "GfnQPdb= \"%s\"", qpdb

    if (c != "")
        printf "              %s", c

    printf "\n"
    next
}

###############################################################################
# Everything else is preserved exactly
###############################################################################

{
    print
}

' "$infile" > "$tmp"

###############################################################################
# Install result
###############################################################################

mv "$tmp" "$outfile"

trap - EXIT

###############################################################################
# Summary
###############################################################################

echo
echo "Patched nonlinear input written to: $outfile"
echo
echo "  Parallelization"
echo "    NL_CPU        = ${NL_CPU_LAYOUT}"
echo "    NL_ROLEs      = ${NL_ROLES}"
echo "    NL_Threads    = ${NL_THREADS}"
echo
echo "  Bands"
echo "    NLBands       = ${nl_band_i}-${nl_band_f}"
echo
echo "  Time propagation"
echo "    NLstep        = ${nl_step} fs"
echo "    NLtime        = ${nl_time} fs"
echo "    NLintegrator  = ${nl_integrator}"
echo "    NLCorrelation = ${nl_correlation}"
echo
echo "  Floquet"
echo "    NTDiv         = ${NTDiv}"
echo "    FLOrder       = ${fl_order}"
echo
echo "  NL frequency range"
echo "    NLEnRange     = ${nl_en_min}-${nl_en_max} eV"
echo "    NLEnSteps     = ${nl_esteps}"
echo
echo "  Damping"
echo "    NLDamping     = ${nl_damp} eV"
echo
echo "  Pump field"
echo "    Field1_Freq   = ${field1_freq} eV"
echo "    Field1_Int    = ${field1_int} kWLm2"
echo "    Field1_kind   = ${field1_kind}"
echo "    Field1_pol    = ${field1_pol}"
echo "    Field1_Dir    = ${efield1_x} ${efield1_y} ${efield1_z}"
echo