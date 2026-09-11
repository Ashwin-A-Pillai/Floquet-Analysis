#!/usr/bin/env bash
#
# patch_NLR.sh
#
# Patch ONLY existing values in a Yambo/Lumen nonlinear-optics input.
# No variables or blocks are added if they are absent.
#
# Generate template with:
#   yambo_nl -u n -V all -Input ynl_temp.in
#
# Usage:
#   ./patch_NLR.sh ynl_temp.in YNL.in
#

set -euo pipefail

###############################################################################
# Defaults with override behavior
###############################################################################

# Parallelization
: "${NL_CPU_LAYOUT:=10 1}"
: "${NL_ROLES:=w k}"

# Bands
: "${nl_band_i:=3}"
: "${nl_band_f:=6}"

# Time propagation
: "${nl_time:=-1.000000}"
: "${nl_integrator:=INVINT}"
: "${nl_correlation:=IPA}"
: "${nl_lrc_alpha:=0.000000}"

# Frequency range
: "${nl_en_min:=1.000000}"
: "${nl_en_max:=5.000000}"
: "${nl_esteps:=10}"

# Damping
: "${nl_damp:=0.150000}"

# Field frequency
: "${field1_freq_i:=0.100000}"
: "${field1_freq_f:=0.100000}"
: "${field1_nfreqs:=1}"

# Field
: "${field1_int:=1000.00}"
: "${field1_width:=0.000000}"
: "${field1_kind:=SOFTSIN}"
: "${field1_pol:=linear}"

# Field direction
: "${efield1_x:=1.000000}"
: "${efield1_y:=1.000000}"
: "${efield1_z:=0.000000}"

# QP database
: "${gfnqpdb:=E < SAVE/ndb.QP}"

export NL_CPU_LAYOUT NL_ROLES \
       nl_band_i nl_band_f \
       nl_time nl_integrator nl_correlation nl_lrc_alpha \
       nl_en_min nl_en_max nl_esteps nl_damp \
       field1_freq_i field1_freq_f field1_nfreqs \
       field1_int field1_width field1_kind field1_pol \
       efield1_x efield1_y efield1_z \
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
    -v bi="$nl_band_i" \
    -v bf="$nl_band_f" \
    -v nltime="$nl_time" \
    -v nlint="$nl_integrator" \
    -v nlcorr="$nl_correlation" \
    -v nlalpha="$nl_lrc_alpha" \
    -v emin="$nl_en_min" \
    -v emax="$nl_en_max" \
    -v esteps="$nl_esteps" \
    -v damp="$nl_damp" \
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

    if ($0 !~ /^[[:space:]]*#/ && $0 !~ /^[[:space:]]*$/) {
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

    if ($0 !~ /^[[:space:]]*#/ && $0 !~ /^[[:space:]]*$/) {
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
# Field1_Freq
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

    if ($0 !~ /^[[:space:]]*#/ && $0 !~ /^[[:space:]]*$/) {
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

    if ($0 !~ /^[[:space:]]*#/ && $0 !~ /^[[:space:]]*$/) {
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
# Scalar variables
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

/^[[:space:]]*NLtime[[:space:]]*=/ {
    c = get_comment($0)
    printf "NLtime=%s           fs", nltime
    if (c != "")
        printf "    %s", c
    printf "\n"
    next
}

/^[[:space:]]*NLintegrator[[:space:]]*=/ {
    c = get_comment($0)
    printf "NLintegrator= \"%s\"", nlint
    if (c != "")
        printf "           %s", c
    printf "\n"
    next
}

/^[[:space:]]*NLCorrelation[[:space:]]*=/ {
    c = get_comment($0)
    printf "NLCorrelation= \"%s\"", nlcorr
    if (c != "")
        printf "             %s", c
    printf "\n"
    next
}

/^[[:space:]]*NLLrcAlpha[[:space:]]*=/ {
    c = get_comment($0)
    printf "NLLrcAlpha= %s", nlalpha
    if (c != "")
        printf "             %s", c
    printf "\n"
    next
}

/^[[:space:]]*NLEnSteps[[:space:]]*=/ {
    c = get_comment($0)
    printf "NLEnSteps=%s", esteps
    if (c != "")
        printf "                     %s", c
    printf "\n"
    next
}

/^[[:space:]]*NLDamping[[:space:]]*=/ {
    c = get_comment($0)
    printf "NLDamping= %s        eV", damp
    if (c != "")
        printf "    %s", c
    printf "\n"
    next
}

/^[[:space:]]*Field1_NFreqs[[:space:]]*=/ {
    c = get_comment($0)
    printf "Field1_NFreqs= %s", fnfreq
    if (c != "")
        printf "                 %s", c
    printf "\n"
    next
}

/^[[:space:]]*Field1_Int[[:space:]]*=/ {
    c = get_comment($0)
    printf "Field1_Int= %s       kWLm2", fint
    if (c != "")
        printf "    %s", c
    printf "\n"
    next
}

/^[[:space:]]*Field1_Width[[:space:]]*=/ {
    c = get_comment($0)
    printf "Field1_Width= %s     fs", fwidth
    if (c != "")
        printf "    %s", c
    printf "\n"
    next
}

/^[[:space:]]*Field1_kind[[:space:]]*=/ {
    c = get_comment($0)
    printf "Field1_kind= \"%s\"", fkind
    if (c != "")
        printf "           %s", c
    printf "\n"
    next
}

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
# Since the template is generated with -V all, GfnQPdb may be present.
# Patch it ONLY if it actually exists.
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

echo "Patched nonlinear input written to: $outfile"
echo "  NLBands       = ${nl_band_i}-${nl_band_f}"
echo "  NLEnRange     = ${nl_en_min}-${nl_en_max} eV"
echo "  NLEnSteps     = ${nl_esteps}"
echo "  NLDamping     = ${nl_damp} eV"
echo "  NL_CPU        = ${NL_CPU_LAYOUT}"
echo "  NL_ROLEs      = ${NL_ROLES}"
echo "  Field1_Freq   = ${field1_freq_i}-${field1_freq_f} eV"
echo "  Field1_Dir    = ${efield1_x} ${efield1_y} ${efield1_z}"