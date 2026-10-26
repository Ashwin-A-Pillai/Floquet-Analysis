#!/usr/bin/env bash
# patch_BSE.sh: edit a generated Lumen/Yambo input; do not generate or run it.
# Template: 2D: "$YAMBO" -o b -k sex -X s -y d -r -V all -Input BSE_temp.in -Q
# 3D: "$YAMBO" -o b -k sex -y d -V qp -Input BSE_temp.in -Q
# Usage: bash patch_BSE.sh input.in [output.in]
set -euo pipefail

# Defaults / environment overrides
# This patcher can save full BSE eigenvectors for exciton/BZ analysis.
# BSE_WRITE_WF=1 enables WRbsWF (default). BSE_WRITE_WF=0 explicitly disables it.

# Unset band bounds/cutoffs preserve the generated template; no fixed material/grid.
: "${dim:=${DIM:-3D}}" "${cgeo:=}" "${RQpts:=5000000}" "${RGvec:=29}"
: "${efield1_x:=1.0}" "${efield1_y:=0.0}" "${efield1_z:=0.0}"
: "${MPI_RANKS:=${SLURM_STEP_NUM_TASKS:-${SLURM_NTASKS:-1}}}"
: "${OMP_NUM_THREADS:=${SLURM_CPUS_PER_TASK:-1}}"
export dim cgeo RQpts RGvec efield1_x efield1_y efield1_z MPI_RANKS OMP_NUM_THREADS

: "${xp_band_i:=}" "${xp_band_f:=${SCREEN_BANDS:-}}"
: "${ngs_blk_xp:=${SCREEN_CUTOFF_MHA:-}}" "${ngs_blk_xs:=${ngs_blk_xp}}"
export xp_band_i xp_band_f ngs_blk_xp ngs_blk_xs

# Choose one energy model explicitly. Relative paths refer to the run directory.
: "${QP_MODE:=none}" "${SCISSOR_EV:=${Extqp_E:-${KE_V:-0.0}}}" "${QP_DB:=}"
if [[ ${QP_MODE,,} == gwdb ]]; then
    QP_DB=${QP_DB:-${QPdir:-GW0/ndb.QP}}
    QP_DB=${QP_DB#E < }
    [[ -n $QP_DB && $QP_DB != none ]] || { echo "Set QP_DB for QP_MODE=gwdb" >&2; exit 2; }
    if [[ ${CHECK_QP_DB:-1} == 1 && ! -r $QP_DB ]]; then
        echo "[ERROR] Cannot read QP database: $QP_DB (relative to $PWD)" >&2; exit 2
    fi
fi
export QP_MODE SCISSOR_EV QP_DB

: "${i_band:=${i_corband:-}}" "${f_band:=${f_corband:-}}" "${iq:=1}" "${fq:=1}"
: "${BSE_Ex:=${BSE_EX:-}}" "${BSE_Blk:=}" "${BSE_Blk_unit:=mHa}"
: "${iE:=}" "${fE:=}" "${stepE:=}" "${BSE_DAMP:=}"
: "${BSKmod:=SEX}" "${BSEmod:=resonant}" "${BSE_WRITE_WF:=1}"
export BSE_WRITE_WF
[[ $BSE_WRITE_WF == 0 || $BSE_WRITE_WF == 1 ]] || { echo "BSE_WRITE_WF must be 0 or 1" >&2; exit 2; }
export i_band f_band iq fq BSE_Ex BSE_Blk BSE_Blk_unit iE fE stepE BSE_DAMP BSKmod BSEmod

# Export optional decompositions. 'auto'/unset retains Yambo automatic distribution.
export DIP_CPU_LAYOUT DIP_ROLES X_Q0_CPU_LAYOUT X_Q0_ROLES X_ALL_CPU_LAYOUT X_ALL_ROLES
export SE_CPU_LAYOUT SE_ROLES BS_CPU_LAYOUT BS_ROLES X_CPU_LAYOUT X_ROLES COLL_CPU_LAYOUT COLL_ROLES
export DIP_THREADS X_THREADS SE_THREADS BS_THREADS NL_THREADS COLL_THREADS BS_DIAGO_CPUS BS_INV_CPUS

# Arguments and atomic output (input=output is safe).
[[ $# -ge 1 && $# -le 2 ]] || { echo "Usage: $0 input.in [output.in]" >&2; exit 2; }
infile=$1
outfile=${2:-BSE.in}
[[ -r "$infile" ]] || { echo "Unreadable input: $infile" >&2; exit 2; }
mkdir -p -- "$(dirname -- "$outfile")"
tmp=$(mktemp "${outfile}.tmp.XXXXXX")
trap 'rm -f -- "$tmp"' EXIT
export PATCH_STAGE=patch_BSE

# Index and edit complete entries, preserving unrelated generated settings.
awk '

function trim(v) { sub(/^[ \t]+/, "", v); sub(/[ \t]+$/, "", v); return v }
function fail(msg) { print "[" stage "] ERROR: " msg > "/dev/stderr"; bad=1; exit 2 }
function quote(v) { if (v ~ /["\r\n]/) fail("Invalid quoted value: " v); return "\"" v "\"" }
function number(v) { return v ~ /^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eEdD][+-]?[0-9]+)?$/ }
function positive_int(v) { return v ~ /^[0-9]+$/ && v+0 > 0 }
function uncomment(v) { sub(/^[ \t]*#[ \t]?/, "", v); return v }
function code(v,    j,q,c,out) {
    q=0; out=""
    for (j=1;j<=length(v);j++) { c=substr(v,j,1); if(c=="\"") q=!q; if(c=="#" && !q) break; out=out c }
    return trim(out)
}
function locate(k) {
    if (active_count[k]>1) fail("Duplicate active definition of " k)
    return active_index[k] ? active_index[k] : commented_index[k]
}
function has(k) { return locate(k)>0 }
function value(k,    j,v) {
    j=locate(k); if(!j) return ""
    if(kind[j]=="block") { v=""; for(j=j+1;j<finish[locate(k)];j++) v=v " " code(disabled_node[locate(k)]?uncomment(lines[j]):lines[j]); return trim(v) }
    v=code(lines[j]); sub(/^[^=]*=/,"",v); v=trim(v); gsub(/^"|"$/,"",v); return v
}
function put(k,typ,v,policy,    j,e,new) {
    j=locate(k)
    if(!j && policy=="optional") return
    if(!j && policy=="required") fail("Missing " k "; regenerate the correct runlevel/template (use -V all where supported)")
    if(j && kind[j]!=typ) fail("Wrong input type for " k ": found " kind[j] ", expected " typ)
    if(typ=="block") new="% " k "\n " v "\n%"
    else if(typ=="flag") new=(v ? k : "#" k)
    else new=k "= " v
    if(j) { replacement[j]=new; for(e=j+1;e<=finish[j];e++) skip[e]=1 }
    else { appended[++nappend]=new }
    changes++
}
function scalar(k,v,policy) { put(k,"scalar",v,policy) }
function block(k,v,policy) { put(k,"block",v,policy) }
function flag(k,v,policy) { put(k,"flag",v,policy) }
function optional_env(k,env,units,policy,    v) {
    v=ENVIRON[env]; if(v!="") { if(!number(v)) fail(env " must be numeric"); scalar(k,v (units!=""?" " units:""),policy) }
}
function range(k,lo,hi,    a,n,v) {
    if(lo=="" && hi=="") return
    if(!has(k)) fail("Missing % " k " for requested band/range override")
    n=split(value(k),a,"|"); if(n<2) fail("Malformed % " k)
    if(lo=="") lo=trim(a[1]); if(hi=="") hi=trim(a[2])
    if(!positive_int(lo)||!positive_int(hi)||lo+0>hi+0) fail("Invalid range for " k ": " lo ".." hi)
    block(k,lo " | " hi " |","required")
}
function energy_range(k,lo,hi,units,    a,n) {
    if(lo=="" && hi=="") return
    n=split(value(k),a,"|"); if(n<2) fail("Missing/malformed % " k)
    if(lo=="") lo=trim(a[1]); if(hi=="") hi=trim(a[2])
    if(!number(lo)||!number(hi)||hi+0<=lo+0) fail("Invalid energy range for " k)
    block(k,lo " | " hi " | " units,"required")
}
function field(k,policy,    x,y,z) {
    x=ENVIRON["efield1_x"]; y=ENVIRON["efield1_y"]; z=ENVIRON["efield1_z"]
    if(!number(x)||!number(y)||!number(z)||x*x+y*y+z*z==0) fail("Field direction must be a finite nonzero 3-vector")
    block(k,x " | " y " | " z " |",policy)
}
function cutoff(    d,g) {
    d=toupper(ENVIRON["dim"]); if(d!="2D" && d!="3D") fail("dim must be 2D or 3D")
    g=ENVIRON["cgeo"]; if(g=="") g=(d=="2D"?"slab z":"none")
    if(d=="3D" && tolower(g)!="none") fail("3D bulk requested with a Coulomb cutoff; use cgeo=none")
    scalar("CUTGeo",quote(g),(d=="2D" && tolower(g)!="none")?"required":"optional")
    optional_env("RandQpts","RQpts","","optional"); optional_env("RandGvec","RGvec","RL","optional")
}
function screening(    found) {
    found=0
    if(has("BndsRnXs")) {
        range("BndsRnXs",ENVIRON["xp_band_i"],ENVIRON["xp_band_f"])
        optional_env("NGsBlkXs","ngs_blk_xs","mHa","required"); field("LongDrXs","optional"); found=1
    }
    if(has("BndsRnXp")) {
        range("BndsRnXp",ENVIRON["xp_band_i"],ENVIRON["xp_band_f"])
        optional_env("NGsBlkXp","ngs_blk_xp","mHa","required"); field("LongDrXp","optional"); found=1
    }
    if(!found && (ENVIRON["xp_band_f"]!="" || ENVIRON["ngs_blk_xp"]!="" || ENVIRON["ngs_blk_xs"]!=""))
        print "[" stage "] NOTE: no screening runlevel in template; screening controls were not inserted." > "/dev/stderr"
}
function exciton_eigenvectors(    save_wf) {
    save_wf=ENVIRON["BSE_WRITE_WF"]
    # This flag is essential for YamboPy BS_EIGENSTATES and k-resolved exciton maps.
    if (!has("WRbsWF")) fail("WRbsWF flag absent; regenerate with -V all")
    if (!has("BSSmod")) fail("BSSmod absent; require diagonalization solver")
    scalar("BSSmod",quote("d"),"required")
    flag("WRbsWF",save_wf=="1","required")
}

function qp(prefix,    mode,shift,db) {
    mode=tolower(ENVIRON["QP_MODE"]); shift=ENVIRON["SCISSOR_EV"]
    if(!number(shift)) fail("SCISSOR_EV must be numeric")
    if(mode=="none") { if(shift+0!=0) fail("QP_MODE=none conflicts with a nonzero SCISSOR_EV"); db="none" }
    else if(mode=="scissor") db="none"
    else if(mode=="gwdb") { if(shift+0!=0) fail("GW database + nonzero scissor would double-count corrections"); db="E < " ENVIRON["QP_DB"] }
    else fail("QP_MODE must be none, scissor or gwdb")
    scalar(prefix "db",quote(db),"add")
    block(prefix "_E",shift " | 1.000000 | 1.000000 |","add")
}
function parallel(prefix,envprefix,default_roles,policy,    layout,roles,a,b,n,m,j,prod) {
    layout=trim(ENVIRON[envprefix "_CPU_LAYOUT"]); roles=trim(ENVIRON[envprefix "_ROLES"])
    if(layout=="" || tolower(layout)=="auto") {
        scalar(prefix "_CPU",quote(""),"optional"); scalar(prefix "_ROLEs",quote(""),"optional"); return
    }
    gsub(/[.,]+/," ",layout); layout=trim(layout); n=split(layout,a,/[ \t]+/); prod=1
    if(roles=="") roles=default_roles
    gsub(/[.,]+/," ",roles); roles=trim(roles); m=split(roles,b,/[ \t]+/)
    if(m!=n) fail(envprefix " CPU/role counts differ")
    for(j=1;j<=n;j++) { if(!positive_int(a[j])) fail("Invalid MPI factor in " envprefix); prod*=a[j] }
    if(prod!=ENVIRON["MPI_RANKS"]+0) fail(envprefix " MPI product=" prod ", MPI_RANKS=" ENVIRON["MPI_RANKS"])
    scalar(prefix "_CPU",quote(layout),policy); scalar(prefix "_ROLEs",quote(roles),policy)
}
function common_parallel(    j,n,k,a) {
    if(!positive_int(ENVIRON["MPI_RANKS"])) fail("MPI_RANKS must be a positive integer")
    parallel("DIP","DIP","k c v","add")
    if(stage=="patch_IPA" || stage=="patch_GW" || stage=="patch_BSE" || stage=="patch_COL") {
        parallel("X_q_0","X_Q0","k c v","add")
        parallel("X_all_q","X_ALL","g q k c v","add")
        # Legacy template variable is handled only when actually present.
        parallel("X","X","q k c v","required")
    }
    if(stage=="patch_GW") parallel("SE","SE","q qp b","add")
    if(stage=="patch_BSE") parallel("BS","BS","k eh t","add")
    if(stage=="patch_COL") parallel("COLL","COLL","k q b","required")
    n=split("DIP X SE BS NL COLL",a," ")
    for(j=1;j<=n;j++) {
        k=a[j] "_Threads"
        if(ENVIRON[toupper(a[j]) "_THREADS"]!="" && ENVIRON[toupper(a[j]) "_THREADS"] !~ /^[0-9]+$/) fail("Invalid thread count for " k)
        scalar(k,(ENVIRON[toupper(a[j]) "_THREADS"]!="" ? ENVIRON[toupper(a[j]) "_THREADS"] : "0"),"optional")
    }
    if(stage=="patch_BSE" && ENVIRON["BS_DIAGO_CPUS"]!="") {
        if(has("BS_nCPU_LinAlg_DIAGO")) k="BS_nCPU_LinAlg_DIAGO"; else k="BS_nCPU_diago"
        optional_env(k,"BS_DIAGO_CPUS","","required")
    }
    if(stage=="patch_BSE" && ENVIRON["BS_INV_CPUS"]!="") {
        if(has("BS_nCPU_LinAlg_INV")) k="BS_nCPU_LinAlg_INV"; else k="BS_nCPU_invert"
        optional_env(k,"BS_INV_CPUS","","required")
    }
}
{ sub(/\r$/, ""); lines[NR]=$0 }
END {
    stage=ENVIRON["PATCH_STAGE"]
    # Index complete blocks; never assume a block has only one data row.
    for(i=1;i<=NR;i++) {
        raw=lines[i]; disabled=(raw ~ /^[ \t]*#/); c=disabled?uncomment(raw):raw; c=code(c)
        k=""; typ=""; last=i
        if(c ~ /^%[ \t]*[A-Za-z_][A-Za-z0-9_]*[ \t]*$/) {
            k=c; sub(/^%[ \t]*/,"",k); k=trim(k); typ="block"
            for(last=i+1;last<=NR;last++) {
                close_line=disabled?uncomment(lines[last]):lines[last]
                if(code(close_line)=="%") break
            }
            if(last>NR) fail("Unterminated block % " k)
        } else if(c ~ /^[A-Za-z_][A-Za-z0-9_]*[ \t]*=/) {
            k=c; sub(/[ \t]*=.*/,"",k); typ="scalar"
        } else if(c ~ /^[A-Za-z_][A-Za-z0-9_]*$/) { k=c; typ="flag" }
        if(k!="") {
            kind[i]=typ; finish[i]=last; disabled_node[i]=disabled
            if(disabled) { if(!commented_index[k]) commented_index[k]=i }
            else { active_count[k]++; active_index[k]=i }
        }
        i=last
    }

    cutoff(); screening(); range("BSEBands",ENVIRON["i_band"],ENVIRON["f_band"])
    range("BSEQptR",ENVIRON["iq"],ENVIRON["fq"])
    scalar("BSKmod",quote(ENVIRON["BSKmod"]),"required"); scalar("BSEmod",quote(ENVIRON["BSEmod"]),"optional")
    optional_env("BSENGexx","BSE_Ex","mHa","optional")
    optional_env("BSENGBlk","BSE_Blk",ENVIRON["BSE_Blk_unit"],"optional")
    energy_range("BEnRange",ENVIRON["iE"],ENVIRON["fE"],"eV"); optional_env("BEnSteps","stepE","","optional")
    if(ENVIRON["BSE_DAMP"]!="") {
        if(!number(ENVIRON["BSE_DAMP"])||ENVIRON["BSE_DAMP"]+0<0) fail("BSE_DAMP must be nonnegative")
        block("BDmRange",ENVIRON["BSE_DAMP"] " | " ENVIRON["BSE_DAMP"] " | eV","required")
    }
    field("BLongDir","required"); qp("KfnQP"); exciton_eigenvectors(); common_parallel()

    if(bad) exit 2
    for(i=1;i<=NR;i++) {
        if(skip[i]) continue
        if(i in replacement) print replacement[i]; else print lines[i]
    }
    if(nappend) print "\n# Explicit controls added by " stage
    for(i=1;i<=nappend;i++) print appended[i]
    print "[" stage "] " changes " input entries set; MPI_RANKS=" ENVIRON["MPI_RANKS"] > "/dev/stderr"
}

' "$infile" > "$tmp"
mv -- "$tmp" "$outfile"
trap - EXIT
printf '[%s] Wrote %s\n' "$PATCH_STAGE" "$outfile"
