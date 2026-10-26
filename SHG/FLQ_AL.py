#!/usr/bin/env python3
"""Shared readers, Berry/mesh utilities and audit CLI for the four-file FLQ suite.

Derived from the supplied, validated Lumen/YamboPy postprocessors. No legacy
script is imported. The supported database model is a spin-degenerate,
nonmagnetic insulator, with one SIN/SOFTSIN driving field (hBN / AlAs).
"""
from __future__ import annotations
import argparse
import contextlib
import hashlib
import itertools
import json
import logging
import os
from pathlib import Path
import re
import sys
import tempfile
import traceback
import numpy as np

HA_TO_EV = 27.211386245988
AU_TIME_TO_FS = 0.024188843265857
BOHR_TO_ANG = 0.529177210903
SCHEMA = "FLQ-four-scripts-1.0"

class CompareError(RuntimeError):
    pass

def jsafe(x):
    if isinstance(x, np.ndarray): return jsafe(x.tolist())
    if isinstance(x, (complex, np.complexfloating)):
        return {"real": float(x.real), "imag": float(x.imag)}
    if isinstance(x, (bool, np.bool_)): return bool(x)
    if isinstance(x, (float, np.floating)): return float(x) if np.isfinite(x) else None
    if isinstance(x, (int, np.integer)): return int(x)
    if isinstance(x, Path): return str(x)
    if isinstance(x, dict): return {str(k): jsafe(v) for k,v in x.items()}
    if isinstance(x, (tuple,list)): return [jsafe(v) for v in x]
    return x

def relnorm(delta, ref):
    den=float(np.linalg.norm(ref)); num=float(np.linalg.norm(delta))
    return num/den if den else (0.0 if not num else float("inf"))

def numeric(var):
    return np.asarray(np.ma.asarray(var[:], dtype=float).filled(np.nan))

def need_finite(name, a):
    if not np.all(np.isfinite(a)): raise CompareError(f"{name} contains NaN/inf")

class Tee:
    def __init__(self, *streams): self.streams=streams
    def write(self, text):
        for stream in self.streams: stream.write(text); stream.flush()
        return len(text)
    def flush(self):
        for stream in self.streams: stream.flush()
    def isatty(self): return False

@contextlib.contextmanager
def run_log(path):
    path=Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        with contextlib.redirect_stdout(Tee(sys.stdout, fh)), contextlib.redirect_stderr(Tee(sys.stderr, fh)):
            print(SCHEMA)
            try:
                yield
                print("STATUS: COMPLETE")
            except Exception:
                print("STATUS: FAILED; no successful result is implied")
                traceback.print_exc()
                raise

def atomic_npz(path, **data):
    path=Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    for key,value in data.items():
        if np.asarray(value).dtype.hasobject:
            raise TypeError(f"Object/pickle array forbidden: {key}")
    fd,tmp=tempfile.mkstemp(prefix='.'+path.name, dir=path.parent)
    try:
        with os.fdopen(fd,'wb') as fh: np.savez_compressed(fh, **data)
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)

def atomic_text(path, text):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix='.'+path.name, dir=path.parent)
    try:
        with os.fdopen(fd,'w') as fh: fh.write(text)
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)

def load_npz(path):
    with np.load(path, allow_pickle=False) as z: return {key:z[key] for key in z.files}

def source_stamp(paths):
    out={}
    for p in map(Path,paths):
        st=p.stat(); h=hashlib.sha256()
        with p.open('rb') as fh:
            for block in iter(lambda: fh.read(1024*1024), b''): h.update(block)
        out[str(p.resolve())]={'size':st.st_size,'mtime_ns':st.st_mtime_ns,'sha256':h.hexdigest()}
    return out

def check_source_stamp(stamps):
    # Hash only if metadata changed. Unchanged input bytes remain valid after a copy.
    for name,record in stamps.items():
        p=Path(name)
        if not p.is_file(): raise CompareError(f"Missing source database: {p}")
        st=p.stat()
        if st.st_size!=record['size'] or st.st_mtime_ns!=record['mtime_ns']:
            if source_stamp([p])[str(p.resolve())]['sha256']!=record['sha256']:
                raise CompareError(f"Stale extraction: {p} changed; rerun FLQuasi.py --force")

def calculation_dirs(work, pattern):
    out=sorted(p for p in Path(work).glob(pattern) if p.is_dir())
    if pattern=='FLQ_[0-9][0-9][0-9]':
        out=[p for p in out if re.fullmatch(r'FLQ_[0-9]{3}',p.name)]
    if not out: raise CompareError(f"No calculation directories match {pattern!r}")
    return out

def load_pair(calc, out_name='PyFLQ', verify_sources=False):
    root=Path(calc)/out_name
    q=load_npz(root/'quasienergies.npz'); f=load_npz(root/'fourier.npz')
    if str(q.get('schema',''))!=SCHEMA or str(f.get('schema',''))!=SCHEMA:
        raise CompareError('Unsupported extraction schema; rerun FLQuasi.py')
    if str(q['generation'])!=str(f['generation']):
        raise CompareError('Quasienergy and Fourier files belong to different generations')
    meta=json.loads(str(q['metadata_json']))
    if not meta.get('passed'): raise CompareError('Extraction has not passed its checks')
    if verify_sources: check_source_stamp(meta['sources'])
    return q,f,meta

def run_settings(calc):
    # Reports document effective settings; input text is retained for provenance.
    paths=sorted(Path(calc).glob('r-*nloptics*'),key=lambda p:p.stat().st_mtime_ns)
    text=paths[-1].read_text(errors='replace') if paths else ''
    def last(pattern):
        a=re.findall(pattern,text,re.M|re.I); return a[-1] if a else None
    def scalar(key): return last(r'^[ \t]*(?:\|[ \t]*)?'+re.escape(key)+r'[ \t]*=[ \t]*([^#\n]+)')
    clean=re.sub(r'(?m)^[ \t]*\|[ \t]?', '', text)
    def numeric_block(key):
        hits=re.findall(r'(?ms)^[ \t]*%[ \t]*'+re.escape(key)+r'\b([^%]*?)^[ \t]*%',clean)
        if not hits:return None
        body=' '.join(line.split('#')[0] for line in hits[-1].splitlines())
        return [float(x.replace('D','E').replace('d','e')) for x in re.findall(r'[-+]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eEdD][-+]?[0-9]+)?',body)]
    return {'report':str(paths[-1]) if paths else None,
            'correlation':scalar('NLCorrelation'), 'integrator':scalar('NLintegrator'),
            'qp_database':scalar('GfnQPdb'), 'qp_parameters':numeric_block('GfnQP_E'),
            'band_window':numeric_block('NLBands'), 'damping':scalar('NLDamping'),
            'symmetrize':last(r'Symmetrize current and pol\s*:\s*(yes|no)'),
            'restart_complete':bool(re.search(r'Number of steps done\s*:\s*100\.',text)),
            'game_over':'Game Over' in text}

def divisors_desc(n: int):
    out = []
    for d in range(2, n + 1):
        if n % d == 0:
            out.append(d)
    return sorted(out, reverse=True)

def periodic_tree(kred: np.ndarray):
    from scipy.spatial import cKDTree

    k = np.mod(np.asarray(kred, float), 1.0)
    k[np.isclose(k, 1.0, atol=1e-12)] = 0.0
    return k, cKDTree(k, boxsize=1.0)

def permutation_cycles(perm: np.ndarray):
    n = len(perm)
    seen = np.zeros(n, bool)
    cycles = []

    for i in range(n):
        if seen[i]:
            continue
        cyc = []
        j = i
        while not seen[j]:
            seen[j] = True
            cyc.append(j)
            j = int(perm[j])
        if j != i:
            raise CompareError(
                "Neighbor map contains a non-cycle collision."
            )
        cycles.append(cyc)
    return cycles

def infer_axis_neighbor(kred: np.ndarray, axis: int,
                        tol: float):
    k, tree = periodic_tree(kred)
    nk = len(k)

    for N in divisors_desc(nk):
        delta = np.zeros(3)
        delta[axis] = 1.0 / N
        target = np.mod(k + delta[None, :], 1.0)

        dist, idx = tree.query(target, k=1)
        idx = np.asarray(idx, int)

        if float(np.max(dist)) > tol:
            continue
        if len(np.unique(idx)) != nk:
            continue

        cycles = permutation_cycles(idx)
        lengths = {len(c) for c in cycles}
        if lengths == {N}:
            return {
                "N": int(N),
                "delta_red": delta,
                "neighbor": idx,
                "cycles": cycles,
                "max_match_distance": float(np.max(dist)),
            }

    return {
        "N": 1,
        "delta_red": np.zeros(3),
        "neighbor": np.arange(nk, dtype=int),
        "cycles": [[i] for i in range(nk)],
        "max_match_distance": 0.0,
    }

def infer_regular_grid(kred: np.ndarray, tol: float):
    axes = {}
    shape = []

    for a in range(3):
        info = infer_axis_neighbor(kred, a, tol)
        axes[a] = info
        shape.append(info["N"])

    varying = [a for a in range(3) if axes[a]["N"] > 1]

    product = int(np.prod(shape))
    if product != len(kred):
        raise CompareError(
            "Could not infer a separable regular reciprocal grid. "
            f"Detected shape={tuple(shape)}, product={product}, "
            f"nk={len(kred)}. "
            "Do not continue with Berry polarization until the topology "
            "is understood."
        )

    return tuple(shape), axes, varying

def _decode_bulk_DIP_S(var, nk: int, nb_expected: int):
    a = np.asarray(var[:], float)
    shape = a.shape
    dims = list(var.dimensions)

    if len(shape) != 6:
        raise CompareError(
            f"DIP_S must have 6 dimensions; got shape={shape}, dims={dims}"
        )

    standard = (
        shape[0] == 1
        and shape[1] == nk
        and shape[2] == 6
        and shape[3] == nb_expected
        and shape[4] == nb_expected
        and shape[5] == 2
    )

    if standard:
        ispin = 0
        ik = 1
        idir = 2
        band_axes = [3, 4]
        ic = 5
        detection = "STANDARD_POSITIONAL_1_NK_6_NB_NB_2"

    else:
        ic_candidates = [i for i, n in enumerate(shape) if n == 2]
        ik_candidates = [i for i, n in enumerate(shape) if n == nk]
        ispin_candidates = [i for i, n in enumerate(shape) if n == 1]

        if len(ic_candidates) != 1:
            raise CompareError(
                f"Cannot uniquely identify DIP_S complex axis: "
                f"shape={shape}, dims={dims}, candidates={ic_candidates}"
            )
        if len(ik_candidates) != 1:
            raise CompareError(
                f"Cannot uniquely identify DIP_S k axis: "
                f"shape={shape}, dims={dims}, nk={nk}, "
                f"candidates={ik_candidates}"
            )
        if len(ispin_candidates) != 1:
            raise CompareError(
                f"Cannot uniquely identify DIP_S n_sp_pol=1 axis: "
                f"shape={shape}, dims={dims}, candidates={ispin_candidates}"
            )

        ic = ic_candidates[0]
        ik = ik_candidates[0]
        ispin = ispin_candidates[0]

        used0 = {ic, ik, ispin}
        remaining = [i for i in range(6) if i not in used0]

        if nb_expected == 6:
            raise CompareError(
                "Non-standard DIP_S layout with nb_expected=6 is ambiguous: "
                f"shape={shape}, dims={dims}. "
                "Expected standard layout (1,nk,6,6,6,2)."
            )

        idir_candidates = [i for i in remaining if shape[i] == 6]
        if len(idir_candidates) != 1:
            raise CompareError(
                f"Cannot uniquely identify DIP_S direction axis: "
                f"shape={shape}, dims={dims}, candidates={idir_candidates}"
            )
        idir = idir_candidates[0]

        band_axes = [
            i for i in remaining
            if i != idir and shape[i] == nb_expected
        ]
        if len(band_axes) != 2:
            raise CompareError(
                f"Expected exactly two DIP_S band axes of size {nb_expected}; "
                f"shape={shape}, dims={dims}, candidates={band_axes}"
            )

        detection = "SIZE_BASED_FALLBACK"

    order = [
        ik,
        idir,
        band_axes[0],
        band_axes[1],
        ic,
        ispin,
    ]

    x = a.transpose(order)
    x = x[..., 0]  # spin channel
    z = x[..., 0] + 1j * x[..., 1]

    expected = (nk, 6, nb_expected, nb_expected)
    if z.shape != expected:
        raise CompareError(
            f"Canonical DIP_S shape is {z.shape}, expected {expected}"
        )

    return z, {
        "raw_shape": list(shape),
        "raw_dimensions": dims,
        "axis_detection": detection,
        "complex_axis": int(ic),
        "direction_axis": int(idir),
        "k_axis": int(ik),
        "band_axes": [int(x) for x in band_axes],
        "spin_axis": int(ispin),
        "canonical_shape": list(z.shape),
    }

def load_yambo_DIP_S(calc: Path, nk: int, nb_expected: int):
    from netCDF4 import Dataset

    path = calc / "ndb.Overlap"
    if not path.is_file():
        raise CompareError(
            f"Missing {path}; exact Berry reconstruction requires ndb.Overlap."
        )

    with Dataset(str(path)) as db:
        if "DIP_S" not in db.variables:
            raise CompareError(
                f"{path} does not contain the expected monolithic DIP_S variable."
            )

        DIP, schema = _decode_bulk_DIP_S(
            db.variables["DIP_S"],
            nk=nk,
            nb_expected=nb_expected,
        )

        bands = None
        if "Bands" in db.variables:
            bands = numeric(db.variables["Bands"]).reshape(-1)

    return DIP, {
        "path": str(path),
        "storage_mode": "MONOLITHIC_BULK_DIP_S",
        "Bands": bands,
        "schema": schema,
        "numpy_band_axis_convention": "USE_DIP_k_dir_TRANSPOSE",
    }

def berry_difference_metrics(P: np.ndarray, Pref: np.ndarray):
    dP = np.asarray(P, float) - np.asarray(P, float)[0]
    dR = np.asarray(Pref, float) - np.asarray(Pref, float)[0]

    err = relnorm(dP - dR, dR)

    den = float(np.vdot(dP.ravel(), dP.ravel()).real)
    if den > 0:
        scalar = float(
            np.vdot(dP.ravel(), dR.ravel()).real / den
        )
        scaled = relnorm(
            scalar * dP - dR,
            dR,
        )
    else:
        scalar = float("nan")
        scaled = float("inf")

    return err, scalar, scaled

def harmonic_fit(times_au: np.ndarray, signal: np.ndarray,
                 omega: float, xorder: int,
                 field_initial_time: float):
    T = np.asarray(times_au, float) - float(field_initial_time)
    M = np.zeros(
        (len(T), 2 * xorder + 1),
        complex,
    )
    M[:, 0] = 1.0

    for n in range(1, xorder + 1):
        M[:, n] = np.exp(-1j * n * omega * T)
        M[:, n + xorder] = np.exp(+1j * n * omega * T)

    raw = np.linalg.lstsq(
        M,
        np.asarray(signal, complex),
        rcond=1e-12,
    )[0]

    return raw[:xorder + 1]

def normalize_one_harmonic(raw: np.ndarray, nldb,
                           ana, order: int):
    from yambopy.nl.external_efield import Divide_by_the_Field

    return (
        raw[order]
        * Divide_by_the_Field(
            nldb.Efield[0], order
        )
        * ana.get_Unit_of_Measure(order)
    )


def load_vbands(calc):
    """Decode the uploaded Lumen layout; preserve raw record indices and masks."""
    from netCDF4 import Dataset
    calc=Path(calc)
    main=calc/'ndb.RT_V_bands'; section=calc/'ndb.RT_V_bands_K_section'
    with Dataset(main) as db, Dataset(section) as sec:
        tv=db['IO_TIME_points']; times=numeric(tv).ravel()
        unit=str(getattr(tv,'units','')).strip().lower()
        if unit in ('fs','femtosecond','femtoseconds'): times/=AU_TIME_TO_FS
        elif unit not in ('','au','a.u.','a.u','atomic units','atomic_units'):
            raise CompareError(f'Unknown checkpoint time unit: {unit}')
        var=sec['V_bands']; dims=list(var.dimensions)
        wanted=['V_bands_IO_time_step','n_sp_pol','QP_nk','nb_full','RT_nbands','complex']
        if len(dims)!=6 or set(dims)!=set(wanted):
            raise CompareError(f'Unsupported V_bands dimensions: {dims}')
        if var.shape[dims.index('n_sp_pol')]!=1 or var.shape[dims.index('complex')]!=2:
            raise CompareError('Only one spin channel with real/imaginary storage is supported')
        select=[slice(None)]*6; select[dims.index('n_sp_pol')]=0
        a=np.ma.asarray(var[tuple(select)],dtype=float).filled(np.nan)
        remaining=[d for d in dims if d!='n_sp_pol']
        target=['V_bands_IO_time_step','QP_nk','RT_nbands','nb_full','complex']
        a=np.asarray(a).transpose([remaining.index(d) for d in target])
        C=a[...,0]+1j*a[...,1]
        bands=numeric(db['RT_bands_kpts']).ravel()[:2].round().astype(int)
        nf=int(abs(round(numeric(db['IO_Floquet_order']).ravel()[0])))
    if nf<1 or len(times)!=len(C) or bands[1]-bands[0]+1!=C.shape[2]:
        raise CompareError('Inconsistent band window, Fourier order or record dimensions')
    with Dataset(calc/'ndb.Nonlinear') as db:
        v=db['Field_Freq_1']; omega=float(numeric(v).ravel()[0])
        unit=str(getattr(v,'units','')).strip().lower()
        if unit=='ev': omega/=HA_TO_EV
        elif unit not in ('','au','a.u.','a.u','ha','hartree','atomic units','atomic_units'):
            raise CompareError(f'Unknown frequency unit: {unit}')
    if not np.isfinite(omega) or omega<=0: raise CompareError('Invalid driving frequency')
    raw=np.arange(len(times)); finite=np.isfinite(times)&np.isfinite(C).all(axis=(1,2,3))
    bad=raw[~finite].tolist()
    if bad: print(f'{calc.name}: nonfinite raw records excluded: {bad}')
    times,C,raw=times[finite],C[finite],raw[finite]
    # Only a redundant record AFTER the required 4T checkpoint may be removed.
    last_required=2*nf+3; drop=[]
    for i in range(1,len(times)):
        if times[i]==times[i-1] and raw[i]>last_required:
            if not np.allclose(C[i],C[i-1],rtol=1e-10,atol=1e-13):
                raise CompareError('Duplicate timestamp carries different wavefunctions')
            drop.append(i)
    if drop:
        print(f'{calc.name}: redundant post-checkpoint records excluded: {raw[drop].tolist()}')
        keep=np.ones(len(times),bool);keep[drop]=False
        times,C,raw=times[keep],C[keep],raw[keep]
    if len(times)<3 or np.any(np.diff(times)<=0):
        raise CompareError('V_bands timestamps are not strictly increasing')
    return times,C,raw,omega,{'band_lo':int(bands[0]),'band_hi':int(bands[1]),
        'nf':nf,'decoded_real_imag_shape':list(a.shape),
        'dropped_nonfinite_records':bad,'dropped_duplicate_count':len(drop)}


def exact_samples(times, raw, period, nf, tol=1e-8):
    """Select exact first-period samples and T/4T by their original record numbers."""
    ns=2*nf+1
    loc={int(r):i for i,r in enumerate(raw)}
    required=list(range(1,ns+3))
    missing=[r for r in required if r not in loc]
    if missing: raise CompareError(f'Missing required Floquet records: {missing}')
    fit=np.array([loc[r] for r in range(1,ns+1)])
    i0=loc[1]; iT=loc[ns+1]; i4=loc[ns+2]
    tau=(times-times[i0])/period
    errors={'phase_grid':float(np.max(abs(tau[fit]-np.arange(ns)/ns))),
            'T':float(abs(tau[iT]-1)), 'fourT':float(abs(tau[i4]-4))}
    if max(errors.values())>tol:
        raise CompareError(f'Exact-checkpoint validation failed: {errors}')
    return i0,fit,iT,i4,errors


def load_lattice(save):
    """Use the *driven* SAVE symmetries, never the original equilibrium group."""
    from yambopy import YamboLatticeDB
    from netCDF4 import Dataset
    path=Path(save)/'ns.db1'
    lat=YamboLatticeDB.from_db_file(filename=str(path),Expand=False)
    if int(lat.spinor_components)!=1 or int(lat.spin)!=1:
        raise CompareError('This release supports nonmagnetic, spin-degenerate insulators only')
    if not np.isclose(float(lat.nelectrons)/2,round(float(lat.nelectrons)/2)):
        raise CompareError('Fractional occupations require an occupation-matrix implementation')
    full=YamboLatticeDB.from_db_file(filename=str(path),Expand=True)
    mapping=np.asarray(full.BZ_to_IBZ_indexes,int)
    if mapping.shape!=(len(full.red_kpoints),) or mapping.min()<0 or mapping.max()>=len(lat.red_kpoints):
        raise CompareError('Invalid YamboPy full-BZ to IBZ index map')
    shape,_,varying=infer_regular_grid(np.asarray(full.red_kpoints,float),2e-5)
    with Dataset(path) as db:
        v=db['EIGENVALUES']; en=numeric(v)
        unit=str(getattr(v,'units','')).strip().lower()
    if en.ndim==3 and en.shape[0]==1: en=en[0]
    if en.ndim!=2 or en.shape[0]!=len(lat.red_kpoints):
        raise CompareError(f'Unsupported KS eigenvalue shape {en.shape}')
    if unit not in ('ev','ha','hartree','au','a.u.',''):
        raise CompareError(f'Unknown KS energy unit: {unit}')
    if unit!='ev': en*=HA_TO_EV
    nocc=int(lat.nbandsv)
    if np.min(en[:,nocc])-np.max(en[:,nocc-1])<=0:
        print('WARNING: nonpositive KS fundamental gap; verify insulating integer occupations')
    time_reversal=np.asarray(lat.time_rev_list,bool)
    data={'kpoints_red':np.asarray(lat.red_kpoints,float),
          'kpoints_full_red':np.asarray(full.red_kpoints,float),
          'bz_to_ibz':mapping,'symmetry_indexes':np.asarray(full.symmetry_indexes,int),
          'symmetry_cart':np.asarray(lat.sym_car,float),'symmetry_antiunitary':time_reversal,
          'lattice_bohr':np.asarray(lat.lat,float),'reciprocal_no2pi_bohr':np.asarray(lat.rlat,float),
          'volume_bohr3':float(abs(lat.lat_vol)),'spin_factor':float(lat.nspin),
          'nocc':nocc,'ks_energies_eV':en,'ks_vbm_eV':float(np.max(en[:,nocc-1])),
          'grid_shape':np.asarray(shape,int),'varying_axes':np.asarray(varying,int)}
    return data


def full_occupied(C, nocc, band_lo, band_hi):
    """Frozen lower occupied orbitals remain in the full Berry determinant."""
    C=np.asarray(C,complex)
    if C.shape[-2:]!=(band_hi-band_lo+1,nocc-band_lo+1):
        raise CompareError(f'Dynamic coefficient dimensions incompatible with bands {band_lo}:{band_hi}')
    V=np.zeros(C.shape[:-2]+(band_hi,nocc),complex)
    for b in range(band_lo-1): V[...,b,b]=1
    V[...,band_lo-1:band_hi,band_lo-1:nocc]=C
    return V


def berry_polarization(C, DIP, axes, lattice, volume, spin, nocc, band_lo, band_hi):
    """Validated positive-sign / transposed-DIP convention, vectorized over k.

    C is [time, full_BZ_k, KS_basis, propagated_occupied]. Return cartesian
    polarization in atomic units. Unwrapping follows time on every k string.
    """
    C=np.asarray(C,complex); nt,nk=C.shape[:2]
    if DIP.shape[:2]!=(nk,6): raise CompareError('Overlap k/direction dimensions disagree')
    phases={}; p_red=np.zeros((nt,3))
    for axis,info in axes.items():
        if info['N']<=1: continue
        phases[axis]=np.empty((nt,len(info['cycles'])))
    for it in range(nt):
        V=full_occupied(C[it],nocc,band_lo,band_hi)
        for axis,phase in phases.items():
            info=axes[axis]; vp=V[info['neighbor']]
            link=np.einsum('kbi,kbc,kcj->kij',V.conj(),DIP[:,axis].transpose(0,2,1),vp,optimize=True)
            sign,logabs=np.linalg.slogdet(link)
            if np.any(~np.isfinite(logabs)) or np.any(logabs<np.log(1e-14)):
                raise CompareError(f'Near-singular occupied Berry link at time {it}, axis {axis}')
            cycle=np.asarray(info['cycles'],int)
            phase[it]=np.angle(np.prod(sign[cycle],axis=1))
    for axis,phase in phases.items():
        p_red[:,axis]=spin/(2*np.pi)*np.unwrap(phase,axis=0).mean(axis=1)
    return p_red@np.asarray(lattice,float)/volume


def polar_projector(q, field, setting='none'):
    """Unitary field-preserving polar-vector projector, not a fit to Pstored."""
    e=np.asarray(field,float); e=e/np.linalg.norm(e)
    if setting=='none': return np.eye(3)
    syms=np.asarray(q['symmetry_cart'],float)
    anti=np.asarray(q['symmetry_antiunitary'],bool)
    keep=[]
    for s,isanti in zip(syms,anti):
        if not isanti and np.allclose(s@e,e,atol=2e-6,rtol=0):
            if not np.allclose(s@s.T,np.eye(3),atol=2e-6):
                raise CompareError('Stored Cartesian symmetry is not orthogonal')
            keep.append(s)
    if not keep: raise CompareError('No unitary field-preserving symmetries found')
    proj=np.mean(keep,axis=0)
    if not np.allclose(proj@proj,proj,atol=2e-6):
        raise CompareError('Residual symmetry average is not a projector')
    if setting=='xy':
        target=np.array([[.5,.5,0],[.5,.5,0],[0,0,1]])
        if not np.allclose(abs(e),[2**-.5,2**-.5,0],atol=1e-6) or e[0]*e[1]<=0:
            raise CompareError('xy projection requires E parallel to [110]')
        mirror=np.array([[0,1,0],[1,0,0],[0,0,1]])
        if not any(np.allclose(s,mirror,atol=2e-6) for s in keep):
            raise CompareError('The x/y mirror is not a symmetry of this driven SAVE')
        return target
    return proj


def project_delta(P, projector):
    P=np.asarray(P,float)
    return P[0]+(P-P[0])@projector.T


def load_polarization_at(db, sample_times):
    """Reject sparse restart histories; interpolate only a genuinely dense record."""
    from scipy.interpolate import CubicSpline
    t=np.asarray(db.IO_TIME_points,float).ravel()
    if len(db.Polarization)!=1: raise CompareError('Expected one Floquet frequency per folder')
    P=np.asarray(db.Polarization[0],float).T
    if P.shape!=(len(t),3): raise CompareError(f'Unexpected polarization shape {P.shape}')
    need_finite('stored times',t);need_finite('stored polarization',P)
    if np.any(np.diff(t)<=0): raise CompareError('Stored polarization timestamps are not increasing')
    norms=np.linalg.norm(P-P[0],axis=1); mx=float(norms.max())
    populated=int(np.count_nonzero(norms>max(mx*1e-12,1e-300)))
    if mx==0 or (len(t)>20 and populated<=2):
        raise CompareError(f'Invalid/sparse stored polarization ({populated}/{len(t)} populated records); '
                           'check completed-run restart contamination, do not relax Berry validation')
    targets=np.asarray(sample_times,float)
    tol=max(float(np.ptp(t))*1e-12,1e-9)
    if targets.min()<t.min()-tol or targets.max()>t.max()+tol:
        raise CompareError('V_bands samples fall outside stored polarization times')
    targets=np.clip(targets,t.min(),t.max())
    out=CubicSpline(t,P,axis=0,extrapolate=False)(targets)
    need_finite('interpolated polarization',out)
    return out


def first_bz(reciprocal, axes):
    """Wigner-Seitz cell of reciprocal lattice; handles hexagonal and skew cells."""
    from scipy.spatial import Voronoi, ConvexHull
    axes=np.asarray(axes,int); dim=len(axes)
    B=np.asarray(reciprocal,float)[axes]*2*np.pi/BOHR_TO_ANG
    if dim==2:
        u=B[0]/np.linalg.norm(B[0]); n=np.cross(B[0],B[1]); n/=np.linalg.norm(n)
        frame=np.stack([u,np.cross(n,u)],axis=1)
    elif dim==3: frame=np.eye(3)
    else: raise CompareError('Only 2D/3D BZ plotting is implemented')
    b=B@frame; previous=None
    for span in (1,2,3,4):
        shifts=np.array(list(itertools.product(range(-span,span+1),repeat=dim)))
        points=shifts@b; origin=np.flatnonzero(np.all(shifts==0,axis=1))[0]
        vor=Voronoi(points); reg=vor.regions[vor.point_region[origin]]
        if not reg or -1 in reg: continue
        vertices=vor.vertices[reg]; hull=ConvexHull(vertices)
        measure=hull.volume
        if previous is not None and np.isclose(measure,previous,rtol=1e-10): break
        previous=measure
    expected=abs(np.linalg.det(b))
    if not np.isclose(measure,expected,rtol=1e-7):
        raise CompareError('Reciprocal Voronoi cell has wrong area/volume')
    if dim==2: vertices=vertices[hull.vertices]; hull=ConvexHull(vertices)
    return vertices,hull,frame


def periodic_images(kred, values, reciprocal, axes, span=1):
    shifts=np.zeros(( (2*span+1)**len(axes),3))
    shifts[:,axes]=np.array(list(itertools.product(range(-span,span+1),repeat=len(axes))))
    k=np.asarray(kred,float); k=k-np.floor(k+.5)
    points=(k[None,:,:]+shifts[:,None,:]).reshape(-1,3)
    return points@reciprocal*2*np.pi/BOHR_TO_ANG, np.tile(values,(len(shifts),)+(1,)*(np.ndim(values)-1))


def read_spectrum(path):
    a=np.atleast_2d(np.loadtxt(path,comments='#'))
    if a.shape[1]<7: raise CompareError(f'{path}: expected E, Im/Re x, Im/Re y, Im/Re z')
    need_finite('spectrum',a)
    order=np.argsort(a[:,0]);a=a[order]
    if len(a)>1 and np.any(np.diff(a[:,0])<=0): raise CompareError('Duplicate spectrum energies')
    return a[:,0],a[:,[2,4,6]]+1j*a[:,[1,3,5]]


def spectrum_text(E,X,label):
    import io
    stream=io.StringIO();a=np.column_stack([E,X[:,0].imag,X[:,0].real,X[:,1].imag,X[:,1].real,X[:,2].imag,X[:,2].real])
    np.savetxt(stream,a,header='E_eV Im_x Re_x Im_y Re_y Im_z Re_z\n'+label)
    return stream.getvalue()


def comparison_metrics(Ea,A,Eb,B,component=1,tol=1e-4):
    if len(Ea)!=len(Eb) or not np.allclose(Ea,Eb,atol=tol,rtol=0):
        raise CompareError('RT/Floquet frequency grids differ; no silent interpolation is performed')
    a=A[:,component];b=B[:,component];scale=max(float(np.max(abs(a))),1e-300)
    meaningful=abs(a)>.05*scale
    rel=np.divide(abs(abs(b)-abs(a)),abs(a),out=np.full(a.shape,np.nan),where=abs(a)>scale*1e-10)
    return {'complex_L2_relative':relnorm(b-a,a),
        'peak_normalized_max_abs_magnitude_error':float(np.max(abs(abs(b)-abs(a)))/scale),
        'median_magnitude_relative_error_above_5pct_peak':float(np.nanmedian(rel[meaningful])) if meaningful.any() else None,
        'npoints':len(Ea),'npoints_above_5pct_peak':int(meaningful.sum()),
        'max_energy_difference_eV':float(np.max(abs(Ea-Eb)))}



# ============================================================================
# Degenerate harmonic susceptibility tensors (Cartesian indices, 1-based CLI)
# ============================================================================

def parse_tensor_component(spec, order):
    """'y' means contracted P_y/E^order; '[1,1,2]' means chi_xxy."""
    s=str(spec).strip().lower()
    if s in ('x','y','z'):
        return {'kind':'effective','output_axis':'xyz'.index(s),'label':s,'raw':s}
    m=re.fullmatch(r'\[\s*([123](?:\s*,\s*[123])*)\s*\]',s)
    if not m:
        raise CompareError('Use --component x/y/z for a contracted response, or '
                           f'--component "[i,j,...]" with {order+1} indices 1,2,3')
    inds=tuple(int(v.strip())-1 for v in m.group(1).split(','))
    if len(inds)!=order+1:
        raise CompareError(f'Harmonic order {order} requires {order+1} tensor indices '
                           f'(one output and {order} incident fields), not {len(inds)}')
    key='['+','.join(str(i+1) for i in (inds[0],)+tuple(sorted(inds[1:])))+']'
    return {'kind':'tensor','index':inds,'key':key,'label':'chi_'+''.join('xyz'[i] for i in inds),
            'output_axis':inds[0],'raw':str(spec)}


def tensor_catalog(order, efield, rotations=None, symmetry_source='not_verified', antiunitary=None):
    """Identify tensor elements obtainable from one degenerate-frequency field.

    A homogeneous tensor is symmetric in its incident-field indices. Crystal
    rotations constrain its invariant subspace. A component is 'extractable'
    only if it lies in the row space of the three measured contractions.
    No material names or assumed crystallographic axes are used.
    """
    from math import factorial
    if order<1:raise CompareError('Harmonic order must be positive')
    if order>6:
        e=np.asarray(efield,float)
        if e.shape!=(3,) or not np.isfinite(e).all() or np.linalg.norm(e)==0:
            raise CompareError('Invalid applied field')
        return {'harmonic_order':int(order),'field_unit_xyz':(e/np.linalg.norm(e)).tolist(),
                'symmetry_source':'individual Cartesian tensor elements not audited above order 6',
                'symmetry_operations':None,'independent_tensor_parameters':None,
                'components':{},'measured_effective':['chi_'+c+'_eff' for c in 'xyz'],
                'convention':'field-contracted susceptibility only for harmonic order >6'}
    field=np.asarray(efield,float).reshape(-1)
    if field.shape!=(3,) or not np.isfinite(field).all() or np.linalg.norm(field)==0:
        raise CompareError('The applied Cartesian electric-field direction is invalid')
    field=field/np.linalg.norm(field)
    field[np.abs(field)<1e-9]=0
    field/=np.linalg.norm(field)
    if rotations is None:rotations=np.eye(3)[None,:,:]
    rotations=np.asarray(rotations,float)
    if rotations.ndim!=3 or rotations.shape[1:]!=(3,3) or not np.isfinite(rotations).all():
        raise CompareError('Point-group rotations must have shape (nsym,3,3)')
    # Time reversal reverses k but is NOT the inversion of the polar response.
    # Treating antiunitary -R maps as unitary spatial symmetries would falsely
    # set chi^(2) to zero even in noncentrosymmetric monolayers.
    excluded_antiunitary=0
    if antiunitary is not None:
        antiunitary=np.asarray(antiunitary,bool).reshape(-1)
        if len(antiunitary)!=len(rotations):
            raise CompareError('Unitary/antiunitary symmetry labels have different lengths')
        excluded_antiunitary=int(np.count_nonzero(antiunitary))
        rotations=rotations[~antiunitary]
        if not len(rotations):
            raise CompareError('No unitary spatial point-group operations supplied')
    # Avoid degeneracy from repeated time-reversed symmetry operations.
    unique=[]
    for R in rotations:
        if not np.allclose(R@R.T,np.eye(3),rtol=0,atol=3e-5):
            raise CompareError('The symmetry database contains non-Cartesian/nonorthogonal rotations')
        if not any(np.allclose(R,U,rtol=0,atol=5e-7) for U in unique):unique.append(R)
    if not any(np.allclose(R,np.eye(3),rtol=0,atol=5e-7) for R in unique):unique.append(np.eye(3))
    labels=list(itertools.combinations_with_replacement(range(3),order))
    columns=[(i,inds) for i in range(3) for inds in labels]
    d=len(columns); eye=np.eye(d)
    unit_tensors=[]
    shape=(3,)*(order+1)
    for i,inds in columns:
        T=np.zeros(shape)
        for perm in set(itertools.permutations(inds)):
            T[(i,)+perm]=1.0
        unit_tensors.append(T)
    if len(unique)==1:
        invariant=eye
    else:
        # Cartesian tensor transformation: T'_i...=R_ia ... T_a...
        old='abcdefgh'[:order+1];new='uvwxyzst'[:order+1]
        rule=','.join(n+o for n,o in zip(new,old))+','+old+'->'+new
        restraints=[]
        for R in unique:
            transformed=np.empty((d,d))
            for col,T in enumerate(unit_tensors):
                rot=np.einsum(rule,*([R]*(order+1)),T,optimize=True)
                transformed[:,col]=[rot[(i,)+inds] for i,inds in columns]
            restraints.append(transformed-eye)
        A=np.concatenate(restraints,axis=0)
        _,s,Vh=np.linalg.svd(A,full_matrices=False)
        # ns.db1 Cartesian rotations can have ~1e-7 rounding residuals.
        # Count approximate symmetry null modes with a tolerant relative SVD cutoff.
        tol=max(1e-8,1e-6*(float(s.max()) if len(s) else 1.))
        rank=int(np.count_nonzero(s>tol))
        invariant=Vh[rank:].T
    contractions=np.zeros((3,d))
    for c,(i,inds) in enumerate(columns):
        multiplicity=factorial(order)
        for j in range(3):multiplicity//=factorial(inds.count(j))
        contractions[i,c]=multiplicity*np.prod(field[np.array(inds,dtype=int)])
    measured=contractions@invariant
    catalog={}
    for col,(i,inds) in enumerate(columns):
        key='['+','.join(str(a+1) for a in (i,)+inds)+']'
        t=invariant[col,:]
        if np.linalg.norm(t)<1e-8:
            state='symmetry_forbidden';coeff=np.zeros(3);error=0.
        else:
            coeff=np.linalg.lstsq(measured.T,t,rcond=1e-10)[0]
            error=float(np.linalg.norm(measured.T@coeff-t)/max(np.linalg.norm(t),1e-300))
            if error<=1e-7 and np.linalg.norm(coeff)<=1e7:
                state='extractable'
                coeff[np.abs(coeff)<1e-9]=0.
            else:
                state='not_identifiable';coeff=np.zeros(3)
        catalog[key]={'label':'chi_'+''.join('xyz'[a] for a in (i,)+inds),
                      'indices_1based':[a+1 for a in (i,)+inds],
                      'status':state,'from_effective_xyz':coeff.tolist() if state=='extractable' else None,
                      'identifiability_error':error if state=='not_identifiable' else 0.0}
    report={'harmonic_order':int(order),'field_unit_xyz':field.tolist(),
            'symmetry_source':str(symmetry_source),'symmetry_operations':len(unique),
            'excluded_antiunitary_operations':excluded_antiunitary,
            'independent_tensor_parameters':int(invariant.shape[1]),
            'convention':'P_i(nw)/E0^n = sum_{j1..jn} chi_i,j1..jn e_j1...e_jn; '
                         'degenerate input indices symmetric; mixed terms include multiplicity',
            'measured_effective':['chi_'+c+'_eff' for c in 'xyz'],
            'components':catalog}
    return report


def tensor_catalog_for_work(work, save, order, field, override_db=None):
    """Prefer the *equilibrium*, pre-field-symmetry-reduction SAVE/ns.db1.

    A driven FixSymm/SAVE only contains field-preserving operations. It must
    not be used to assert a full crystal tensor symmetry without verification.
    When equilibrium symmetry is missing, use identity-only (safe algebra).
    """
    work=Path(work)
    if override_db:
        path=Path(override_db).expanduser()
        if not path.is_absolute():path=work/path
        if not path.is_file():raise CompareError(f'No equilibrium symmetry database: {path}')
        source_type='explicit_equilibrium_symmetry_database'
    else:
        path=work.parent/'SAVE'/'ns.db1'
        if not path.is_file():
            return tensor_catalog(order,field,symmetry_source='unavailable: identity only; full point group unknown')
        source_type='parent_equilibrium_SAVE'
    from yambopy import YamboLatticeDB
    lat=YamboLatticeDB.from_db_file(filename=str(path),Expand=False)
    local=work/str(save)/'ns.db1'
    if local.is_file():
        driven=YamboLatticeDB.from_db_file(filename=str(local),Expand=False)
        if not np.allclose(np.asarray(lat.lat,float),np.asarray(driven.lat,float),rtol=1e-6,atol=1e-5):
            raise CompareError('Equilibrium symmetry database and driven SAVE have different lattice vectors')
        if int(lat.nelectrons)!=int(driven.nelectrons):
            raise CompareError('Equilibrium and driven SAVE have different electron counts')
    rotations=np.asarray(lat.sym_car,float)
    anti=np.asarray(lat.time_rev_list,bool)
    return tensor_catalog(order,field,rotations=rotations,antiunitary=anti,
                          symmetry_source=f'{source_type}: {path.resolve()}')


def tensor_response(catalog, spec, values):
    """Convert field-contracted x/y/z response to an identifiable component."""
    c=parse_tensor_component(spec,catalog['harmonic_order'])
    X=np.asarray(values,complex)
    if X.ndim!=2 or X.shape[1]!=3:
        raise CompareError('Expected frequency-by-3 field-contracted susceptibility')
    if c['kind']=='effective':return X[:,c['output_axis']],{'label':'chi_'+c['label']+'_eff',
                                                             'status':'field_contracted','coefficients_xyz':None}
    if c['key'] not in catalog['components']:
        raise CompareError('Individual Cartesian tensor analysis is not supported at this harmonic order; '
                           'select x, y or z for the field-contracted response')
    entry=catalog['components'][c['key']]
    if entry['status']=='symmetry_forbidden':
        raise CompareError(f'{c["label"]} is forbidden by the supplied crystal symmetry: '
                           'do not reinterpret finite numerical leakage as an allowed susceptibility')
    if entry['status']!='extractable':
        raise CompareError(f'{c["label"]} cannot be separated from other components using this '
                           f'field {catalog["field_unit_xyz"]} and known symmetries. '
                           'Use additional independent field polarizations.')
    factors=np.asarray(entry['from_effective_xyz'],float)
    return X@factors,{'label':c['label'],'canonical_label':entry['label'],
                      'status':'symmetry_isolated','coefficients_xyz':factors.tolist()}


def describe_tensor_catalog(report, requested='y'):
    """Log every allowed and identifiable tensor element without relabeling contractions."""
    print(f'TENSOR DIAGNOSTICS: harmonic order={report["harmonic_order"]}, '
          f'unit field={report["field_unit_xyz"]}')
    print('Point-group source:',report['symmetry_source'])
    print('Unitary spatial symmetry operations:',report['symmetry_operations'],
          '; excluded antiunitary operations:',report.get('excluded_antiunitary_operations',0),
          '; invariant tensor parameters:',report['independent_tensor_parameters'])
    print('Measured quantities: chi_x_eff, chi_y_eff, chi_z_eff (field contractions)')
    for state,heading in [('extractable','Recoverable tensor components'),
                          ('not_identifiable','Not individually identifiable from this field/symmetry information'),
                          ('symmetry_forbidden','Symmetry-forbidden components')]:
        entries=[(k,v) for k,v in report['components'].items() if v['status']==state]
        print(f'{heading} ({len(entries)}):')
        for key,v in entries:
            if state=='extractable':
                parts=['{:+.6g} chi_{}_eff'.format(a,'xyz'[i]) for i,a in enumerate(v['from_effective_xyz']) if abs(a)>1e-8]
                print(f'  {key:13s} {v["label"]:11s} = '+(' '.join(parts) if parts else '0'))
            else:print(f'  {key:13s} {v["label"]:11s}')
    request=parse_tensor_component(requested,report['harmonic_order'])
    print('Requested:',requested,'(vector contraction)' if request['kind']=='effective' else '(tensor component)')
    return request


def adaptive_frequency_grid(E, X, budget=48, component='y', max_step_eV=0.35,
                            min_step_eV=0.0, peak_prominence=0.05,
                            max_peak_anchors=8, slope_weight=0.08,
                            curvature_weight=0.08, relative_floor=0.03):
    """Pick original RT energies, emphasizing actual complex interpolation error.

    Keeps endpoints; respects a maximum allowed gap; seeds prominent extrema;
    greedily spends the remaining budget where the normalized complex RT
    spectrum is hardest to interpolate. Does not invent off-RT frequencies.
    """
    from scipy.signal import find_peaks

    E=np.asarray(E,float)
    X=np.asarray(X,complex)
    if E.ndim!=1 or X.shape!=(len(E),3):
        raise CompareError('Adaptive source requires E[n], complex response[n,3]')
    if len(E)<3 or not np.all(np.isfinite(E)) or np.any(E<=0) or np.any(np.diff(E)<=0):
        raise CompareError('Adaptive source needs >=3 positive increasing finite energies')
    if not np.all(np.isfinite(X)):
        raise CompareError('Non-finite complex response in adaptive source')
    if not (2<=budget<=len(E)):
        raise CompareError(f'Point budget must be in [2,{len(E)}], got {budget}')
    if not (max_step_eV>0 and min_step_eV>=0 and max_step_eV>=min_step_eV):
        raise CompareError('Invalid minimum/maximum frequency spacing')
    if not (0<=peak_prominence<1 and max_peak_anchors>=0 and slope_weight>=0
            and curvature_weight>=0 and 0<relative_floor<=1):
        raise CompareError('Adaptive weights/prominence/floor must be nonnegative and valid')
    if component not in ('x','y','z','all'):
        raise CompareError(f'Unknown adaptive component: {component!r}')

    Z=X if component=='all' else X[:,['x','y','z'].index(component),None]
    amplitude=np.linalg.norm(Z,axis=1)
    global_scale=float(np.max(amplitude))
    if global_scale<=0:
        raise CompareError('RT spectrum is identically zero in selected component(s)')
    floor=relative_floor*global_scale
    deriv1=np.gradient(Z,E,axis=0,edge_order=2)
    deriv2=np.gradient(deriv1,E,axis=0,edge_order=2)
    selected={0,len(E)-1}
    reason={0:'endpoint',len(E)-1:'endpoint'}
    min_eps=max(1e-10*float(np.ptp(E)),1e-11)

    def legal(k,left,right):
        return (E[k]-E[left]>=min_step_eV-min_eps and
                E[right]-E[k]>=min_step_eV-min_eps)

    def add(k,why):
        if k in selected:return
        i=np.searchsorted(sorted(selected),k)
        sel=sorted(selected)
        left=sel[i-1];right=sel[i]
        if not legal(k,left,right):
            raise CompareError(f'Cannot add {E[k]:.8f} eV ({why}): min-step {min_step_eV:g} eV')
        if len(selected)>=budget:
            raise CompareError(f'Budget {budget} too small for required {why} anchors')
        selected.add(k);reason[k]=why

    # A coverage constraint matters even in weak-signal portions of the spectrum.
    while True:
        intervals=list(zip(sorted(selected)[:-1],sorted(selected)[1:]))
        too_wide=[(E[j]-E[i],i,j) for i,j in intervals
                  if E[j]-E[i]>max_step_eV+min_eps]
        if not too_wide:break
        _,i,j=max(too_wide)
        candidates=[k for k in range(i+1,j) if legal(k,i,j)]
        if not candidates:
            raise CompareError('Source RT grid is too coarse to satisfy --max-step-eV; '
                               'compute additional RT frequencies or relax the maximum step')
        target=.5*(E[i]+E[j])
        k=min(candidates,key=lambda k:(abs(E[k]-target),k))
        add(k,'coverage')

    # Prominent extrema are retained even when their local derivative vanishes.
    peaks,properties=find_peaks(amplitude,prominence=peak_prominence*global_scale)
    peak_candidates=sorted(zip(peaks,properties['prominences']),
                           key=lambda item:(-item[1],E[item[0]]))
    chosen_peaks=[]
    for k,prom in peak_candidates[:max_peak_anchors]:
        if k in selected:
            chosen_peaks.append((int(k),float(prom)))
            continue
        sel=sorted(selected);j=np.searchsorted(sel,k)
        if 0<j<len(sel) and legal(int(k),sel[j-1],sel[j]):
            add(int(k),'peak')
            chosen_peaks.append((int(k),float(prom)))
        else:
            # Record ignored peaks; do not misrepresent the spacing constraint.
            pass

    def interval_candidates(left,right):
        k=np.arange(left+1,right,dtype=int)
        if not len(k):return k,np.array([],float)
        mask=((E[k]-E[left]>=min_step_eV-min_eps)&
              (E[right]-E[k]>=min_step_eV-min_eps))
        k=k[mask]
        if not len(k):return k,np.array([],float)
        h=E[right]-E[left]
        frac=(E[k]-E[left])/h
        interp=Z[left]+frac[:,None]*(Z[right]-Z[left])
        local_scale=np.maximum(floor,np.maximum.reduce((
            amplitude[k],np.full(len(k),amplitude[left]),
            np.full(len(k),amplitude[right]))))
        # Exact source-point error of linear interpolation, plus user-inspired
        # 1st- and 2nd-derivative change indicators (tapered at endpoints).
        interp_error=np.linalg.norm(Z[k]-interp,axis=1)/local_scale
        tent=4*frac*(1-frac)
        slope=(slope_weight*(.5*h)*np.linalg.norm(deriv1[k],axis=1)
               *tent/local_scale)
        curvature=(curvature_weight*(h*h/8)*np.linalg.norm(deriv2[k],axis=1)
                   *tent/local_scale)
        return k,interp_error+slope+curvature

    # Max-score bisection is deterministic and uses only actual RT frequencies.
    while len(selected)<budget:
        sel=sorted(selected)
        best=None
        for left,right in zip(sel[:-1],sel[1:]):
            candidates,score=interval_candidates(left,right)
            if len(candidates):
                q=int(np.argmax(score))
                item=(float(score[q]),-int(candidates[q]),int(candidates[q]))
                if best is None or item>best:best=item
        if best is None:
            raise CompareError('No legal additional frequencies at the requested minimum step')
        selected.add(best[2]);reason[best[2]]='adaptive'

    result=np.asarray(sorted(selected),int)
    prediction=np.empty_like(Z)
    for left,right in zip(result[:-1],result[1:]):
        k=np.arange(left,right+1)
        frac=(E[k]-E[left])/(E[right]-E[left])
        prediction[k]=Z[left]+frac[:,None]*(Z[right]-Z[left])
    local_scale=np.maximum(floor,amplitude)
    err=np.linalg.norm(Z-prediction,axis=1)/local_scale
    gap=np.diff(E[result]);hmin=float(np.min(gap));hmax=float(np.max(gap))
    return {
        'indices':result,'energies_eV':E[result],
        'source_energy_eV':E,'source_magnitude':amplitude,
        'pointwise_relative_interpolation_error':err,
        'max_relative_interpolation_error':float(err.max()),
        'rms_relative_interpolation_error':float(np.sqrt(np.mean(err**2))),
        'min_gap_eV':hmin,'max_gap_eV':hmax,
        'reasons':[reason[int(i)] for i in result],
        'detected_peaks':[{'E_eV':float(E[int(k)]), 'prominence':float(p),
                           'anchored':int(k) in result} for k,p in peak_candidates],
        'component':component,
        'source_count':len(E), 'selected_count':len(result),
        'parameters':{'budget':budget,'max_step_eV':max_step_eV,
          'min_step_eV':min_step_eV,'peak_prominence':peak_prominence,
          'max_peak_anchors':max_peak_anchors,'slope_weight':slope_weight,
          'curvature_weight':curvature_weight,'relative_floor':relative_floor}
    }


def reference_on_floquet_grid(E_rt, X_rt, E_fl, tolerance=1e-4, method='exact'):
    """Align complex RT response to FL energies. Interpolation is opt-in."""
    E_rt=np.asarray(E_rt,float);E_fl=np.asarray(E_fl,float)
    X_rt=np.asarray(X_rt,complex)
    if E_rt.ndim!=1 or E_fl.ndim!=1 or X_rt.shape!=(len(E_rt),3):
        raise CompareError('Invalid RT/Floquet energies or susceptibility arrays')
    if len(E_rt)==0 or len(E_fl)==0 or np.any(np.diff(E_rt)<=0):
        raise CompareError('No RT/FL energies or RT grid is not increasing')
    j=np.searchsorted(E_rt,E_fl)
    hi=np.clip(j,0,len(E_rt)-1);lo=np.clip(j-1,0,len(E_rt)-1)
    index=np.where(np.abs(E_rt[hi]-E_fl)<np.abs(E_rt[lo]-E_fl),hi,lo)
    difference=np.abs(E_rt[index]-E_fl)
    exact=difference<=tolerance
    if method=='exact':
        if not np.all(exact):
            k=int(np.argmax(difference))
            raise CompareError(f'Floquet energy {E_fl[k]:.8f} eV has no RT point '
                               f'within {tolerance:g} eV (closest {E_rt[index[k]]:.8f}). '
                               'Use --rt-reference pchip only with an explicitly off-grid frequency plan')
        if len(set(map(int,index)))!=len(index):
            raise CompareError('Multiple Floquet energies map to the same RT frequency')
        out=X_rt[index]
    elif method=='pchip':
        from scipy.interpolate import PchipInterpolator
        if len(E_rt)<3 or np.any(E_fl<E_rt[0]) or np.any(E_fl>E_rt[-1]):
            raise CompareError('Cannot interpolate RT reference outside sampled interval')
        out=PchipInterpolator(E_rt,X_rt.real,axis=0)(E_fl)+1j*PchipInterpolator(
            E_rt,X_rt.imag,axis=0)(E_fl)
        out[exact]=X_rt[index[exact]]
    else:
        raise CompareError(f'Unknown reference method {method!r}')
    return out,index,float(np.max(difference))


def self_test():
    """Synthetic algebra/geometry tests. These are not an HPC/physical convergence test."""
    from types import SimpleNamespace
    from FLQuasi import diagonalize, fourier_modes, reconstruct_states
    from FLQ_HG import berry_validation
    from FLQPlot import periodic_interpolate
    tests=[]
    def passed(name):tests.append(name);print('PASS:',name)
    omega=.21;period=2*np.pi/omega;nf=4;ns=2*nf+1;t0=100.
    tau=np.r_[np.arange(ns)*period/ns,period,4*period]
    xi=np.array([-.37*omega,.2*omega]);nk=4;mix=np.array([[.8,-.6],[.6,.8]],complex)
    def truth(t):
        t=np.atleast_1d(t);D=np.zeros((len(t),nk,4,2),complex)
        for k in range(nk):
            a=.12+.01*k;b=.17+.01*k
            D[:,k,0,0]=np.cos(a);D[:,k,2,0]=np.sin(a)*np.exp(1j*omega*t)
            D[:,k,1,1]=np.cos(b);D[:,k,3,1]=np.sin(b)*np.exp(-2j*omega*t)
        psi=D*np.exp(-1j*t[:,None,None]*xi[None,None,:])[:,:,None,:]
        return np.einsum('tkba,aj->tkbj',psi,mix)
    C=truth(tau);q=diagonalize(C,0,ns,ns+1,period)
    assert q['closure_T'].max()<1e-12 and q['semigroup_4T'].max()<1e-12
    passed('mixed two-state direct monodromy and independent 4T consistency')
    f=fourier_modes(C[:ns],tau[:ns],q['Gamma'],q['epsilon_principal_au'],omega,nf)
    f.update(omega_au=omega,period_au=period,t0_au=t0)
    dense=np.linspace(0,4*period,117)
    assert relnorm(reconstruct_states(q,f,dense)-truth(dense),truth(dense))<1e-12
    passed('compact Fourier gauge plus physical-state mixing at independent off-grid times')
    raw=np.arange(1,len(tau)+1);i0,fit,iT,i4,errs=exact_samples(t0+tau,raw,period,nf)
    assert (i0,iT,i4)==(0,ns,ns+1) and max(errs.values())<1e-12
    passed('exact phase/T/4T record selection')
    V=full_occupied(C[:1],4,3,6)
    assert V.shape==(1,nk,6,4) and np.all(V[...,0,0]==1) and np.all(V[...,1,1]==1)
    assert np.allclose(V[...,2:6,2:4],C[:1])
    passed('NLBands=3..6 includes frozen occupied bands 1 and 2 in full determinant')
    for shape in ((4,3,1),(3,3,2)):
        k=np.array(list(itertools.product(*[np.arange(n)/n for n in shape])))
        shape2,axes,vary=infer_regular_grid(k,2e-8)
        assert tuple(shape2)==shape
        lat=np.eye(3)*7.;N=len(k)
        dyn=np.broadcast_to(C[:,0:1],(len(tau),N,4,2)).copy()
        dip=np.broadcast_to(np.eye(6),(N,6,6,6)).astype(complex).copy()
        for axis in range(3):
            # Small unitary link mixes an active occupied and empty basis orbital.
            phase=.03*(axis+1)
            block=np.array([[np.cos(phase),1j*np.sin(phase)],[1j*np.sin(phase),np.cos(phase)]])
            for kk in range(N):dip[kk,axis][np.ix_([2,4],[2,4])]=block
        P=berry_polarization(dyn,dip,axes,lat,343.,2.,4,3,6)
        # Independent literal loop implementation of the uploaded determinant convention.
        V=full_occupied(dyn,4,3,6);red=np.zeros((len(tau),3))
        for a,info in axes.items():
            if info['N']==1:continue
            phase=np.empty((len(tau),len(info['cycles'])))
            for it in range(len(tau)):
                for j,cycle in enumerate(info['cycles']):
                    product=1.+0j
                    for kk in cycle:
                        kp=info['neighbor'][kk];det=np.linalg.det(V[it,kk].conj().T@dip[kk,a].T@V[it,kp])
                        product*=det/abs(det)
                    phase[it,j]=np.angle(product)
            red[:,a]=2/(2*np.pi)*np.unwrap(phase,axis=0).mean(axis=1)
        assert np.allclose(P,red@lat/343.,rtol=1e-11,atol=1e-15)
        qq={'kpoints_full_red':k,'grid_shape':np.array(shape),'varying_axes':np.array(vary)}
        vals=np.cos(2*np.pi*k[:,0])[:,None]
        assert np.allclose(periodic_interpolate(qq,vals,k),vals,atol=1e-12)
        passed(f'{len(vary)}D Berry determinant and full periodic-mesh interpolation')
    A=np.array([[5.,0,0],[-2.5,5*np.sqrt(3)/2,0],[0,0,20.]])
    B=np.linalg.inv(A).T
    bz,hull,frame=first_bz(B,[0,1]);assert len(bz)==6
    passed('hexagonal first-BZ geometry (six vertices, correct reciprocal area)')
    A=np.array([[0,3,3],[3,0,3],[3,3,0]],float);B=np.linalg.inv(A).T
    bz,hull,frame=first_bz(B,[0,1,2]);assert len(bz)==24
    passed('FCC first-BZ geometry (truncated-octahedron vertices and volume)')
    mirror=np.array([[0,1,0],[1,0,0],[0,0,1]])
    qq={'symmetry_cart':np.array([np.eye(3),mirror]),'symmetry_antiunitary':np.zeros(2,bool)}
    proj=polar_projector(qq,[1,1,0],'residual')
    rawP=np.column_stack([np.sin(tau)*1.01,np.sin(tau)*.99,np.cos(tau)*.001])
    symP=project_delta(rawP,proj);berry_validation(symP,symP,1e-8,1e-10)
    assert np.allclose(symP[:,2],rawP[:,2]) and np.allclose(np.diff(symP[:,:2],axis=1),np.diff(symP[0,:2]))
    passed('residual AlAs [110] projector preserves z and static Berry offsets')
    time=np.linspace(.3,period+.3,301,endpoint=False);known=np.array([.2+.1j,.03-.04j,.005+.006j])
    signal=np.column_stack([sum(2*np.real(c*np.exp(-1j*(j+1)*omega*time)) for j,c in enumerate(known))]*3)
    fit=harmonic_fit(time,signal,omega,6,0.)
    assert np.allclose(fit[1:4],known[:,None],atol=1e-13)
    passed('complex harmonic 1/2/3 extraction at fit order 6')
    fake=SimpleNamespace(IO_TIME_points=np.arange(100.),Polarization=[np.zeros((3,100))])
    fake.Polarization[0][1,-1]=.001
    try:load_polarization_at(fake,np.array([10.,50.,99.]))
    except CompareError:pass
    else:raise AssertionError('Sparse history was accepted')
    passed('completed-run sparse polarization is rejected, not zero-filled or bypassed')
    with tempfile.TemporaryDirectory() as tmp:
        p=Path(tmp)/'test.npz';atomic_npz(p,data=np.array([1+2j]),meta='safe')
        assert load_npz(p)['data'][0]==1+2j
        src=Path(tmp)/'source';src.write_text('old');st=source_stamp([src]);src.write_text('different')
        try:check_source_stamp(st)
        except CompareError:pass
        else:raise AssertionError('Stale source accepted')
    passed('pickle-free NPZ round trip and stale-source rejection')
    # Adaptive-frequency regression: complex response, peaks, spacing and bounds.
    e=np.linspace(.5,8.,601)
    yy=(.15+0j + (.04+.1j)/(e-2.6+.10j) + (.025-.04j)/(e-5.15+.20j))
    xx=np.zeros((len(e),3),complex);xx[:,1]=yy
    plan=adaptive_frequency_grid(e,xx,budget=72,component='y',max_step_eV=.35)
    chosen=plan['energies_eV']
    assert len(chosen)==72 and np.all(np.diff(chosen)>0)
    assert chosen[0]==e[0] and chosen[-1]==e[-1]
    assert np.max(np.diff(chosen))<=.35+1e-10
    assert np.any(np.abs(chosen-2.6)<.025) and np.any(np.abs(chosen-5.15)<.025)
    assert np.count_nonzero(abs(chosen-2.6)<.3)>np.count_nonzero(abs(np.linspace(.5,8.,72)-2.6)<.3)
    passed('complex adaptive frequency selection preserves peaks, endpoints, budget and maximum gap')
    same,indices,_=reference_on_floquet_grid(e,xx,e[::5])
    assert np.array_equal(indices,np.arange(0,len(e),5)) and np.allclose(same,xx[::5])
    offgrid=np.array([2.615,4.121])
    try:reference_on_floquet_grid(e,xx,offgrid,tolerance=1e-6)
    except CompareError:pass
    else:raise AssertionError('Off-grid reference was silently interpolated')
    pi,_,_=reference_on_floquet_grid(e,xx,offgrid,method='pchip')
    assert np.all(np.isfinite(pi))
    passed('RT subset matching is exact; off-grid complex interpolation requires explicit opt-in')
    # Cartesian SHG index selection, intrinsic j/k symmetry, and point-group constraints.
    s=parse_tensor_component('[1,2,1]',2)
    assert s['key']=='[1,1,2]' and s['label']=='chi_xyx'
    try:parse_tensor_component('[1,2]',2)
    except CompareError:pass
    else:raise AssertionError('Mismatched harmonic tensor index count accepted')
    passed('1-based tensor indices; input-index permutation; order validation')
    trivial=tensor_catalog(2,[0,1,0])
    assert trivial['components']['[2,2,2]']['status']=='extractable'
    assert trivial['components']['[1,1,2]']['status']=='not_identifiable'
    try:tensor_response(trivial,'[1,1,2]',np.ones((1,3),complex))
    except CompareError:pass
    else:raise AssertionError('Unmeasured chi_xxy was mislabeled as measurable')
    passed('One E||y trajectory cannot generically isolate chi_xxy')
    angle=2*np.pi/3
    C3=np.array([[np.cos(angle),-np.sin(angle),0],
                 [np.sin(angle),np.cos(angle),0],[0,0,1.]])
    Myz=np.diag([-1.,1.,1.]);Mxy=np.diag([1.,1.,-1.])
    ops=[np.linalg.matrix_power(C3,n)@m@z for n in range(3)
         for m in [np.eye(3),Myz] for z in [np.eye(3),Mxy]]
    d3=tensor_catalog(2,[0,1,0],np.array(ops),'synthetic D3h')
    assert d3['independent_tensor_parameters']==1
    sample=np.array([[.0,.3+.5j,.0]],complex)
    derived,_=tensor_response(d3,'[1,1,2]',sample)
    assert np.allclose(derived,-sample[:,1])
    assert d3['components']['[1,2,2]']['status']=='symmetry_forbidden'
    passed('Synthetic D3h: chi_xxy=-chi_yyy for y-directed monolayer SHG')
    td=[]
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((-1,1),repeat=3):
            if np.prod(signs)!=1:continue
            R=np.zeros((3,3))
            for i,j in enumerate(perm):R[i,j]=signs[i]
            td.append(R)
    cubic=tensor_catalog(2,[1,1,0],np.array(td),'synthetic Td')
    assert cubic['independent_tensor_parameters']==1
    sample=np.array([[0,0,.23-.4j]],complex)
    derived,_=tensor_response(cubic,'[3,2,1]',sample)
    assert np.allclose(derived,sample[:,2])
    assert cubic['components']['[1,1,2]']['status']=='symmetry_forbidden'
    passed('Synthetic Td [110]: chi_zxy=chi_z_eff, chi_xxy forbidden')
    inversion=tensor_catalog(2,[0,1,0],np.array([np.eye(3),-np.eye(3)]))
    assert inversion['independent_tensor_parameters']==0
    assert all(x['status']=='symmetry_forbidden' for x in inversion['components'].values())
    passed('Centrosymmetric point group forbids all electric-dipole SHG components')
    with_time_reversal=tensor_catalog(2,[0,1,0],np.array([np.eye(3),-np.eye(3)]),
                                     antiunitary=[False,True])
    assert with_time_reversal['independent_tensor_parameters']==18
    assert with_time_reversal['components']['[2,2,2]']['status']=='extractable'
    passed('Antiunitary time reversal is not mistaken for spatial inversion')
    print(f'{len(tests)} synthetic checks passed. Raw Lumen/HPC integration remains a separate test.')
    return tests


def main(argv=None):
    p=argparse.ArgumentParser(description='FLQ shared utilities: audit, frequencies, compare and self-test')
    sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('audit');a.add_argument('--workdir',default='.');a.add_argument('--calc-glob',default='FLQ_[0-9][0-9][0-9]')
    a.add_argument('--out-name',default='PyFLQ');a.add_argument('--expected-count',type=int);a.add_argument('--verify-sources',action='store_true')
    f=sub.add_parser('frequencies');f.add_argument('--workdir',default='.');f.add_argument('--rt-calc',default='NL')
    f.add_argument('--output',default='FLQ_frequencies.dat');f.add_argument('--expected-count',type=int)
    f.add_argument('--adaptive',action='store_true',help='Select original RT spectrum energies adaptively')
    f.add_argument('--off-grid',action='store_true',help='Allow new Floquet energies predicted by PCHIP interpolation of RT')
    f.add_argument('--candidate-step-eV',type=float,default=0.02,help='Candidate spacing when --off-grid is enabled')
    f.add_argument('--spectrum',help='RT 7-column YamboPy/FLQ_HG spectrum; default PyFLQ/RT_order_2.dat')
    f.add_argument('--points',type=int,default=48,help='Floquet frequency budget, including endpoints')
    f.add_argument('--component',default='y',
                   help='x/y/z/all; or [1,1,2] when a matching FLQ_HG NPZ supplies the symmetry tensor catalog')
    f.add_argument('--harmonic-order',type=int,default=2)
    f.add_argument('--max-step-eV',type=float,default=0.35)
    f.add_argument('--min-step-eV',type=float,default=0.0)
    f.add_argument('--peak-prominence',type=float,default=0.05)
    f.add_argument('--max-peak-anchors',type=int,default=8)
    f.add_argument('--slope-weight',type=float,default=0.08)
    f.add_argument('--curvature-weight',type=float,default=0.08)
    f.add_argument('--relative-floor',type=float,default=0.03)
    f.add_argument('--force',action='store_true',help='Allow changed frequency list only if no existing FLQ directories')
    c=sub.add_parser('compare');c.add_argument('--rt',required=True);c.add_argument('--fl',required=True)
    c.add_argument('--component',choices=['x','y','z'],default='y');c.add_argument('--energy-tol',type=float,default=1e-4)
    sub.add_parser('self-test')
    args=p.parse_args(argv)
    if args.command=='self-test':self_test();return
    if args.command=='compare':
        Er,R=read_spectrum(args.rt);Ef,F=read_spectrum(args.fl)
        print(json.dumps(jsafe(comparison_metrics(Er,R,Ef,F,'xyz'.index(args.component),args.energy_tol)),indent=2));return
    work=Path(args.workdir).expanduser().resolve()
    if args.command=='frequencies':
        if args.adaptive:
            if args.spectrum:
                source=Path(args.spectrum).expanduser()
                if not source.is_absolute():source=work/source
            else:
                possibilities=[work/'PyFLQ/RT_order_2.dat',work/'o-NL.YamboPy-X_probe_order_2']
                source=next((p for p in possibilities if p.is_file()),possibilities[0])
            if not source.is_file():
                raise CompareError(f'Missing RT source spectrum {source}; first run FLQ_HG.py --source rt')
            E,X=read_spectrum(source)
            if args.expected_count is not None and len(E)!=args.expected_count:
                raise CompareError(f'Expected {args.expected_count} RT source points, found {len(E)}')
            adaptive_component=args.component
            if args.component not in ('x','y','z','all'):
                selection=parse_tensor_component(args.component,args.harmonic_order)
                if selection['kind']!='tensor':raise CompareError('Invalid adaptive component')
                derived=source.parent/f'HG_order_{args.harmonic_order}.npz'
                if not derived.is_file():
                    raise CompareError(f'Named tensor sampling needs {derived}, produced by FLQ_HG.py')
                data=load_npz(derived)
                if 'tensor_catalog_json' not in data or 'rt_chi' not in data:
                    raise CompareError('The harmonic cache predates tensor diagnostics or has no RT spectrum')
                if not np.allclose(E,np.asarray(data['rt_energy_eV'],float),atol=1e-7,rtol=0):
                    raise CompareError('RT spectrum and harmonic cache have different energy grids')
                if not np.allclose(X,np.asarray(data['rt_chi'],complex),atol=1e-12,rtol=1e-7):
                    raise CompareError('RT spectrum and tensor cache disagree; rerun FLQ_HG.py')
                catalog=json.loads(str(data['tensor_catalog_json']))
                selected,identification=tensor_response(catalog,args.component,X)
                print('Adaptive tensor channel:',identification)
                X=np.zeros_like(X);X[:,0]=selected
                adaptive_component='x'
            source_E,source_X=E,X
            if args.off_grid:
                from scipy.interpolate import PchipInterpolator
                if args.candidate_step_eV<=0:
                    raise CompareError('--candidate-step-eV must be positive')
                n_candidates=int(np.ceil((E[-1]-E[0])/args.candidate_step_eV))+1
                if n_candidates>20000:raise CompareError('Candidate spacing generates >20000 points')
                E=np.unique(np.r_[E,np.linspace(E[0],E[-1],n_candidates)])
                X=PchipInterpolator(source_E,source_X.real,axis=0)(E)+1j*PchipInterpolator(
                    source_E,source_X.imag,axis=0)(E)
            plan=adaptive_frequency_grid(E,X,budget=args.points,component=adaptive_component,
                max_step_eV=args.max_step_eV,min_step_eV=args.min_step_eV,
                peak_prominence=args.peak_prominence,max_peak_anchors=args.max_peak_anchors,
                slope_weight=args.slope_weight,curvature_weight=args.curvature_weight,
                relative_floor=args.relative_floor)
            plan['component']=args.component
            frequencies=plan['energies_eV']
            summary={k:v for k,v in plan.items()
                     if k not in ('source_energy_eV','source_magnitude','pointwise_relative_interpolation_error')}
            dist_to_rt=np.min(np.abs(frequencies[:,None]-source_E[None,:]),axis=1)
            summary.update(schema=SCHEMA,source_spectrum=str(source.resolve()),
                           source_sha256=source_stamp([source])[str(source.resolve())]['sha256'],
                           method='PCHIP_surrogate_off_grid' if args.off_grid else 'actual_RT_points_only',
                           rt_reference_requirements='--rt-reference pchip' if args.off_grid else '--rt-reference exact',
                           actual_RT_point_count=int(len(source_E)),
                           offgrid_selected_count=int(np.sum(dist_to_rt>1e-6)),
                           original_RT_energy_eV=source_E,
                           original_RT_magnitude=np.linalg.norm(source_X if adaptive_component=='all' else source_X[:,['x','y','z'].index(adaptive_component),None],axis=1),
                           interpolation_reference='complex PCHIP surrogate' if args.off_grid else 'actual RT points only',
                           source_energy_eV=plan['source_energy_eV'],
                           source_magnitude=plan['source_magnitude'],
                           pointwise_relative_interpolation_error=plan['pointwise_relative_interpolation_error'])
        else:
            from yambopy import YamboNLDB
            db=YamboNLDB(folder=str(work),calc=args.rt_calc)
            frequencies=np.array([float(field['freq_range'][0])*HA_TO_EV for field in db.Efield])
            if not np.isfinite(frequencies).all() or np.any(frequencies<=0) or (
               len(frequencies)>1 and np.any(np.diff(frequencies)<=0)):
                raise CompareError('Invalid/duplicate RT frequencies')
            if args.expected_count is not None and len(frequencies)!=args.expected_count:
                raise CompareError('Unexpected RT frequency count')
            summary=None
        if len(frequencies)>1000:raise CompareError('FLQ_### folder naming supports at most 1000 points')
        output=work/args.output
        content=''.join(f'{i:03d} {e:.12f}\n' for i,e in enumerate(frequencies))
        existing=output.read_text() if output.exists() else None
        if existing!=content:
            prior=sorted(work.glob('FLQ_[0-9][0-9][0-9]'))
            if prior:
                raise CompareError('Cannot change FLQ_frequencies.dat while indexed FLQ_### directories exist. '
                                   'Their names may refer to OLD energies. Preserve/move those runs first.')
            if existing is not None and not args.force:
                raise CompareError(f'{output} already contains a different frequency grid; '
                                   'use --force after verifying no FLQ_### runs exist')
            atomic_text(output,content)
        if summary is not None:
            summary_file=output.with_suffix('.json')
            atomic_text(summary_file,json.dumps(jsafe(summary),indent=2,allow_nan=False)+'\n')
            print(f'Adaptive RT grid: {len(E)} source -> {len(frequencies)} Floquet points; '
                  f'complex interpolation residual max={plan["max_relative_interpolation_error"]:.4g}, '
                  f'frequency steps={plan["min_gap_eV"]:.5g}–{plan["max_gap_eV"]:.5g} eV; '
                  f'off-grid={summary["offgrid_selected_count"]}')
            print(f'Wrote plan metadata to {summary_file}')
        print(f'Wrote {len(frequencies)} frequencies to {output}');return
    with run_log(work/args.out_name/'FLQ_AL_audit.log'):
        calcs=calculation_dirs(work,args.calc_glob)
        if args.expected_count is not None and len(calcs)!=args.expected_count:raise CompareError('Unexpected calculation count')
        energies=[];failed=[]
        for calc in calcs:
            try:
                q,f,m=load_pair(calc,args.out_name,args.verify_sources)
                e=float(q['omega_au'])*HA_TO_EV;energies.append(e)
                print(f'{calc.name}: {e:.8f} eV; bands {m["band_window"]}; '
                      f'4T closure={np.max(q["closure_4T"]):.3e}; edge={np.max(f["edge_fraction"]):.3e}')
            except Exception as exc:failed.append(calc.name);print(f'{calc.name}: FAILED: {exc}')
        if failed:raise CompareError(f'{len(failed)} missing/invalid extractions: {failed}')
        if len(np.unique(np.round(energies,8)))!=len(energies):raise CompareError('Duplicate cached photon energies')
        print(f'Validated {len(calcs)} cache pairs; this is not physical-convergence certification')

if __name__=='__main__':main()
