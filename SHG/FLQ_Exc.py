#!/usr/bin/env python3
"""FLQ_Exc: TD-HSEX SHG / BSE exciton resonance and k-space analysis.

Two independent levels of evidence:
  (1) Sorted ypp exciton energies + ONE-PHOTON oscillator strengths identify
      candidate resonance windows (E_exc/2 and E_exc) in chi^(2)(2w;w,w).
  (2) Saved BSE EIGENSTATES give gauge-invariant k-resolved probability weights:
      W_group(k) = mean_{lambda in group} sum_{v,c,spin}|A_lambda(k,v,c,spin)|^2.
      Averaging/summing a degenerate subspace avoids basis-dependent maps.

Neither (1) nor (2) gives exciton-resolved SHG oscillator strengths or a
state-by-state decomposition of chi^(2). That requires additional coupling data.

Requires NumPy, Matplotlib, SciPy. Optional BZ stage requires YamboPy/netCDF4.
The BS eigenvectors must be stored as BS_EIGENSTATES in ndb.BS_diago_Q1.
No Floquet, QE, BSE, GW, or Lumen calculations are modified by this script.
"""
from __future__ import annotations
import argparse
import csv
import itertools
import json
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPLPath
from scipy.spatial import Voronoi

BOHR_ANG = 0.529177210903


def read_sorted(filename):
    path = Path(filename).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f'Sorted exciton file missing: {path}')
    a = np.atleast_2d(np.loadtxt(path, comments='#'))
    if a.shape[1] < 3 or len(a) < 2:
        raise ValueError('Expected ypp columns: E_exc[eV], normalized one-photon strength, 1-based state index')
    e, s, state_float = a[:, 0], a[:, 1], a[:, 2]
    if not np.all(np.isfinite(a[:, :3])) or np.any(e <= 0) or np.any(s < 0):
        raise ValueError('Nonfinite, nonpositive-energy, or negative-strength exciton data')
    idx = np.rint(state_float).astype(int)
    if np.any(abs(state_float-idx) > 1e-3) or np.any(idx < 1) or len(set(idx)) != len(idx):
        raise ValueError('Invalid or repeated 1-based exciton state indices')
    order = np.argsort(idx)
    return dict(energy=e[order], strength=s[order], index=idx[order], path=str(path))


def compare_sorted(energy, intensity):
    if len(energy['index']) != len(intensity['index']) or not np.array_equal(energy['index'], intensity['index']):
        raise ValueError('Energy- and intensity-sorted files identify different BSE states')
    if not np.allclose(energy['energy'], intensity['energy'], rtol=0, atol=5e-5):
        raise ValueError('Energy- and intensity-sorted files have different BSE eigenvalues')
    if not np.allclose(energy['strength'], intensity['strength'], rtol=1e-5, atol=1e-9):
        raise ValueError('Energy- and intensity-sorted files have different oscillator strengths')


def find_sorted(root, supplied, kind):
    if supplied:
        p = Path(supplied).expanduser()
        return p.resolve() if p.is_absolute() else (root / p).resolve()
    suffix = '_E_sorted' if kind == 'energy' else '_I_sorted'
    candidates = set()
    for folder in (root,root.parent):
        for p in folder.glob('*.exc_qpt1'+suffix): candidates.add(p.resolve())
        for p in folder.glob('BSE_files/*.exc_qpt1'+suffix): candidates.add(p.resolve())
    found = sorted(p for p in candidates if p.is_file())
    if len(found) != 1:
        raise FileNotFoundError(f'Found {len(found)} {kind}-sorted files in the standard calculation directories; specify --{kind}-sorted explicitly (avoid mixing different scissor runs).')
    return found[0]


def auto_groups(data, tol_eV, max_groups, max_states):
    """Build contiguous nearly-degenerate multiplets from sorted BSE eigenvalues."""
    order = np.argsort(data['energy'])
    energies = data['energy'][order]
    indices = data['index'][order]
    groups = []
    i = 0
    while i < min(len(order), max_states):
        j = i + 1
        while j < len(order) and energies[j]-energies[i] <= tol_eV:
            j += 1
        groups.append(indices[i:j].tolist())
        if len(groups) >= max_groups:
            break
        i = j
    return groups


def parse_groups(spec, data, tol_eV, max_groups, max_states):
    if spec.strip().lower() == 'auto':
        return auto_groups(data, tol_eV, max_groups, max_states)
    groups = []
    all_index = set(map(int, data['index']))
    seen = set()
    for item in spec.split(';'):
        g = [int(x.strip()) for x in item.split(',') if x.strip()]
        if not g or not set(g) <= all_index or any(x in seen for x in g):
            raise ValueError(f'Invalid or overlapping exciton group: {item!r}')
        groups.append(g)
        seen.update(g)
    # Do not accidentally split a near-degenerate physical subspace.
    emap = dict(zip(data['index'], data['energy']))
    for group in groups:
        for n in group:
            close = [int(other) for other, e in emap.items() if abs(e - emap[n]) <= tol_eV]
            if not set(close) <= set(group):
                raise ValueError(f'Group {group} splits exciton {n} from degenerate partner(s) {close}; '
                                 'include the full multiplet or reduce --degeneracy-tol-eV')
    return groups


def write_csv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def select_shg(path, component, source):
    if path is None:
        return None
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(p)
    axis = 'xyz'.index(component)
    units = 'unknown'
    if p.suffix.lower() == '.npz':
        with np.load(p, allow_pickle=False) as db:
            kinds = ('rt', 'fl') if source == 'auto' else (source,)
            names = next(((k+'_energy_eV', k+'_chi', k) for k in kinds
                          if k+'_energy_eV' in db and k+'_chi' in db), None)
            if names is None:
                raise ValueError('NPZ needs rt_energy_eV/rt_chi or fl_energy_eV/fl_chi')
            E = np.asarray(db[names[0]], float)
            chi = np.asarray(db[names[1]], complex)
            if chi.shape != (len(E), 3):
                raise ValueError('SHG NPZ susceptibility must have shape (Nfreq,3)')
            value = chi[:, axis]
            origin = names[2].upper()
            units = str(db['units']) if 'units' in db else 'unknown'
    else:
        a = np.atleast_2d(np.loadtxt(p, comments='#'))
        if a.shape[1] != 7:
            raise ValueError('SHG text must have 7 columns: E, Im_x,Re_x, Im_y,Re_y, Im_z,Re_z')
        E = a[:, 0]
        value = a[:, 2+2*axis] + 1j*a[:, 1+2*axis]
        origin = 'tabulated'
        units = 'cm/statV' # expected for this SHG 7-column output; verify file header
    mask = np.argsort(E)
    E, value = E[mask], value[mask]
    if not np.all(np.isfinite(E)) or not np.all(np.isfinite(value)) or len(np.unique(E)) != len(E):
        raise ValueError('SHG spectrum has invalid/duplicate energies or nonfinite values')
    return E, value, origin, units


def plot_catalog(data, groups, out, dark_threshold, shg=None, show_count=40):
    E = data['energy']; s = data['strength']; index = data['index']
    sorted_order = np.argsort(E)[:show_count]
    brightest = np.argsort(-s)[:10]
    fig, ax = plt.subplots(figsize=(10.5, 4.5))
    for j in sorted_order:
        ax.vlines(E[j], 0, s[j], linewidth=1.6 if s[j] > dark_threshold else 0.8,
                  linestyles='-' if s[j] > dark_threshold else ':')
    for j in brightest:
        if E[j] <= max(E[sorted_order]) and s[j] > 1e-3:
            ax.annotate(str(int(index[j])), (E[j], s[j]), xytext=(0, 4),
                        textcoords='offset points', ha='center', fontsize=8)
    ax.set(xlabel=r'Exciton energy $E_{\lambda}$ (eV)',
           ylabel='Normalized one-photon oscillator strength',
           title='BSE excitons at optical momentum Q = 0 (polarization-specific)')
    ax.grid(alpha=0.18); fig.tight_layout()
    fig.savefig(out/'01_BSE_exciton_strengths.png', dpi=230)
    plt.close(fig)

    fig, (ax, bx) = plt.subplots(2, 1, figsize=(10.7, 8.0))
    for group in groups:
        loc = [np.flatnonzero(index == x)[0] for x in group]
        Eg = float(np.mean(E[loc])); bright = float(np.sum(s[loc]))
        description = ','.join(map(str, group))
        ax.scatter(Eg, bright, s=32)
        ax.annotate(description, (Eg, bright), xytext=(5, 3), textcoords='offset points', fontsize=8)
    ax.set(xlabel=r'Exciton energy $E_{\lambda}$ (eV)',
           ylabel='Group summed linear oscillator strength',
           title='BSE exciton multiplets (nearly degenerate states grouped)')
    ax.grid(alpha=.2)

    if shg is None:
        # The markers indicate only possible RESONANCE LOCATIONS. No fabricated chi^(2).
        for g in groups:
            loc = [np.flatnonzero(index == x)[0] for x in g]
            Eg = float(np.mean(E[loc])); bright = float(np.max(s[loc]))
            bx.scatter([Eg/2, Eg], [bright, bright], marker='|' if bright > dark_threshold else 'x', s=95)
        bx.set(ylabel='Linear strength (resonance guide only)',
               title=r'Possible SHG driving-energy resonances (no $\chi^{(2)}$ spectrum supplied)')
    else:
        omega, chi, origin, units = shg
        bx.plot(omega, abs(chi), label=rf'$|\chi^{{(2)}}_{{\mathrm{{eff}}}}|$ ({origin})')
        bx.set(ylabel=rf'$|\chi^{{(2)}}_{{\mathrm{{eff}}}}|$ ({units})',
               title='Measured/calculated SHG and possible exciton resonances')
        bx.legend()
    for q in (.5,1.):
        for g in groups:
            loc = [np.flatnonzero(index == x)[0] for x in g]
            Eg = float(np.mean(E[loc]))
            bx.axvline(Eg*q, linewidth=.8, linestyle='--' if q==.5 else ':', alpha=.35)
    # Dashed guides: 2*w=E_exc; dotted guides: w=E_exc.
    from matplotlib.lines import Line2D
    handles = [Line2D([0],[0],ls='--',lw=1,label=r'$2\hbar\omega=E_\lambda$'),
               Line2D([0],[0],ls=':',lw=1,label=r'$\hbar\omega=E_\lambda$')]
    bx.legend(handles=handles,loc='upper right',fontsize=9)
    bx.set_xlabel(r'Driving photon energy $\hbar\omega$ (eV)')
    bx.grid(alpha=.18)
    fig.tight_layout();fig.savefig(out/'02_SHG_resonance_guide.png', dpi=230);plt.close(fig)



def shg_peak_matches(data, groups, shg, out, window_eV=0.25, min_prominence_relative=0.02):
    """Flag spectral coincidences; NOT a quantitative exciton contribution."""
    from scipy.signal import find_peaks, peak_prominences
    E,chi,origin,units=shg
    amplitude=np.abs(chi)
    if len(E)<3: raise ValueError('At least 3 SHG energy nodes required for peak detection')
    min_prom=float(np.max(amplitude))*min_prominence_relative
    ids, properties=find_peaks(amplitude,prominence=min_prom)
    emap=dict(zip(map(int,data['index']),data['energy']))
    rows=[]
    for states in groups:
        avg=float(np.mean([emap[n] for n in states]))
        for channel,factor in (('2omega=exciton',.5),('omega=exciton',1.)):
            expected=avg*factor
            eligible=[i for i in ids if abs(E[i]-expected)<=window_eV]
            selected=min(eligible,key=lambda i: abs(E[i]-expected)) if eligible else None
            rows.append([
                ','.join(map(str,states)), avg,channel,expected,
                float(E[selected]) if selected is not None else '',
                float(E[selected]-expected) if selected is not None else '',
                float(amplitude[selected]) if selected is not None else '',
                origin,units,
                'candidate_only; not BSE-resolved SHG'
            ])
    write_csv(out/'SHG_peak_alignment_candidates.csv',
        ['exciton_states','exciton_E_eV','resonance_channel','predicted_pump_E_eV',
         'closest_SHG_peak_E_eV','detuning_eV','chi_magnitude_at_peak',
         'SHG_source','SHG_units','interpretation'],rows)
    print(f'SHG peak guide: {len(ids)} local peaks exceed {min_prominence_relative*100:g}% global prominence; '
          f'search window +/-{window_eV:g} eV; matches are NOT causal attributions')
    return rows


def first_bz_geometry(reciprocal_bohr_inv):
    """2D Wigner-Seitz reciprocal cell in Å^-1, generated from the actual lattice."""
    B = np.asarray(reciprocal_bohr_inv, float) * (2*np.pi/BOHR_ANG)
    b1, b2 = B[0], B[1]
    e1 = b1 / np.linalg.norm(b1)
    perpendicular = np.cross(b1, b2)
    if np.linalg.norm(perpendicular) < 1e-8:
        raise ValueError('Reciprocal basis not two-dimensional')
    perpendicular /= np.linalg.norm(perpendicular)
    e2 = np.cross(perpendicular, e1)
    frame = np.column_stack([e1, e2])
    primitive = np.stack([b1, b2]) @ frame
    shifts = np.array(list(itertools.product(range(-3,4), repeat=2)), dtype=int)
    p = shifts @ primitive
    center = np.flatnonzero(np.all(shifts == 0, axis=1))
    vor = Voronoi(p)
    region = vor.regions[vor.point_region[int(center[0])]]
    if -1 in region or len(region) < 3:
        raise ValueError('Reciprocal Wigner-Seitz construction failed')
    vertices = vor.vertices[region]
    angles = np.arctan2(vertices[:,1],vertices[:,0])
    vertices = vertices[np.argsort(angles)]
    return B,frame,vertices


def first_bz_data(red_kpoints, weights, reciprocal_bohr_inv):
    B, frame, poly = first_bz_geometry(reciprocal_bohr_inv)
    red = np.asarray(red_kpoints, float)
    nk = len(red)
    if len(weights) != nk or red.shape != (nk,3):
        raise ValueError('Exciton weight array and expanded lattice k-points disagree')
    translations = np.array(list(itertools.product(range(-2,3), repeat=2)), int)
    shift = np.zeros((len(translations),3))
    shift[:, :2] = translations
    candidates_red = (red[None,:,:] + shift[:,None,:]).reshape(-1,3)
    candidates_xy = (candidates_red @ B) @ frame
    parent = np.tile(np.arange(nk), len(shift))
    inside = MPLPath(np.vstack([poly,poly[0]]), closed=True).contains_points(candidates_xy, radius=2e-6)
    xy = candidates_xy[inside]
    ids = parent[inside]
    # Boundary translations can show the same point twice. Keep unique coordinates.
    _, keep = np.unique(np.round(xy,8),axis=0,return_index=True)
    keep.sort()
    xy,ids = xy[keep],ids[keep]
    return xy, np.asarray(weights)[ids], poly, ids


def state_weight_and_transitions(exc, states, nk):
    if exc.eigenvectors is None or exc.table is None:
        raise RuntimeError('No BS_EIGENSTATES or BS_TABLE in ndb.BS_diago_Q1. '
                           'Check with ncdump -h BSE/ndb.BS_diago_Q1. '
                           'To generate missing eigenvectors, rerun BSE with WRbsWF (new tag).')
    table = np.asarray(exc.table)
    if table.ndim != 2 or table.shape[1] < 3:
        raise ValueError(f'Unexpected BS_TABLE shape {table.shape}')
    k = np.rint(table[:,0]).astype(int)-1
    v = np.rint(table[:,1]).astype(int)
    c = np.rint(table[:,2]).astype(int)
    if np.any(k < 0) or np.any(k >= nk):
        raise ValueError(f'BS_TABLE k indices extend outside expanded SAVE k-list (N={nk}): '
                         f'kmin={k.min()+1}, kmax={k.max()+1}')
    eig = np.asarray(exc.eigenvectors)
    if eig.ndim != 2 or eig.shape[1] != len(k) or max(states)>len(eig):
        raise ValueError(f'BS eigenvectors have shape {eig.shape}; expected (>=maxstate,{len(k)})')
    result = np.zeros(nk)
    vc = {}
    for state in states:
        coefficients = np.asarray(eig[state-1], complex)
        wt = abs(coefficients)**2
        norm = float(np.sum(wt))
        if not 0.99 <= norm <= 1.01:
            raise ValueError(f'Exciton {state} eigenvector norm {norm:.8g}; unexpected resonant-BSE normalization')
        wt /= norm
        result += np.bincount(k, weights=wt, minlength=nk)
        for (vv,cc) in sorted(set(zip(v,c))):
            mask = (v == vv) & (c == cc)
            vc[(int(vv),int(cc))] = vc.get((int(vv),int(cc)),0.)+float(np.sum(wt[mask]))
    # Mean over degenerate states: gauge independent and normalized sum_k W = 1.
    result /= len(states)
    vc = {key:val/len(states) for key,val in vc.items()}
    if not np.isclose(sum(result),1,atol=1e-8):
        raise ValueError('k-resolved spectral weight did not normalize to 1')
    return result,vc



def parse_hs_csv(root, spec):
    """Read the user-specified high-symmetry vertices, no assumed hexagonal axes."""
    p=Path(spec).expanduser()
    p=(p if p.is_absolute() else root/p).resolve()
    if not p.is_file():
        raise FileNotFoundError(f'High-symmetry path does not exist: {p}; pass --hs with its correct location')
    with p.open(newline='') as fh:
        rd=csv.DictReader(fh)
        if not {'Label','k1','k2','k3'} <= set(rd.fieldnames or []):
            raise ValueError(f'{p} must contain Label,k1,k2,k3 columns')
        rows=list(rd)
    if len(rows)<2:raise ValueError('High-symmetry path needs >=2 vertices')
    labels=['Γ' if x['Label'].strip().upper() in ('G','GAMMA','Γ') else x['Label'].strip() for x in rows]
    points=np.array([[float(x[k]) for k in ('k1','k2','k3')] for x in rows],float)
    if not np.all(np.isfinite(points)):
        raise ValueError('Nonfinite high-symmetry k coordinates')
    return labels,points,p


def group_record(g,data,threshold):
    lookup={int(i):(float(e),float(s)) for i,e,s in zip(data['index'],data['energy'],data['strength'])}
    energy=np.array([lookup[i][0] for i in g]);s=np.array([lookup[i][1] for i in g])
    return {'states':list(map(int,g)), 'label':','.join(map(str,g)),
            'mean_energy_eV':float(np.mean(energy)), 'splitting_eV':float(np.ptp(energy)),
            'max_linear_strength':float(np.max(s)), 'sum_linear_strength':float(np.sum(s)),
            'weak_in_one_photon_y':bool(np.all(s<threshold)),
            'possible_2ph_pump_eV':float(np.mean(energy)*.5),
            'possible_1ph_pump_eV':float(np.mean(energy))}


def summarize_shg(shg,args,out,data,groups):
    """Peak comparison is evidence for detuning only, never exciton-resolved χ(2)."""
    from scipy.signal import find_peaks
    E,chi,origin,units=shg
    mag=np.abs(chi);peak_ids,_=find_peaks(mag,prominence=args.shg_peak_prominence_rel*mag.max())
    write_csv(out/'SHG_local_peaks.csv', ['pump_energy_eV','Re_chi','Im_chi','abs_chi','source','units'],
              [[float(E[i]),float(chi[i].real),float(chi[i].imag),float(mag[i]),origin,units] for i in peak_ids])
    matches=[]
    for g in groups:
        r=group_record(g,data,args.dark_threshold)
        for kind,target in [('two_photon',r['possible_2ph_pump_eV']),('one_photon',r['possible_1ph_pump_eV'])]:
            if target<E[0] or target>E[-1]:continue
            nearest=int(np.argmin(abs(E-target)))
            in_window=[int(i) for i in peak_ids if abs(E[i]-target)<=args.peak_window_eV]
            closest=min(in_window,key=lambda i:abs(E[i]-target)) if in_window else None
            matches.append([r['label'],r['mean_energy_eV'],kind,target,
                r['sum_linear_strength'],int(r['weak_in_one_photon_y']),
                float(E[nearest]),float(mag[nearest]),
                float(E[closest]) if closest is not None else '',
                float(E[closest]-target) if closest is not None else '',
                float(mag[closest]) if closest is not None else '',
                'coincidence only; NOT exciton-resolved SHG'])
    write_csv(out/'SHG_exciton_resonance_matches.csv',
              ['states','exciton_energy_eV','resonance_condition','predicted_pump_eV',
               'sum_linear_oscillator_strength','one_photon_weak_group','nearest_sample_eV',
               'SHG_magnitude_nearest_sample','nearest_detected_peak_eV','detuning_peak_minus_prediction_eV',
               'SHG_magnitude_detected_peak','interpretation'],matches)
    # Search a narrower physically motivated energy interval for nominally dark states.
    dark=[]
    for i,e,s in zip(data['index'],data['energy'],data['strength']):
        if s>=args.dark_threshold or not args.dark_energy_min<=e<=args.dark_energy_max:continue
        predicted=.5*float(e)
        if not E[0]<=predicted<=E[-1]:continue
        nearest=int(np.argmin(abs(E-predicted)))
        if len(peak_ids):
            ip=min(peak_ids,key=lambda j:abs(E[j]-predicted))
            distance=float(E[ip]-predicted)
        else:ip=None;distance=float('nan')
        dark.append([int(i),float(e),float(s),predicted,float(E[nearest]),float(mag[nearest]),
                     float(E[ip]) if ip is not None else '',distance if ip is not None else '',
                     int(ip is not None and abs(distance)<=args.peak_window_eV)])
    dark.sort(key=lambda row: (not row[-1], abs(row[-2]) if row[-2]!='' else np.inf,row[1]))
    write_csv(out/'dark_candidates_SHG.csv',
              ['state','exciton_energy_eV','one_photon_strength','two_photon_pump_eV',
               'nearest_sample_pump_eV','SHG_magnitude_at_pump','closest_SHG_peak_pump_eV',
               'peak_minus_prediction_eV','within_peak_window'],dark)
    # Small separation from a bright multiplet relative to a broadening scale means
    # spectral overlap; no unique dark-state attribution can be made.
    bright_records=[group_record(g,data,args.dark_threshold) for g in groups]
    bright_records=[r for r in bright_records if r['max_linear_strength']>=args.bright_reference_threshold]
    resolving=[]
    for d in dark:
        pump=d[3]
        if not bright_records:break
        b=min(bright_records,key=lambda r:abs(r['possible_2ph_pump_eV']-pump))
        separation=abs(pump-b['possible_2ph_pump_eV'])
        resolving.append([d[0],d[1],pump,b['label'],b['possible_2ph_pump_eV'],
                          separation,separation/args.broadening_reference_eV,
                          int(separation<args.broadening_reference_eV),
                          'Comparison to a user-supplied broadening scale; not a fitted linewidth'])
    write_csv(out/'dark_vs_bright_resolvability.csv',
              ['dark_state','dark_exciton_E_eV','dark_two_photon_pump_eV',
               'nearest_bright_states','bright_two_photon_pump_eV',
               'two_photon_pump_separation_eV','separation_over_broadening',
               'closer_than_reference_broadening','note'],resolving)
    print(f'SHG source {origin}: {len(E)} samples, {len(peak_ids)} resolved magnitude peaks; '
          f'{len(dark)} one-photon-weak candidates in [{args.dark_energy_min:g},{args.dark_energy_max:g}] eV')
    return peak_ids,dark


def plot_dark_spectrum(data,groups,shg,out,args,fl=None):
    E,chi,origin,units=shg
    fig,ax=plt.subplots(2,1,figsize=(10.4,7.8),sharex=True,
                         gridspec_kw={'height_ratios':[2,1.1]})
    a,b=ax
    a.plot(E,abs(chi),lw=1.6,label=f'{origin} |χ|')
    if fl is not None:
        ef,xf,_,_=fl
        a.plot(ef,abs(xf),'.',ms=3.6,label='Floquet |χ|')
    b.plot(E,chi.real,lw=1.25,label='Re χ')
    b.plot(E,chi.imag,lw=1.25,label='Im χ')
    if fl is not None:
        b.plot(ef,xf.real,'.',ms=2.7,alpha=.6)
        b.plot(ef,xf.imag,'.',ms=2.7,alpha=.6)
    for group in groups:
        info=group_record(group,data,args.dark_threshold)
        weak=info['weak_in_one_photon_y']
        for f,factor in [('2ω',.5),('1ω',1.)]:
            x=info['mean_energy_eV']*factor
            if not E[0]<=x<=E[-1]:continue
            for axis in ax:
                axis.axvline(x,color='tab:purple' if weak else 'tab:orange',alpha=.34,
                             ls=':' if weak else '--',lw=1.)
    a.set(ylabel=rf'$|\chi^{{(2)}}_{{\mathrm{{eff}},y}}|$ ({units})',
          title='TD-HSEX SHG and candidate BSE resonance positions')
    b.set(xlabel=r'Driving photon energy $\hbar\omega$ (eV)',
          ylabel=rf'$\chi^{{(2)}}_{{\mathrm{{eff}},y}}$ ({units})')
    a.legend(fontsize=9);b.legend(fontsize=9)
    for axis in ax:axis.grid(alpha=.15)
    fig.tight_layout();fig.savefig(out/'02_SHG_BSE_resonances.png',dpi=240);plt.close(fig)

    relevant=[group_record(g,data,args.dark_threshold) for g in groups
              if args.dark_zoom_min<=group_record(g,data,args.dark_threshold)['possible_2ph_pump_eV']<=args.dark_zoom_max]
    fig,(ax,bx)=plt.subplots(2,1,figsize=(11,7.5),sharex=True)
    mask=(E>=args.dark_zoom_min)&(E<=args.dark_zoom_max)
    if np.count_nonzero(mask)<3:
        plt.close(fig);return
    ax.plot(E[mask],abs(chi[mask]),lw=1.8,label=f'{origin} |χ|')
    bx.plot(E[mask],chi[mask].real,lw=1.5,label='Re χ')
    bx.plot(E[mask],chi[mask].imag,lw=1.5,label='Im χ')
    for record in relevant:
        pos=record['possible_2ph_pump_eV'];weak=record['weak_in_one_photon_y']
        color='tab:purple' if weak else 'tab:orange'
        for axis in (ax,bx):axis.axvline(pos,alpha=.65,color=color,ls=':' if weak else '--',lw=1.1)
        ax.annotate(f"{record['label']}"+(' (weak)' if weak else ' (bright)'),
                    (pos,1),xycoords=('data','axes fraction'),rotation=90,
                    xytext=(3,-7),textcoords='offset points',va='top',fontsize=8,color=color)
    ax.set_ylim(bottom=0);ax.set_ylabel(r'$|\chi^{(2)}_{\mathrm{eff},y}|$'+f' ({units})')
    bx.set_ylabel(r'Complex $\chi^{(2)}_{\mathrm{eff},y}$'+f' ({units})')
    bx.set_xlabel(r'Driving photon energy $\hbar\omega$ (eV)')
    ax.set_xlim(args.dark_zoom_min,args.dark_zoom_max)
    ax.set_title('Near-resonant SHG: purple = weak in one-photon BSE; orange = bright')
    for axis in (ax,bx):axis.grid(alpha=.14)
    ax.legend();bx.legend();fig.tight_layout()
    fig.savefig(out/'03_dark_exciton_SHG_zoom.png',dpi=240);plt.close(fig)


def plot_ibz_group(ibz_red,ibz_weights,states,out,stars,hs_labels,hs_vertices):
    descriptor='_'.join(f'{s:04d}' for s in states)
    fig,ax=plt.subplots(figsize=(7,6.3))
    sc=ax.scatter(ibz_red[:,0],ibz_red[:,1],c=ibz_weights,marker='h',s=100,
                  cmap='viridis',edgecolors='none')
    # The user-provided symmetry path is an explicit guide, never derived from an assumed crystal orientation.
    ax.plot(hs_vertices[:,0],hs_vertices[:,1],color='0.25',ls='--',lw=1,alpha=.75)
    for name,vert in zip(hs_labels,hs_vertices):
        ax.annotate(name,(vert[0],vert[1]),xytext=(4,4),textcoords='offset points',fontsize=9)
    ax.set_aspect('equal')
    ax.set_xlabel(r'Irreducible reduced $k_1$')
    ax.set_ylabel(r'Irreducible reduced $k_2$')
    ax.set_title(f'IBZ star-integrated BSE weight, states {states}\n'
                 r'$W_\lambda(K)=\sum_{k\in\mathrm{star}(K)}\sum_{vc}|A^\lambda_{vck}|^2$')
    fig.colorbar(sc,ax=ax,label='Fraction of total exciton weight per symmetry star')
    ax.grid(alpha=.12)
    fig.tight_layout();fig.savefig(out/f'04_IBZ_states_{descriptor}.png',dpi=240);plt.close(fig)


def nearest_periodic_kmesh(kfull,path_red,reciprocal_bohr_inv):
    from scipy.spatial import cKDTree
    b=np.asarray(reciprocal_bohr_inv,float)* (2*np.pi/BOHR_ANG)
    base=np.mod(np.asarray(kfull,float),1.)
    planar=np.array(list(itertools.product([-1,0,1],repeat=2)),float)
    shifts=np.zeros((len(planar),3),float);shifts[:,:2]=planar
    mesh=(base[None]+shifts[:,None]).reshape(-1,3)@b
    ref=np.tile(np.arange(len(base)),len(shifts))
    index=cKDTree(mesh)
    dist,j=index.query(np.asarray(path_red,float)@b)
    return ref[j],dist


def sampled_hs_path(lat,weights,labels,verts,nseg):
    b=np.asarray(lat.rlat,float)*(2*np.pi/BOHR_ANG)
    pieces=[];distance=[];ticks=[0.];cumulative=0.
    for va,vb in zip(verts[:-1],verts[1:]):
        line=np.linspace(va,vb,nseg,endpoint=False)
        edge=float(np.linalg.norm((vb-va)@b))
        pieces.append(line)
        distance.extend(cumulative+edge*np.arange(nseg)/nseg)
        cumulative+=edge;ticks.append(cumulative)
    path=np.vstack([*pieces,verts[-1:]])
    d=np.array([*distance,cumulative])
    indexes,misfit=nearest_periodic_kmesh(lat.red_kpoints,path,lat.rlat)
    return {'path_red':path,'path_distance':d,'sample_index':indexes,
            'sample_misfit_Ainv':misfit,'values':np.asarray(weights)[indexes],
            'labels':labels,'ticks':ticks}


def plot_sampled_path(obj,states,out):
    descriptor='_'.join(f'{s:04d}' for s in states)
    fig,(ax,bx)=plt.subplots(2,1,figsize=(9,6.4),sharex=True,
                              gridspec_kw={'height_ratios':[3,1]})
    ax.step(obj['path_distance'],obj['values'],where='mid',lw=1.5)
    ax.scatter(obj['path_distance'],obj['values'],s=7)
    bx.plot(obj['path_distance'],obj['sample_misfit_Ainv'],lw=1.2)
    for a in (ax,bx):
        for t in obj['ticks']:a.axvline(t,color='0.65',ls=':',lw=.8)
        a.grid(alpha=.13)
    ax.set_ylabel('Exciton $W(k)$ (per sampled full-BZ k)')
    bx.set_ylabel('Nearest-k\ndistance ($\\AA^{-1}$)')
    bx.set_xlabel('High-symmetry path');bx.set_xticks(obj['ticks'],obj['labels'])
    ax.set_title(f'States {states}: nearest-mesh BSE weights along provided HS path\n'
                 'Discrete sampling, not interpolated exciton dispersion')
    fig.tight_layout();fig.savefig(out/f'06_HS_path_states_{descriptor}.png',dpi=240);plt.close(fig)


def plot_first_bz_group(lat,weights,states,labels,verts,out,marker_size=95):
    xy,w,polygon,_=first_bz_data(lat.red_kpoints,weights,lat.rlat)
    descriptor='_'.join(f'{s:04d}' for s in states)
    fig,ax=plt.subplots(figsize=(7.7,6.5))
    sc=ax.scatter(xy[:,0],xy[:,1],c=w,marker='h',s=marker_size,
                  vmin=0,cmap='viridis',rasterized=True)
    bound=np.vstack([polygon,polygon[0]]);ax.plot(bound[:,0],bound[:,1],lw=1.3,color='0.3')
    # Locate HS vertices on the Wigner-Seitz cell (not arbitrary polygon vertices).
    B,frame,_=first_bz_geometry(lat.rlat)
    planar=np.array(list(itertools.product(range(-2,3),repeat=2)),float)
    shifts=np.zeros((len(planar),3));shifts[:,:2]=planar
    for name,k in zip(labels,verts):
        images=((k[None]+shifts)@B)@frame
        # Prefer the image inside WS, then the smallest cartesian momentum.
        inside=MPLPath(bound,closed=True).contains_points(images,radius=1e-6)
        valid=images[inside] if inside.any() else images
        point=valid[np.argmin(np.linalg.norm(valid,axis=1))]
        ax.scatter(*point,color='black',s=18,zorder=5)
        ax.annotate(name,point,xytext=(5,4),textcoords='offset points',fontsize=10)
    ax.set_aspect('equal');ax.set(xlabel=r'$k_x$ ($\AA^{-1}$)',ylabel=r'$k_y$ ($\AA^{-1}$)',
                                   title=f'Exciton states {states}: full first BZ transition weight')
    fig.colorbar(sc,ax=ax,label=r'$W(k)=\langle\sum_{vc}|A_{vck}|^2\rangle_{\lambda}$')
    fig.tight_layout();fig.savefig(out/f'05_BZ_states_{descriptor}.png',dpi=240);plt.close(fig)
    return xy,w,polygon


def plot_dark_bright_ibz(maps,out):
    """Compare one bright multiplet and one one-photon-weak candidate on the same scale."""
    bright=next((m for m in maps if 2 in m['states']),next((m for m in maps if not m['weak']),None))
    dark=next((m for m in maps if m['weak'] and 5 in m['states']),
              next((m for m in maps if m['weak']),None))
    if bright is None or dark is None:return None
    if len(bright['ibz'])!=len(dark['ibz']):return None
    a=bright['ibz'];b=dark['ibz'];XY=bright['ibz_red']
    fig,axes=plt.subplots(1,3,figsize=(15.2,4.7),constrained_layout=True)
    scale=max(float(np.max(a)),float(np.max(b)),1e-300)
    for ax,values,title in zip(axes[:2],[a,b],
                               [f'Bright states {bright["states"]}',f'Weak states {dark["states"]}']):
        img=ax.scatter(XY[:,0],XY[:,1],c=values,vmin=0,vmax=scale,cmap='viridis',marker='h',s=95)
        ax.set_title(title);ax.set_aspect('equal')
        ax.set_xlabel('IBZ reduced k1');ax.set_ylabel('IBZ reduced k2')
    fig.colorbar(img,ax=axes[:2],shrink=.8,label='Star-integrated exciton probability')
    delta=b-a;d=max(float(np.max(abs(delta))),1e-300)
    img2=axes[2].scatter(XY[:,0],XY[:,1],c=delta,vmin=-d,vmax=d,cmap='coolwarm',marker='h',s=95)
    axes[2].set_aspect('equal');axes[2].set_xlabel('IBZ reduced k1');axes[2].set_ylabel('IBZ reduced k2')
    axes[2].set_title('Weak minus bright (normalized state weights)')
    fig.colorbar(img2,ax=axes[2],shrink=.8,label='Difference in star weight')
    fig.savefig(out/'07_dark_vs_bright_IBZ.png',dpi=240);plt.close(fig)
    return (bright['states'],dark['states'])


def plot_bz_groups(root,args,data,groups,out,hs):
    try:
        from yambopy import YamboLatticeDB, YamboExcitonDB
    except ImportError as err:
        raise RuntimeError('YamboPy not importable; activate yambo-lelphc and its PYTHONPATH') from err
    save=resolve_at_root(root,args.save)
    bse=resolve_at_root(root,args.bse_dir)
    db=save/'ns.db1';wfdb=bse/args.bs_file
    if not db.is_file() or not wfdb.is_file():
        raise FileNotFoundError(f'Need {db} and {wfdb}; the sorted ypp text files cannot provide k-resolved eigenvectors.')
    lat_ibz=YamboLatticeDB.from_db_file(filename=str(db),Expand=False)
    lat=YamboLatticeDB.from_db_file(filename=str(db),Expand=True)
    nfull=len(lat.red_kpoints);nibz=len(lat_ibz.red_kpoints)
    if args.expected_kpoints and nfull!=args.expected_kpoints:
        raise ValueError(f'Full electronic BZ has {nfull} k points, expected {args.expected_kpoints}; '
                         'wrong SAVE or mesh, do not mix IPA with TD-HSEX')
    mapping=np.asarray(lat.BZ_to_IBZ_indexes,int)
    if len(mapping)!=nfull or mapping.min()<0 or mapping.max()>=nibz:
        raise ValueError('Corrupt BZ-to-IBZ mapping in ns.db1')
    stars=np.bincount(mapping,minlength=nibz)
    if np.any(stars==0):raise ValueError('One or more irreducible k points have empty stars')
    requested=max(max(g) for g in groups)
    exc=YamboExcitonDB.from_db_file(lat,filename=args.bs_file,folder=str(bse),
                                    Load_WF=True,neigs=requested)
    if exc.eigenvectors is None:
        raise RuntimeError(f'{wfdb} lacks BS_EIGENSTATES: enable WRbsWF in a new BSE run, '
                           'then rerun postprocessing. Do not delete current BSE/TD-HSEX/Floquet outputs.')
    if np.asarray(exc.eigenvectors).ndim!=2:
        raise ValueError('Only resonant/Tamm-Dancoff BSE eigenvectors in (state,transition) layout are supported')
    if len(exc.eigenvectors)<requested:raise ValueError('BSE DB contains fewer excitonic eigenvectors than requested')
    energies=dict(zip(map(int,data['index']),map(float,data['energy'])))
    for state in set(n for g in groups for n in g):
        observed=float(np.real(exc.eigenvalues[state-1]))
        if abs(energies[state]-observed)>args.energy_match_tol_eV:
            raise ValueError(f'Exciton {state}: sorted E={energies[state]:.8f} eV but DB E={observed:.8f} eV. '
                             'BSE scissor/database mismatch; STOP.')
    print(f'BSE eigenvectors verified against sorted tables: {nfull} full / {nibz} IBZ points; '
          f'{len(exc.table)} transitions; max required exciton={requested}')
    records=[];maps=[]
    labels,vertices,hs_path=hs
    for group in groups:
        result=group_record(group,data,args.dark_threshold)
        weight,pairs=state_weight_and_transitions(exc,group,nfull)
        ibz_w=np.bincount(mapping,weights=weight,minlength=nibz)
        if not np.isclose(ibz_w.sum(),1.,rtol=0,atol=1e-8):
            raise ValueError('Full-zone to irreducible-zone weights do not conserve normalization')
        descriptor='_'.join(f'{s:04d}' for s in group)
        np.savez_compressed(out/f'weights_{descriptor}.npz',states=np.asarray(group,int),
                            IBZ_kred=np.asarray(lat_ibz.red_kpoints,float),IBZ_star_size=stars,
                            IBZ_star_weight=ibz_w,full_kred=np.asarray(lat.red_kpoints,float),
                            full_kweight=weight)
        write_csv(out/f'weights_{descriptor}_IBZ.csv',
                  ['IBZ_k_index_1based','k1','k2','k3','star_size','star_integrated_weight','weight_per_point_in_star'],
                  [[i+1,*xyz,int(stars[i]),float(ibz_w[i]),float(ibz_w[i]/stars[i])]
                   for i,xyz in enumerate(lat_ibz.red_kpoints)])
        write_csv(out/f'weights_{descriptor}_transitions.csv',
                  ['valence_band','conduction_band','mean_probability'],
                  [[v,c,p] for (v,c),p in sorted(pairs.items(),key=lambda x:-x[1])])
        plot_ibz_group(np.asarray(lat_ibz.red_kpoints),ibz_w,group,out,stars,labels,vertices)
        plot_first_bz_group(lat,weight,group,labels,vertices,out,args.marker_size)
        path=sampled_hs_path(lat,weight,labels,vertices,args.path_samples_per_edge)
        plot_sampled_path(path,group,out)
        write_csv(out/f'weights_{descriptor}_HS_path.csv',
                  ['path_distance_Ainv','k1','k2','k3','nearest_full_k_index_1based',
                   'distance_to_nearest_mesh_Ainv','weight_per_k'],
                  [[path['path_distance'][j],*p,int(path['sample_index'][j])+1,
                    path['sample_misfit_Ainv'][j],path['values'][j]]
                   for j,p in enumerate(path['path_red'])])
        weightsorted=sorted(pairs.items(),key=lambda x:-x[1])[:8]
        result.update(norm_full_k_weight=float(weight.sum()),norm_IBZ_star_weight=float(ibz_w.sum()),
                      participation_k=1/float(np.sum(weight**2)),
                      effective_num_IBZ_stars=1/float(np.sum(ibz_w**2)),
                      max_path_sampling_misfit_Ainv=float(np.max(path['sample_misfit_Ainv'])),
                      top_vc=[{'v':v,'c':c,'fraction':float(w)} for (v,c),w in weightsorted])
        records.append(result)
        maps.append({'states':list(group),'weak':result['weak_in_one_photon_y'],
                     'ibz':ibz_w,'ibz_red':np.asarray(lat_ibz.red_kpoints)})
        print(f'group {group}: E={result["mean_energy_eV"]:.6f} eV; '
              f'k_participation={result["participation_k"]:.1f} / {nfull}; '
              f'HS nearest-mesh max distance={result["max_path_sampling_misfit_Ainv"]:.3g} Å⁻¹')
    contrasted=plot_dark_bright_ibz(maps,out)
    return records,contrasted


def resolve_at_root(root,what):
    p=Path(what).expanduser()
    return p.resolve() if p.is_absolute() else (root/p).resolve()


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--root',default='.',help='Run from FS_TDSHF, or set --root explicitly')
    ap.add_argument('--out',default='ExcFLQ',help='Output directory inside --root')
    ap.add_argument('--hs',default='../../../IBZ/hBN.csv',help='HS path CSV relative to --root')
    ap.add_argument('--mode',choices=['catalog','bz','full'],default='full')
    ap.add_argument('--save',default='../SAVE',help='EQUILIBRIUM SAVE, not field-reduced FS_TDSHF/SAVE')
    ap.add_argument('--bse-dir',default='../BSE',help='Directory containing ndb.BS_diago_Q1')
    ap.add_argument('--bs-file',default='ndb.BS_diago_Q1')
    ap.add_argument('--energy-sorted',help='Path relative to --root, e.g. ../BSE_files/o-BSE.exc_qpt1_E_sorted')
    ap.add_argument('--intensity-sorted',help='Path relative to --root, e.g. ../BSE_files/o-BSE.exc_qpt1_I_sorted')
    ap.add_argument('--shg-spectrum',default='PyFLQ/HG_order_2.npz',
                    help='RT/Floquet SHG NPZ under root, or seven-column RT_order_2.dat')
    ap.add_argument('--no-shg',action='store_true',help='Do not read SHG file, catalogue-only resonance guide')
    ap.add_argument('--shg-source',choices=['auto','rt','fl'],default='auto')
    ap.add_argument('--component',choices=['x','y','z'],default='y')
    ap.add_argument('--groups',default='1,2;3,4;5;6;7,8;9;10,11',
                    help='1-based groups; or auto. Do NOT split near-degenerate pairs')
    ap.add_argument('--degeneracy-tol-eV',type=float,default=.0005)
    ap.add_argument('--max-groups',type=int,default=14)
    ap.add_argument('--max-states',type=int,default=30)
    ap.add_argument('--dark-threshold',type=float,default=1e-6,
                    help='Relative strength threshold: weak in ONE-PHOTON BSE at its computed polarization')
    ap.add_argument('--expected-kpoints',type=int,default=324,help='18x18x1 electronic BZ; set 0 to disable')
    ap.add_argument('--energy-match-tol-eV',type=float,default=.002)
    ap.add_argument('--peak-window-eV',type=float,default=.15)
    ap.add_argument('--shg-peak-prominence-rel',type=float,default=.02)
    ap.add_argument('--bright-reference-threshold',type=float,default=1e-3)
    ap.add_argument('--broadening-reference-eV',type=float,default=.20,
                    help='Reference broadening for a resolvability warning; not extracted from SHG')
    ap.add_argument('--dark-energy-min',type=float,default=5.5)
    ap.add_argument('--dark-energy-max',type=float,default=8.)
    ap.add_argument('--dark-zoom-min',type=float,default=2.75)
    ap.add_argument('--dark-zoom-max',type=float,default=3.65)
    ap.add_argument('--path-samples-per-edge',type=int,default=40)
    ap.add_argument('--marker-size',type=float,default=80)
    args=ap.parse_args(argv)
    if args.broadening_reference_eV<=0:ap.error('--broadening-reference-eV must be >0')
    if args.dark_threshold<=0 or args.dark_energy_min>=args.dark_energy_max:
        ap.error('Invalid dark threshold or excitonic energy window')
    if args.path_samples_per_edge<2:ap.error('--path-samples-per-edge must be at least 2')
    root=Path(args.root).expanduser().resolve()
    if not root.is_dir():raise FileNotFoundError(root)
    out=resolve_at_root(root,args.out)
    if out==root or root in out.parents and any(p.name.startswith('FLQ_') for p in out.parents):
        raise ValueError('Unsafe output location: output must be a dedicated ExcFLQ directory, not an FLQ input')
    out.mkdir(parents=True,exist_ok=True)
    efile=find_sorted(root,args.energy_sorted,'energy')
    ifile=find_sorted(root,args.intensity_sorted,'intensity')
    data=read_sorted(efile);compare_sorted(data,read_sorted(ifile))
    groups=parse_groups(args.groups,data,args.degeneracy_tol_eV,args.max_groups,args.max_states)
    shg=None;fl=None
    if not args.no_shg:
        spectrum=resolve_at_root(root,args.shg_spectrum)
        shg=select_shg(spectrum,args.component,args.shg_source)
        if spectrum.suffix.lower()=='.npz' and args.shg_source!='fl':
            try:fl=select_shg(spectrum,args.component,'fl')
            except ValueError:pass
    labels,vertices,hspath=parse_hs_csv(root,args.hs)
    print(f'Loading {len(data["energy"])} excitons from {efile}')
    print(f'HS path {hspath}: ' + ' → '.join(labels))
    print('Groups:',groups)
    print('No two-photon transition amplitude is inferred from one-photon BSE residuals.')
    write_csv(out/'exciton_catalogue.csv',
              ['state_1based','exciton_E_eV','normalized_one_photon_strength',
               'candidate_two_photon_pump_eV','candidate_one_photon_pump_eV',
               'weak_in_one_photon_for_BSE_polarization'],
              [[int(i),float(e),float(s),float(.5*e),float(e),int(s<args.dark_threshold)]
               for i,e,s in zip(data['index'],data['energy'],data['strength'])])
    records=[group_record(g,data,args.dark_threshold) for g in groups]
    write_csv(out/'multiplet_resonances.csv',
              ['states','E_mean_eV','splitting_eV','linear_strength_sum','linear_strength_max',
               'weak_in_one_photon','two_photon_pump_eV','one_photon_pump_eV'],
              [[r['label'],r['mean_energy_eV'],r['splitting_eV'],r['sum_linear_strength'],
                r['max_linear_strength'],int(r['weak_in_one_photon_y']),
                r['possible_2ph_pump_eV'],r['possible_1ph_pump_eV']] for r in records])
    plot_catalog(data,groups,out,args.dark_threshold,shg,show_count=35)
    dark=[];detected=[]
    if shg is not None:
        detected,dark=summarize_shg(shg,args,out,data,groups)
        plot_dark_spectrum(data,groups,shg,out,args,fl)
    bzdata=[];comparison=None
    if args.mode in ('bz','full'):
        bzdata,comparison=plot_bz_groups(root,args,data,groups,out,(labels,vertices,hspath))
    report={'script':'FLQ_Exc.py','output_dir':str(out),
            'electronic_data':'TD-HSEX BSE; DO NOT use independent IPA SAVE',
            'energy_sorted':str(efile),'intensity_sorted':str(ifile),
            'equilibrium_save':str(resolve_at_root(root,args.save)),
            'BSE_database':str(resolve_at_root(root,args.bse_dir)/args.bs_file),
            'high_symmetry_path':str(hspath),'shg_spectrum':None if shg is None else args.shg_spectrum,
            'shg_source':None if shg is None else shg[2],
            'exciton_count':len(data['index']),'groups':records,
            'k_resolved_records':bzdata,'dark_vs_bright_IBZ_states':comparison,
            'one_photon_dark_threshold':args.dark_threshold,
            'selected_dark_candidate_count':len(dark),'SHG_magnitude_local_peaks':len(detected),
            'interpretation':'E_exc and E_exc/2 are candidate SHG denominator positions. Linear BSE strength is NOT two-photon absorption, and k weights are NOT a decomposition of chi2. A quantitative dark-state SHG mechanism requires exciton-exciton dipoles and coherent sums. No anticrossing is inferred from a single calculated susceptibility spectrum.',
            'FLOrder':'UNCHANGED; no propagation or BSE computation'}
    (out/'analysis_summary.json').write_text(json.dumps(report,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n')
    print('Finished:',out)
    return 0

if __name__=='__main__':
    try:sys.exit(main())
    except Exception as exc:
        print(f'ERROR [{type(exc).__name__}]: {exc}',file=sys.stderr)
        raise
