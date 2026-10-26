#!/usr/bin/env python3
"""Extract direct-monodromy quasienergies and Fourier components once per run.

Input: SAVE/ns.db1 and FLQ_###/ndb.RT_V_bands{,_K_section}, ndb.Nonlinear.
Output: FLQ_###/PyFLQ/{quasienergies.npz,fourier.npz,FLQuasi.log}.
No polarization reconstruction or plotting occurs in this stage.

Convention: psi_a(t0+tau) = exp(-i epsilon_a tau) sum_h Lambda[a,h] exp(+i h omega tau).
The propagated initial states are psi(tau) @ inverse(Gamma), NOT the individual
Floquet eigenmodes. A projected occupied-subspace monodromy is not a complete
conduction+valence Floquet Hamiltonian; closure tests are mandatory.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import uuid
import numpy as np
import FLQ_AL as al


def diagonalize(C, i0, iT, i4, period):
    """Project U(T), diagonalize and independently validate at 4T."""
    nk,nb,nst=C.shape[1:]
    matrix_shape=(nk,nst,nst)
    q={name:np.empty(matrix_shape,complex) for name in ('L','L4','Gamma','Gram')}
    q.update({name:np.empty((nk,nst),complex) for name in ('multipliers','epsilon_principal_au')})
    checks={name:np.empty(nk) for name in ('gram_condition','eigenvector_condition',
        'closure_T','closure_4T','semigroup_4T','predicted_4T','offdiag_4T','gram_unitarity','gauge_T_error','gauge_4T_error')}
    for k in range(nk):
        C0,CT,C4=C[i0,k],C[iT,k],C[i4,k]
        G=C0.conj().T@C0
        checks['gram_condition'][k]=np.linalg.cond(G)
        if not np.isfinite(checks['gram_condition'][k]) or checks['gram_condition'][k]>1e12:
            raise al.CompareError(f'Ill-conditioned occupied Gram matrix at k={k+1}')
        L=np.linalg.solve(G,C0.conj().T@CT);L4=np.linalg.solve(G,C0.conj().T@C4)
        z,V=np.linalg.eig(L)
        if np.any(abs(z)<1e-300): raise al.CompareError('Zero monodromy eigenvalue')
        # Retain multiplier modulus rather than pretending a projected map is unitary.
        eps=-np.angle(z)/period+1j*np.log(abs(z))/period
        order=np.argsort(eps.real);z=z[order];eps=eps[order];V=V[:,order]
        for a in range(nst):
            norm=np.vdot(V[:,a],G@V[:,a]).real
            if norm<=0: raise al.CompareError('Nonpositive eigenvector metric norm')
            V[:,a]/=np.sqrt(norm)
        checks['eigenvector_condition'][k]=np.linalg.cond(V)
        if not np.isfinite(checks['eigenvector_condition'][k]) or checks['eigenvector_condition'][k]>1e12:
            raise al.CompareError(f'Defective/ill-conditioned Floquet eigenvectors at k={k+1}')
        Lpow=np.linalg.matrix_power(L,4);rot=np.linalg.solve(V,L4@V)
        checks['closure_T'][k]=al.relnorm(C0@L-CT,CT)
        checks['closure_4T'][k]=al.relnorm(C0@L4-C4,C4)
        checks['semigroup_4T'][k]=al.relnorm(L4-Lpow,L4)
        checks['predicted_4T'][k]=al.relnorm(C0@Lpow-C4,C4)
        checks['offdiag_4T'][k]=al.relnorm(rot-np.diag(np.diag(rot)),np.diag(np.diag(rot)))
        checks['gram_unitarity'][k]=al.relnorm(L.conj().T@G@L-G,G)
        shifts=np.arange(-4,5)[:,None]*2*np.pi/period
        checks['gauge_T_error'][k]=float(np.max(abs(np.exp(-1j*(eps[None,:]+shifts)*period)-z[None,:])))
        checks['gauge_4T_error'][k]=float(np.max(abs(np.exp(-4j*(eps[None,:]+shifts)*period)-z[None,:]**4)))
        for name,value in [('L',L),('L4',L4),('Gamma',V),('Gram',G),('multipliers',z),('epsilon_principal_au',eps)]:
            q[name][k]=value
    q.update(checks)
    q['xi_principal_eV']=q['epsilon_principal_au'].real*al.HA_TO_EV
    q['projected_imaginary_energy_eV']=q['epsilon_principal_au'].imag*al.HA_TO_EV
    return q


def fourier_modes(Cfit, tau, Gamma, eps, omega, nf):
    """Exact sampled DFT with compact integer gauge; no invented extra harmonics."""
    nk,nb,nst=Cfit.shape[1:];harm=np.arange(-nf,nf+1);ns=len(harm)
    if len(tau)!=ns: raise al.CompareError('Fourier sample count must be 2*NF+1')
    Bminus=np.exp(-1j*np.outer(tau,harm*omega))
    modes=np.einsum('tkbi,kia->tkba',Cfit,Gamma,optimize=True)
    Lam=np.empty((nk,nst,ns,nb),complex);used=np.empty_like(eps)
    shifts=np.empty((nk,nst),int);edge=np.empty((nk,nst));moment=np.empty((nk,nst))
    for k in range(nk):
        for a in range(nst):
            candidates=[]
            for m in range(-nf,nf+1):
                e=eps[k,a]+m*omega
                D=modes[:,k,:,a]*np.exp(1j*e*tau)[:,None]
                coef=Bminus.T@D/ns
                weights=np.sum(abs(coef)**2,axis=1);total=max(float(weights.sum()),1e-300)
                ed=float((weights[0]+weights[-1])/total)
                mo=float(np.sum((harm/max(nf,1))**2*weights)/total)
                candidates.append((ed,mo,abs(m),m,e,coef))
            candidates.sort(key=lambda x:x[:3]);ed,mo,_,m,e,coef=candidates[0]
            Lam[k,a]=coef;used[k,a]=e;shifts[k,a]=m;edge[k,a]=ed;moment[k,a]=mo
    return {'Lambda':Lam,'epsilon_compact_au':used,'gauge_shift_m':shifts,
            'harmonics':harm,'edge_fraction':edge,'spectral_moment':moment}


def reconstruct_modes(f, tau):
    """Return physical Floquet eigenmodes [time,k,basis,mode]."""
    tau=np.atleast_1d(np.asarray(tau,float));omega=float(f['omega_au'])
    B=np.exp(1j*np.outer(tau,f['harmonics']*omega))
    D=np.einsum('th,kahb->tkba',B,f['Lambda'],optimize=True)
    return D*np.exp(-1j*tau[:,None,None]*f['epsilon_compact_au'][None,:,:])[:,:,None,:]


def reconstruct_states(q, f, tau):
    """Recover propagated occupied orbitals, including the initial-state mixing."""
    modes=reconstruct_modes(f,tau)
    return np.einsum('tkba,kaj->tkbj',modes,np.linalg.inv(q['Gamma']),optimize=True)


def mode_band_assignment(Gamma):
    from scipy.optimize import linear_sum_assignment
    mapping=np.empty(Gamma.shape[:2],int);confidence=np.empty(Gamma.shape[:2])
    for k,g in enumerate(Gamma):
        weights=abs(g)**2;weights/=np.maximum(weights.sum(axis=0,keepdims=True),1e-300)
        rows,cols=linear_sum_assignment(-weights)
        mapping[k,rows]=cols;confidence[k,rows]=weights[rows,cols]
    return mapping,confidence


def extract_one(calc, lattice, save_path, args):
    dest=calc/args.out_name;dest.mkdir(parents=True,exist_ok=True)
    outputs=[dest/'quasienergies.npz',dest/'fourier.npz']
    if any(p.exists() for p in outputs) and not args.force:
        raise al.CompareError(f'{dest}: extraction exists; use --force to replace only derived files')
    # Remove only old derived outputs on explicit --force: a failed replacement cannot look current.
    if args.force:
        for p in outputs: p.unlink(missing_ok=True)
    with al.run_log(dest/'FLQuasi.log'):
        print(f'Calculation: {calc.name}; SAVE: {save_path}')
        t,rawC,raw,omega,info=al.load_vbands(calc);T=2*np.pi/omega;nf=info['nf']
        i0,fit,iT,i4,timing=al.exact_samples(t,raw,T,nf,args.time_tol)
        if args.expect_bands and list(args.expect_bands)!=[info['band_lo'],info['band_hi']]:
            raise al.CompareError(f'Actual NLBands {info["band_lo"]}:{info["band_hi"]} != requested expectation')
        if len(lattice['kpoints_red'])!=rawC.shape[1]: raise al.CompareError('SAVE/V_bands k-point counts differ')
        from yambopy import YamboNLDB
        nl=YamboNLDB(folder=str(calc.parent),calc=calc.name)
        if len(nl.Efield)!=1:raise al.CompareError('Quasienergy extraction needs one driving frequency per folder')
        for other in getattr(nl,'Efield_general',[])[1:]:
            if str(other.get('name','none')).lower()!='none' and float(other.get('intensity',0))>0:
                raise al.CompareError('Multiple active fields are outside this single-frequency extractor')
        field=nl.Efield[0];direction=np.asarray(field['versor'],float)
        if direction.shape!=(3,) or not np.isfinite(direction).all() or np.linalg.norm(direction)<=0:
            raise al.CompareError('Invalid field direction')
        direction=direction/np.linalg.norm(direction)
        for isym in np.unique(lattice['symmetry_indexes']):
            if lattice['symmetry_antiunitary'][isym]:
                raise al.CompareError('Antiunitary BZ expansion is not supported for this driven state; use field-fixed SAVE')
            if not np.allclose(lattice['symmetry_cart'][isym]@direction,direction,atol=2e-6,rtol=0):
                raise al.CompareError('SAVE expansion uses a symmetry broken by the applied field; regenerate FixSymm')
        nocc=int(lattice['nocc']);lo=info['band_lo'];hi=info['band_hi']
        if not (1<=lo<=nocc<hi): raise al.CompareError('Active window must include occupied and empty states')
        norms=np.linalg.norm(rawC[i0],axis=1);mx=np.max(norms,axis=0)
        active=np.flatnonzero(mx>args.active_rtol*max(float(mx.max()),1.0))
        expected=np.arange(lo-1,nocc)
        if np.array_equal(active,expected): slot_layout='absolute_occupied_band_index'
        elif rawC.shape[-1]==len(expected) and np.array_equal(active,np.arange(len(expected))):
            slot_layout='local_dynamic_occupied_index'
        else: raise al.CompareError(f'Active raw slots {active.tolist()} inconsistent with NLBands and nocc')
        C=rawC[:,:,:,active]
        q=diagonalize(C,i0,iT,i4,T)
        for name in ('gram_condition','eigenvector_condition','closure_T','closure_4T','semigroup_4T','predicted_4T','offdiag_4T','gram_unitarity','gauge_T_error','gauge_4T_error'):
            print(f'{name:28s} min/mean/max = {q[name].min():.6e} {q[name].mean():.6e} {q[name].max():.6e}')
        print('Exact timing errors / T:',timing)
        print('Max ||multiplier|-1|:',float(np.max(abs(abs(q['multipliers'])-1))))
        print('Imaginary projected energies are diagnostics, NOT measured lifetimes.')
        for key in ('closure_T','closure_4T','semigroup_4T','predicted_4T'):
            if q[key].max()>args.reconstruction_tol:
                raise al.CompareError(f'{key} exceeds {args.reconstruction_tol:g}; no valid Floquet cache written')
        f=fourier_modes(C[fit],t[fit]-t[i0],q['Gamma'],q['epsilon_principal_au'],omega,nf)
        f.update(omega_au=omega,period_au=T,t0_au=float(t[i0]))
        idx=np.flatnonzero(t>=t[i0])
        rec=reconstruct_states(q,f,t[idx]-t[i0])
        errors=np.array([al.relnorm(a-b,b) for a,b in zip(rec,C[idx])])
        al.need_finite('Fourier reconstruction errors',errors)
        print('Physical-state reconstruction errors at stored post-t0 records:',errors)
        if errors.max()>args.reconstruction_tol: raise al.CompareError('Fourier state reconstruction failed')
        mapping,confidence=mode_band_assignment(q['Gamma'])
        print('Compact integer gauges:',np.unique(f['gauge_shift_m'],return_counts=True))
        print('Maximum Fourier edge fraction:',float(f['edge_fraction'].max()))
        print('WARNING: small edge weights and interpolation-node agreement do not prove harmonic convergence.')
        print(f'Only {len(active)} represented occupied modes extracted; no unpropagated conduction modes invented.')
        print('Frozen occupied bands:',list(range(1,lo)))
        files=[save_path/'ns.db1',calc/'ndb.RT_V_bands',calc/'ndb.RT_V_bands_K_section',calc/'ndb.Nonlinear']
        settings=al.run_settings(calc)
        if settings['restart_complete']:
            print('WARNING: latest report is a completed-run restart. HG will require a valid dense polarization history.')
        meta={'schema':al.SCHEMA,'passed':True,'calc':calc.name,'photon_energy_eV':omega*al.HA_TO_EV,
              'sources':al.source_stamp(files),'settings':settings,'band_window':[lo,hi],
              'active_bands':list(range(lo,nocc+1)),'raw_slot_layout':slot_layout,
              'field_direction':direction.tolist(),'field_intensity':float(field['intensity']),
              'BZ_expansion':'unitary field-preserving symmetries from the actual driven SAVE',
              'timing_errors':timing,'reconstruction_error_max':float(errors.max()),'reader':info,
              'scope':'projected represented occupied subspace; conditional Floquet reconstruction',
              'fourier_convention':'psi=exp(-i epsilon tau)*sum_eta Lambda_eta exp(+i eta omega tau)',
              'ks_energy_note':'KS energies are read unchanged from SAVE; no GW/scissor is inferred or applied'}
        generation=str(uuid.uuid4());common={'schema':al.SCHEMA,'generation':generation,
            'metadata_json':json.dumps(al.jsafe(meta),allow_nan=False),'omega_au':omega,'period_au':T,
            't0_au':float(t[i0]),'nf':nf,'active_bands':np.arange(lo,nocc+1),
            'basis_bands':np.arange(lo,hi+1),'active_raw_slots':active,
            'field_direction':direction,'field_intensity':float(field['intensity'])}
        q.update(lattice);q.update(common)
        q.update(mode_for_band=mapping,mode_assignment_weight=confidence)
        f.update(common);f.update(raw_times_au=t,raw_indices=raw,C_raw=C,
            phase_sample_indices=fit,T_index=iT,fourT_index=i4,
            stored_reconstruction_indices=idx,stored_reconstruction_errors=errors)
        al.atomic_npz(outputs[0],**q);al.atomic_npz(outputs[1],**f)
        print('Written:',*outputs,sep='\n  ')
        return meta


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workdir',default='.');p.add_argument('--save',default='SAVE')
    p.add_argument('--calc-glob',default='FLQ_[0-9][0-9][0-9]');p.add_argument('--calc')
    p.add_argument('--out-name',default='PyFLQ');p.add_argument('--force',action='store_true')
    p.add_argument('--time-tol',type=float,default=1e-8)
    p.add_argument('--reconstruction-tol',type=float,default=1e-6)
    p.add_argument('--active-rtol',type=float,default=1e-10)
    p.add_argument('--expect-bands',nargs=2,type=int);p.add_argument('--expected-count',type=int)
    args=p.parse_args(argv);work=Path(args.workdir).expanduser().resolve()
    save=work/args.save;calcs=al.calculation_dirs(work,args.calc or args.calc_glob)
    if args.expected_count is not None and len(calcs)!=args.expected_count:
        raise al.CompareError(f'Expected {args.expected_count} frequency folders; found {len(calcs)}')
    with al.run_log(work/args.out_name/'FLQuasi.log'):
        lattice=al.load_lattice(save)
        print('Stored/full BZ k-points:',len(lattice['kpoints_red']),len(lattice['kpoints_full_red']))
        print('Full grid:',lattice['grid_shape'])
        for calc in calcs: extract_one(calc,lattice,save,args)
        print(f'Completed {len(calcs)} extractions')

if __name__=='__main__': main()
