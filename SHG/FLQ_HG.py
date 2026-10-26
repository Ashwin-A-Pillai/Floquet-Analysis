#!/usr/bin/env python3
"""RT and Floquet harmonic generation with order-aware normalization.

--harmonic-order 2 selects SHG; 3 selects THG. --fit-order is a distinct
polarization-fit cutoff (default max(4,2*harmonic_order)). Stored FLOrder is
read from the databases, never increased by postprocessing. One monochromatic
SIN/SOFTSIN field only: this is not a two-frequency SFG/DFG implementation.

A finite-field Fourier coefficient divided by E**n is an effective contracted
response. An intensity-convergence study is needed to call it perturbative chi(n).
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import FLQ_AL as al
from FLQuasi import reconstruct_states


def check_single_field(db):
    for field in getattr(db,'Efield_general',[])[1:]:
        if str(field.get('name','none')).strip().lower()!='none' and float(field.get('intensity',0))>0:
            raise al.CompareError('Multiple active fields require multi-frequency response analysis')
    for field in db.Efield:
        if str(field.get('name','')).strip().upper() not in ('SIN','SOFTSIN'):
            raise al.CompareError('Supported driving fields are SIN and SOFTSIN')
        if float(field.get('intensity',0))<=0:raise al.CompareError('Nonpositive field intensity')
        e=np.asarray(field['versor'],float)
        if e.shape!=(3,) or not np.isfinite(e).all() or np.linalg.norm(e)==0:
            raise al.CompareError('Invalid field direction')


def rt_response(work,args):
    from yambopy import YamboNLDB
    from yambopy.nl.sin_analysis import Xn_from_sine
    from yambopy.nl.external_efield import Divide_by_the_Field
    db=YamboNLDB(folder=str(work),calc=args.rt_calc);check_single_field(db)
    times=np.asarray(db.IO_TIME_points,float)
    al.need_finite('RT times',times)
    if np.any(np.diff(times)<=0):raise al.CompareError('RT polarization timestamps are not increasing')
    for i,P in enumerate(db.Polarization):
        P=np.asarray(P,float);al.need_finite(f'RT polarization {i}',P)
        if P.shape!=(3,len(times)):raise al.CompareError('Unsupported RT polarization layout')
        norm=np.linalg.norm(P-P[:,0,None],axis=0);mx=norm.max()
        if mx and len(times)>20 and np.count_nonzero(norm>mx*1e-12)<=2:
            raise al.CompareError(f'RT frequency {i} has a sparse/restart-corrupted history')
    # Fresh mutable T_range per analysis avoids shared-state dependence in older YamboPy.
    ana=Xn_from_sine(db,X_order=args.fit_order,T_range=[-1.,-1.])
    raw=ana.perform_analysis();order=args.harmonic_order
    chi=np.array([raw[order,i]*Divide_by_the_Field(db.Efield[i],order)*ana.get_Unit_of_Measure(order)
                  for i in range(len(db.Polarization))],complex)
    E=np.asarray(ana.freqs,float)*al.HA_TO_EV
    if len(E)>1 and np.any(np.diff(E)<=0):raise al.CompareError('RT frequencies are not strictly increasing')
    print(f'RT spectrum: {len(E)} frequencies; YamboPy Xn_from_sine, fit order {args.fit_order}')
    return E,chi,raw,db


def berry_validation(P,reference,tol,floor):
    dP=P-P[0];dR=reference-reference[0];refnorm=np.linalg.norm(dR,axis=0)
    largest=float(refnorm.max())
    if largest<=0: raise al.CompareError('Stored polarization has no induced signal to validate')
    active=refnorm>floor*largest
    errors=np.linalg.norm(dP-dR,axis=0)/np.maximum(refnorm,1e-300)
    combined=al.relnorm(dP[:,active]-dR[:,active],dR[:,active])
    leakage=float(np.linalg.norm((dP-dR)[:,~active])/max(np.linalg.norm(dR[:,active]),1e-300))
    metrics={'active_components':[c for c,yes in zip('xyz',active) if yes],
        'relative_errors_per_component':[float(v) if a else None for v,a in zip(errors,active)],
        'combined_active_error':combined,'inactive_leakage':leakage}
    if combined>tol or np.any(errors[active]>tol):
        raise al.CompareError(f'Raw Berry validation failed: {metrics}; tolerance {tol:g}')
    return metrics


def sampling_change(coarse, fine, projector=None):
    """Convergence of the symmetry-allowed *induced* harmonic response.

    Tiny forbidden components can have large relative errors despite being
    numerically irrelevant (e.g., DC offset leaking into chi_x for E along y).
    Reject those directions using the field-preserving polar-vector projector.
    The physical relative tolerance and the original amplitude floor are unchanged.
    """
    coarse = np.asarray(coarse, complex)
    fine = np.asarray(fine, complex)
    if coarse.shape != fine.shape or coarse.shape != (3,):
        raise al.CompareError('Expected two x/y/z harmonic response vectors')
    if projector is None:
        allowed = np.ones(3, dtype=bool)
    else:
        projector = np.asarray(projector, float)
        if projector.shape != (3, 3):
            raise al.CompareError('Expected 3x3 polar-vector projector')
        # Stored symmetry matrices have finite numeric precision. A row of the
        # projector smaller than 1e-6 represents a forbidden output direction.
        allowed = np.linalg.norm(projector, axis=1) > 1e-6
    if not allowed.any():
        raise al.CompareError('No symmetry-allowed output direction for response')
    scale = max(float(np.max(abs(fine[allowed]))), 1e-300)
    active = allowed & (abs(fine) > 1e-8 * scale)
    if active.any():
        return float(np.max(abs(coarse[active]-fine[active]) /
                            np.maximum(abs(fine[active]), scale*1e-8)))
    # A zero response is converged only when consecutive evaluations also agree.
    return float(np.max(abs(coarse[allowed]-fine[allowed])) / scale)


def floquet_response(calc,work,args,rt_db=None):
    from yambopy import YamboNLDB
    from yambopy.nl.sin_analysis import Xn_from_sine
    q,f,meta=al.load_pair(calc,args.out_name,verify_sources=True)
    nf=int(q['nf']);H=args.harmonic_order
    if H>nf:
        raise al.CompareError(f'Requested harmonic {H} exceeds stored FLOrder={nf}; saved data cannot support this harmonic')
    db=YamboNLDB(folder=str(work),calc=calc.name);check_single_field(db)
    if len(db.Polarization)!=1:raise al.CompareError('Each FLQ directory must contain one frequency')
    field=db.Efield[0];energy=float(q['omega_au'])*al.HA_TO_EV
    if not np.isclose(float(field['freq_range'][0])*al.HA_TO_EV,energy,atol=1e-5,rtol=0):
        raise al.CompareError('Cached frequency and ndb.Nonlinear disagree')
    if rt_db is not None:
        rtenergy=np.array([float(x['freq_range'][0])*al.HA_TO_EV for x in rt_db.Efield])
        i=int(np.argmin(abs(rtenergy-energy)))
        if abs(rtenergy[i]-energy)>args.energy_tol:
            if args.rt_reference!='pchip':
                raise al.CompareError('Floquet frequency has no matching RT frequency; '
                                      'off-grid comparison requires --rt-reference pchip')
            if energy<rtenergy.min() or energy>rtenergy.max():
                raise al.CompareError('Floquet energy outside RT reference frequency range')
        ref=rt_db.Efield[i]
        e1=np.asarray(ref['versor'],float);e2=np.asarray(field['versor'],float)
        if not np.allclose(e1/np.linalg.norm(e1),e2/np.linalg.norm(e2),atol=1e-7,rtol=0):
            raise al.CompareError('RT/Floquet field directions differ')
        if not np.isclose(float(ref['intensity']),float(field['intensity']),rtol=1e-6):
            raise al.CompareError('RT/Floquet intensities differ')
        rtset=al.run_settings(work/args.rt_calc)
        for key in ('correlation','integrator','qp_database','damping','qp_parameters','band_window'):
            a,b=rtset.get(key),meta['settings'].get(key)
            if a is not None and b is not None and (a.strip()!=b.strip() if isinstance(a,str) and isinstance(b,str) else a!=b):
                raise al.CompareError(f'RT/Floquet {key} differ: {a} vs {b}')
    mapping=q['bz_to_ibz'].astype(int);shape,axes,_=al.infer_regular_grid(q['kpoints_full_red'],args.k_match_tol)
    lo=int(q['basis_bands'][0]);hi=int(q['basis_bands'][-1]);nocc=int(q['nocc'])
    DIP,overlap_meta=al.load_yambo_DIP_S(calc,len(mapping),hi)
    def berry(C):
        return al.berry_polarization(C[:,mapping],DIP,axes,q['lattice_bohr'],float(q['volume_bohr3']),
            float(q['spin_factor']),nocc,lo,hi)
    Praw=berry(f['C_raw']);reference=al.load_polarization_at(db,f['raw_times_au'])
    sym=args.symmetry
    if sym=='auto':
        sym='residual' if meta['settings'].get('symmetrize')=='yes' else 'none'
        if meta['settings'].get('symmetrize') is None:
            print('WARNING: report does not specify symmetrization; auto uses none. Specify --symmetry only with evidence.')
    projector=al.polar_projector(q,field['versor'],sym)
    P=al.project_delta(Praw,projector)
    dest=calc/args.out_name;debug=dest/f'HG_order_{H}_Berry_debug.npz'
    al.atomic_npz(debug,times_au=f['raw_times_au'],raw_unprojected=Praw,raw_projected=P,
        stored=reference,projector=projector)
    print(f'{calc.name}: photon={energy:.8f} eV, Berry grid={shape}, symmetry={sym}')
    print('Unprojected Berry error:',al.berry_difference_metrics(Praw,reference)[0])
    print('Applied polar-vector projector:\n',projector)
    allowed_axes = np.linalg.norm(projector, axis=1) > 1e-6
    print('Quadrature convergence tests symmetry-allowed outputs:',
          [axis for axis, active in zip('xyz', allowed_axes) if active])
    # The projector must be consistent with what Lumen actually stored.
    stored_projection=al.project_delta(reference,projector)
    if al.relnorm(stored_projection-reference,reference-reference[0])>args.berry_tol:
        raise al.CompareError('Stored polarization is not invariant under the requested projector')
    metrics=berry_validation(P,reference,args.berry_tol,args.component_floor)
    print('Raw Berry checks:',json.dumps(al.jsafe(metrics)))
    if metrics['inactive_leakage']>args.leakage_tol:print('WARNING: transverse/inactive Berry leakage exceeds diagnostic threshold')
    nstate=len(q['active_bands']);n=max(args.polarization_points,4*nf*nstate+1,2*args.fit_order+1)
    if n%2==0:n+=1
    ana=Xn_from_sine(db,X_order=args.fit_order,T_range=[-1.,-1.])
    def sample(number):
        tau=np.arange(number)*float(q['period_au'])/number
        states=reconstruct_states(q,f,tau)
        pol=al.project_delta(berry(states),projector)
        # Fit the induced polarization, not the large static Berry offset.
        # A constant cannot contribute at n>0; subtracting it avoids roundoff
        # leakage into a forbidden harmonic component. Keep 'pol' intact for output.
        induced=pol-pol[0]
        coeff=al.harmonic_fit(float(q['t0_au'])+tau,induced,float(q['omega_au']),args.fit_order,float(field['initial_time']))
        chi=np.asarray(al.normalize_one_harmonic(coeff,db,ana,H),complex)
        return tau,pol,coeff,chi
    tau,pol,coeff,chi=sample(n);changes=[];converged=False
    for _ in range(args.max_refinements):
        newn=2*n+1;tau2,pol2,coeff2,chi2=sample(newn)
        change=sampling_change(chi,chi2,projector=projector);changes.append(change)
        print(f'Polarization quadrature {n}->{newn}: harmonic-{H} relative change {change:.6e}')
        n=newn;tau,pol,coeff,chi=tau2,pol2,coeff2,chi2
        if change<=args.sampling_tol:converged=True;break
    if not converged:
        raise al.CompareError('Requested harmonic did not converge with denser reconstructed polarization sampling')
    result={'calc':calc.name,'energy_eV':energy,'harmonic_order':H,'fit_order':args.fit_order,
        'field_versor':np.asarray(field['versor'],float),
        'stored_FLOrder':nf,'symmetry':sym,'berry':metrics,'quadrature_points':n,'quadrature_changes':changes,
        'overlap':overlap_meta,'extraction_generation':str(q['generation']),
        'note':'Convergence checked on symmetry-allowed induced harmonics only; '
               'this does not establish FLOrder, k-grid, band or intensity convergence'}
    al.atomic_npz(dest/f'HG_order_{H}.npz',schema=al.SCHEMA,generation=q['generation'],
        energy_eV=energy,harmonic_order=H,fit_order=args.fit_order,chi=chi,raw_harmonics=coeff,
        times_au=float(q['t0_au'])+tau,polarization=pol,projector=projector,
        metadata_json=json.dumps(al.jsafe(result),allow_nan=False))
    return energy,chi,result,q


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workdir',default='.');p.add_argument('--save',default='SAVE')
    p.add_argument('--out-name',default='PyFLQ');p.add_argument('--rt-calc',default='NL')
    p.add_argument('--calc-glob',default='FLQ_[0-9][0-9][0-9]')
    p.add_argument('--source',choices=['rt','floquet','both'],default='both')
    p.add_argument('--harmonic-order',type=int,default=2);p.add_argument('--fit-order',type=int)
    p.add_argument('--component',default='y',
                   help='x/y/z = field-contracted output; "[1,1,2]" = Cartesian chi_xxy (SHG).')
    p.add_argument('--tensor-symmetry-db',
                   help='Full equilibrium SAVE/ns.db1 for symmetry-derived tensor extraction; '
                        'defaults to ../SAVE/ns.db1, never the field-reduced SAVE')
    p.add_argument('--expected-count',type=int,help='Legacy: expect this count in every requested source')
    p.add_argument('--rt-reference',choices=['exact','pchip'],default='exact',
                   help='Exact RT frequency matching (default), or explicit PCHIP interpolation for off-grid plans')
    p.add_argument('--expected-rt-count',type=int)
    p.add_argument('--expected-fl-count',type=int)
    p.add_argument('--force',action='store_true')
    p.add_argument('--symmetry',choices=['auto','none','residual','xy'],default='auto')
    p.add_argument('--berry-tol',type=float,default=1e-5);p.add_argument('--component-floor',type=float,default=1e-10)
    p.add_argument('--leakage-tol',type=float,default=2e-3);p.add_argument('--energy-tol',type=float,default=1e-4)
    p.add_argument('--k-match-tol',type=float,default=2e-5)
    p.add_argument('--polarization-points',type=int,default=65)
    p.add_argument('--sampling-tol',type=float,default=5e-4);p.add_argument('--max-refinements',type=int,default=3)
    p.add_argument('--effective-thickness-ang',type=float,
        help='Optional 2D effective-bulk rescaling; default is raw supercell convention')
    args=p.parse_args(argv)
    if args.harmonic_order<1:p.error('--harmonic-order must be >=1')
    al.parse_tensor_component(args.component,args.harmonic_order)
    if args.fit_order is None:args.fit_order=max(4,2*args.harmonic_order)
    if args.fit_order<args.harmonic_order:p.error('--fit-order must include the requested harmonic')
    if args.max_refinements<1:p.error('At least one sampling refinement is required')
    work=Path(args.workdir).expanduser().resolve();dest=work/args.out_name;dest.mkdir(parents=True,exist_ok=True)
    H=args.harmonic_order;outfile=dest/f'HG_order_{H}.npz'
    if outfile.exists() and not args.force:
        raise al.CompareError(f'{outfile} exists; use --force to replace derived harmonic results')
    with al.run_log(dest/f'FLQ_HG_order_{H}.log'):
        print(f'Source={args.source}; harmonic={H}; fit-order={args.fit_order}')
        payload={'schema':al.SCHEMA,'harmonic_order':H,'fit_order':args.fit_order}
        rt_db=None;validations=[];lattice=None;field_vec=None;catalog=None
        if args.source in ('rt','both'):
            Er,Xr,raw,rt_db=rt_response(work,args)
            field_vec=np.asarray(rt_db.Efield[0]['versor'],float)
            for ef in rt_db.Efield[1:]:
                v=np.asarray(ef['versor'],float)
                if not np.allclose(v/np.linalg.norm(v),field_vec/np.linalg.norm(field_vec),rtol=0,atol=1e-7):
                    raise al.CompareError('RT source contains different field directions; '
                                          'a single tensor contraction cannot describe its spectrum')
        if args.source=='floquet':
            # Detect tensor availability BEFORE expensive Floquet/Berry processing.
            from yambopy import YamboNLDB
            first=al.calculation_dirs(work,args.calc_glob)[0]
            first_db=YamboNLDB(folder=str(work),calc=first.name)
            check_single_field(first_db)
            field_vec=np.asarray(first_db.Efield[0]['versor'],float)
        catalog=al.tensor_catalog_for_work(work,args.save,H,field_vec,args.tensor_symmetry_db)
        al.describe_tensor_catalog(catalog,args.component)
        # Refuse to pretend that an unconstrained single-field measurement is an
        # individual tensor entry. The full x/y/z contracted spectrum is still saved.
        al.tensor_response(catalog,args.component,np.zeros((1,3),complex))
        if args.source in ('rt','both'):
            expected_rt=args.expected_rt_count if args.expected_rt_count is not None else args.expected_count
            if expected_rt is not None and len(Er)!=expected_rt:
                raise al.CompareError(f'Expected {expected_rt} RT frequencies, got {len(Er)}')
            payload.update(rt_energy_eV=Er,rt_chi=Xr,rt_raw_harmonics=raw)
        if args.source in ('floquet','both'):
            calcs=al.calculation_dirs(work,args.calc_glob)
            expected_fl=args.expected_fl_count if args.expected_fl_count is not None else args.expected_count
            if expected_fl is not None and len(calcs)!=expected_fl:
                raise al.CompareError(f'Expected {expected_fl} Floquet folders, got {len(calcs)}')
            rows=[]
            for calc in calcs:
                with al.run_log(calc/args.out_name/f'FLQ_HG_order_{H}.log'):
                    e,x,validation,q=floquet_response(calc,work,args,rt_db)
                v=np.asarray(validation['field_versor'],float)
                if not np.allclose(v/np.linalg.norm(v),field_vec/np.linalg.norm(field_vec),rtol=0,atol=1e-7):
                    raise al.CompareError(f'{calc.name}: field direction differs across the Floquet grid')
                rows.append((e,x));validations.append(validation);lattice=q
            rows.sort(key=lambda r:r[0]);Ef=np.array([r[0] for r in rows]);Xf=np.array([r[1] for r in rows])
            if len(Ef)>1 and np.any(np.diff(Ef)<=0):raise al.CompareError('Duplicate Floquet frequencies')
            payload.update(fl_energy_eV=Ef,fl_chi=Xf)
        factor=1.;convention='supercell volume; no 2D thickness rescaling'
        if args.effective_thickness_ang is not None:
            if args.effective_thickness_ang<=0:raise al.CompareError('Thickness must be positive')
            if lattice is None:lattice=al.load_lattice(work/args.save)
            axes=lattice['varying_axes'].astype(int)
            if len(axes)!=2:raise al.CompareError('Effective-thickness normalization applies only to a 2D mesh')
            A=lattice['lattice_bohr'];area=np.linalg.norm(np.cross(A[axes[0]],A[axes[1]]))
            height=float(lattice['volume_bohr3'])/area*al.BOHR_TO_ANG
            factor=height/args.effective_thickness_ang
            convention=f'effective bulk: Lperp={height:.10g} Å / deff={args.effective_thickness_ang:.10g} Å'
        units='dimensionless' if H==1 else ('cm/statV' if H==2 else f'(cm/statV)^{H-1}')
        for name in ('rt','fl'):
            if name+'_chi' in payload:payload[name+'_chi']*=factor
        metrics={}
        if args.source=='both':
            # Adaptive Floquet grids are subsets of the denser RT grid.
            # Match actual RT energies one by one: no spectral interpolation.
            rt_match, rt_idx, max_deviation=al.reference_on_floquet_grid(
                Er,payload['rt_chi'],Ef,args.energy_tol,method=args.rt_reference)
            payload.update(rt_matched_energy_eV=Ef,rt_matched_chi=rt_match,
                           rt_matching_indices=rt_idx,rt_reference_method=args.rt_reference)
            for i,comp in enumerate('xyz'):
                metrics[comp]=al.comparison_metrics(Ef,rt_match,Ef,payload['fl_chi'],i,args.energy_tol)
                metrics[comp].update(reference_matching='exact_RT_subset' if args.rt_reference=='exact' else 'PCHIP_interpolated_complex_RT',
                                     reference_points=len(Er),floquet_points=len(Ef),
                                     max_matching_deviation_eV=max_deviation)
        # One named tensor element is a symmetry-constrained linear combination
        # of the x/y/z contractions, not a fourth independent measured output.
        selected_meta=None
        for name in ('rt','fl'):
            if name+'_chi' in payload:
                vals,selected_meta=al.tensor_response(catalog,args.component,payload[name+'_chi'])
                payload[name+'_selected_chi']=vals
        payload.update(selected_component=args.component,
                       selected_label=selected_meta['label'],
                       selected_kind=selected_meta['status'],
                       tensor_catalog_json=json.dumps(al.jsafe(catalog),allow_nan=False))
        if args.source=='both':
            rt_sel,_=al.tensor_response(catalog,args.component,payload['rt_matched_chi'])
            fl_sel=payload['fl_selected_chi']
            metrics['selected']=al.comparison_metrics(Ef,rt_sel[:,None],Ef,fl_sel[:,None],0,args.energy_tol)
            metrics['selected']['label']=selected_meta['label']
        payload.update(units=units,thickness_factor=factor,normalization=convention,
            metadata_json=json.dumps(al.jsafe({'arguments':vars(args),'validation':validations,
                'comparison':metrics,'quantity':'field-contracted effective harmonic susceptibility',
                'normalization':'YamboPy Divide_by_the_Field and get_Unit_of_Measure',
                'note':'Not a complete Cartesian tensor; no finite-field perturbative extrapolation performed'}),allow_nan=False))
        al.atomic_npz(outfile,**payload)
        for name,label in [('rt','RT'),('fl','Floquet')]:
            if name+'_chi' in payload:
                path=dest/f'{label}_order_{H}.dat'
                al.atomic_text(path,al.spectrum_text(payload[name+'_energy_eV'],payload[name+'_chi'],
                    f'Harmonic {H}; {units}; {convention}'))
        for name,label in [('rt','RT'),('fl','Floquet')]:
            if name+'_selected_chi' in payload:
                energies=payload[name+'_energy_eV'];values=payload[name+'_selected_chi']
                rows=['# Energy_eV Im_response Re_response; '+selected_meta['label']+'; '+selected_meta['status']]
                rows += [f'{e:.12g} {z.imag:.15e} {z.real:.15e}' for e,z in zip(energies,values)]
                al.atomic_text(dest/f'{label}_order_{H}_selected.dat','\n'.join(rows)+'\n')
        print('Selected component:',selected_meta)
        print('Normalization:',convention);print('Harmonic units:',units)
        if metrics:print('RT/Floquet comparison:',json.dumps(al.jsafe(metrics),indent=2))
        print('Written:',outfile)

if __name__=='__main__':main()
