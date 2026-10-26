#!/usr/bin/env python3
"""Plot cached Floquet bands, local F-KS/RT projections and full first-BZ weights.

No eigensystem or Fourier fit is recomputed. PNG is standard; --html adds a 3D
interactive BZ. Scalar IBZ weights use YamboPy's saved BZ->IBZ map. Individual
band coefficients need symmetry-covariant labels; complete sums are the default.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import LinearNDInterpolator, RegularGridInterpolator
import FLQ_AL as al


def select_calculation(work, args):
    if args.calc: return work/args.calc
    choices=[]
    for calc in al.calculation_dirs(work,args.calc_glob):
        p=calc/args.out_name/'quasienergies.npz'
        if p.is_file():
            q=al.load_npz(p);e=float(q['omega_au'])*al.HA_TO_EV
            choices.append((abs(e-args.target_eV),e,calc))
    if not choices: raise al.CompareError('No cached extractions; run FLQuasi.py first')
    return min(choices,key=lambda x:x[:2])[2]


def csv_path(path):
    with Path(path).open() as fh:
        reader=csv.DictReader(fh)
        if not {'Label','k1','k2','k3'}.issubset(reader.fieldnames or []):
            raise al.CompareError('Path CSV requires Label,k1,k2,k3 columns')
        rows=list(reader)
    if len(rows)<2: raise al.CompareError('Path CSV needs at least two vertices')
    labels=[r['Label'].strip().replace('GAMMA','Γ') for r in rows]
    labels=['Γ' if x=='G' else x for x in labels]
    red=np.array([[float(r['k1']),float(r['k2']),float(r['k3'])] for r in rows])
    return labels,red


def periodic_interpolate(q,values,points):
    """Linear interpolation on the full periodic mesh; never extrapolate an IBZ."""
    k=q['kpoints_full_red'];shape=np.asarray(q['grid_shape'],int);axes=q['varying_axes'].astype(int)
    ref=k[0];ix=np.rint((k-ref)*shape).astype(int)%shape
    if len(np.unique(ix,axis=0))!=len(k): raise al.CompareError('Full mesh index collision')
    values=np.asarray(values);cube=np.empty(tuple(shape)+values.shape[1:],dtype=values.dtype)
    cube[tuple(ix.T)]=values
    for a in sorted(set(range(3))-set(axes),reverse=True): cube=np.take(cube,0,axis=a)
    pads=[(0,1)]*len(axes)+[(0,0)]*(values.ndim-1)
    periodic=np.pad(cube,pads,mode='wrap')
    query=np.mod(np.asarray(points)-ref,1.0)[:,axes]*shape[axes]
    for a in set(range(3))-set(axes):
        delta=(np.asarray(points)[:,a]-ref[a]+.5)%1-.5
        if np.max(abs(delta))>1e-6: raise al.CompareError('Path leaves the sampled nonperiodic plane')
    interp=RegularGridInterpolator(tuple(np.arange(shape[a]+1) for a in axes),periodic,bounds_error=True)
    return interp(query)


def local_plots(q,f,args,out):
    k=args.kpt-1;active=list(q['active_bands'].astype(int));basis=list(q['basis_bands'].astype(int))
    if not 0<=k<len(q['kpoints_red']): raise al.CompareError('kpt outside the stored IBZ')
    initial=args.initial_band if args.initial_band is not None else active[-1]
    if initial not in active or args.basis_band not in basis:
        raise al.CompareError(f'Active initial bands: {active}; KS basis bands: {basis}')
    mode=int(q['mode_for_band'][k,active.index(initial)]);ib=basis.index(args.basis_band)
    selected=basis if args.projection_bands=='all' else [int(v) for v in args.projection_bands.split(',')]
    if not set(selected).issubset(basis): raise al.CompareError('Projection band outside the actual NLBands window')
    amp=abs(f['Lambda'][k,mode][:,[basis.index(b) for b in selected]].T)
    vmax=max(float(amp.max()),1e-16);floor=max(vmax*1e-12,1e-300)
    fig,ax=plt.subplots(figsize=(8,4))
    im=ax.imshow(np.maximum(amp,floor),origin='lower',aspect='auto',norm=LogNorm(vmin=floor,vmax=max(vmax,floor*10)))
    ax.set_xticks(range(len(f['harmonics'])),f['harmonics']);ax.set_yticks(range(len(selected)),selected)
    ax.set_xlabel('Floquet harmonic η');ax.set_ylabel('KS basis band')
    ax.set_title(f'F-KS projection; k={args.kpt}, mode associated with initial band {initial}')
    fig.colorbar(im,ax=ax,label='|Λ| (amplitude; compact Floquet gauge)')
    fig.tight_layout();fig.savefig(out/'fig-FKS_projection.png',dpi=240);plt.close(fig)
    rawt=f['raw_times_au'];valid=rawt>=float(f['t0_au']);tau=rawt[valid]-float(f['t0_au'])
    rawmode=np.einsum('tbi,ia->tba',f['C_raw'][valid,k],q['Gamma'][k])[:,:,mode]
    e=f['epsilon_compact_au'][k,mode]
    draw=rawmode[:,ib]*np.exp(1j*e*tau)
    T=float(f['period_au']);dense=np.linspace(0,4*T,1201)
    D=np.exp(1j*np.outer(dense,f['harmonics']*float(f['omega_au'])))@f['Lambda'][k,mode,:,ib]
    Ds=np.exp(1j*np.outer(tau,f['harmonics']*float(f['omega_au'])))@f['Lambda'][k,mode,:,ib]
    Dall_sample=np.exp(1j*np.outer(tau,f['harmonics']*float(f['omega_au'])))@f['Lambda'][k,mode]
    Dall_raw=rawmode*np.exp(1j*e*tau)[:,None]
    full_error=al.relnorm(Dall_sample-Dall_raw,Dall_raw)
    print('Periodic-mode reconstruction error (all basis components):',full_error)
    if full_error>1e-6:raise al.CompareError('Cached periodic-mode reconstruction is inconsistent')
    local_scale=float(np.linalg.norm(draw));global_scale=float(np.linalg.norm(Dall_raw))
    print('Selected-component absolute error:',float(np.linalg.norm(Ds-draw)))
    if local_scale>global_scale*1e-12:
        print('Selected-component relative error:',al.relnorm(Ds-draw,draw))
    else: print('Selected-component relative error: not meaningful near numerical zero')
    for part,transform in [('Re',np.real),('Im',np.imag)]:
        fig,ax=plt.subplots(figsize=(8,4))
        ax.plot((dense+float(f['t0_au']))*al.AU_TIME_TO_FS,transform(D),label='Fourier reconstructed mode')
        ax.plot(rawt[valid]*al.AU_TIME_TO_FS,transform(draw),'o',label='Transformed raw checkpoints')
        ax.plot(rawt[valid]*al.AU_TIME_TO_FS,transform(Ds),'x',label='Reconstructed checkpoints')
        for per in (1,4): ax.axvline((float(f['t0_au'])+per*T)*al.AU_TIME_TO_FS,ls='--',lw=.8)
        ax.set(xlabel='Time [fs]',ylabel=f'{part} periodic-mode coefficient',
            title=f'k={args.kpt}, initial-band association={initial}, KS basis={args.basis_band}')
        ax.legend();ax.grid(alpha=.2);fig.tight_layout()
        fig.savefig(out/f'fig-real_time_projection_{part}.png',dpi=240);plt.close(fig)
    print('Mode association weight at selected k:',q['mode_assignment_weight'][k,active.index(initial)])
    # Also show the genuinely propagated initial orbital, not just a demodulated mode.
    j=active.index(initial)
    B=np.exp(1j*np.outer(dense,f['harmonics']*float(f['omega_au'])))
    modes=np.einsum('th,ahb->tba',B,f['Lambda'][k])
    modes*=np.exp(-1j*dense[:,None]*f['epsilon_compact_au'][k][None,:])[:,None,:]
    propagated=np.einsum('tba,aj->tbj',modes,np.linalg.inv(q['Gamma'][k]))[:,ib,j]
    raw_state=f['C_raw'][valid,k,ib,j]
    for part,transform in [('Re',np.real),('Im',np.imag)]:
        fig,ax=plt.subplots(figsize=(8,4))
        ax.plot((dense+float(f['t0_au']))*al.AU_TIME_TO_FS,transform(propagated),label='Floquet-reconstructed propagated orbital')
        ax.plot(rawt[valid]*al.AU_TIME_TO_FS,transform(raw_state),'o',label='Raw propagated orbital')
        ax.set(xlabel='Time [fs]',ylabel=f'{part} C(t)',title=f'Physical state: k={args.kpt}, initial orbital={initial}, KS basis={args.basis_band}')
        ax.legend();ax.grid(alpha=.2);fig.tight_layout()
        fig.savefig(out/f'fig-propagated_state_{part}.png',dpi=240);plt.close(fig)



def band_plot(q,f,labels,vertices,args,out):
    """Compare cached occupied Floquet modes with KS bands on one common k path.

    The physically defined Floquet eigenphases are periodic modulo hbar*omega.
    Interpolate their UNIT-CIRCLE multipliers, not k-dependent compact integer
    Floquet representatives (which create artificial spikes/crossings). For the
    weak-field band comparison, lift those phases to the representative closest
    to each mode's associated KS band. An optional SINGLE GLOBAL energy offset
    then fixes the relative quasienergy/KS zero without modifying splittings.
    """
    B=q['reciprocal_no2pi_bohr']*2*np.pi/al.BOHR_TO_ANG
    pieces=[];dist=[];ticks=[0.0];offset=0.
    for a,b in zip(vertices[:-1],vertices[1:]):
        n=args.path_points;segment=np.linspace(a,b,n,endpoint=False)
        length=float(np.linalg.norm((b-a)@B))
        pieces.append(segment);dist.extend(offset+np.arange(n)*length/n)
        offset+=length;ticks.append(offset)
    points=np.vstack([*pieces,vertices[-1:]]);dist=np.r_[dist,offset]

    # All SAVE KS eigenvalues are available. Only the occupied modes explicitly
    # reconstructed from ndb.RT_V_bands are genuine Floquet states.
    active=np.asarray(q['active_bands'],int)
    mapping=np.asarray(q['mode_for_band'],int)
    nk,nst=mapping.shape
    if len(active)!=nst:
        raise al.CompareError('Floquet mode / propagated occupied-band count mismatch')
    if args.ks_band_min<1 or args.ks_band_max<args.ks_band_min or args.ks_band_max>q['ks_energies_eV'].shape[1]:
        raise al.CompareError('Invalid KS overlay range')
    ix=np.asarray(q['bz_to_ibz'],int)
    omega=float(q['omega_au'])*al.HA_TO_EV
    period=float(q['period_au'])
    ksbands=np.arange(args.ks_band_min-1,args.ks_band_max)
    ks_full=np.asarray(q['ks_energies_eV'])
    ks_plot=periodic_interpolate(q,ks_full[ix][:,ksbands],points)
    ks_assoc=periodic_interpolate(q,ks_full[ix][:,active-1],points)

    raw_principal=np.take_along_axis(np.asarray(q['epsilon_principal_au']).real,mapping,axis=1)
    e_principal_grid=raw_principal*al.HA_TO_EV
    # A global phase-reference difference between Lumen and SAVE may be present.
    # Estimate it once across the BZ; never perform a separate per-k-point offset.
    residual=np.mod(e_principal_grid-ks_full[:,active-1]+omega/2.,omega)-omega/2.
    med_offset=float(np.median(residual))
    spread=float(np.percentile(residual,95)-np.percentile(residual,5))
    reference_offset=med_offset if args.floquet_align=='global' else 0.

    if args.band_gauge=='compact':
        # For reproducibility only; k-dependent integer gauges should NOT be
        # interpolated for a publication-quality KS/Floquet band comparison.
        print('WARNING: compact Floquet gauge varies with k. Linear interpolation '
              'may create artificial sharp spikes. Use --band-gauge ks-connected.')
        compact=np.take_along_axis(f['epsilon_compact_au'].real,mapping,axis=1)
        fl_energy=periodic_interpolate(q,(compact*al.HA_TO_EV)[ix],points)
        min_phase_modulus=float('nan')
    else:
        # Principal Floquet multipliers are invariant under integer gauge
        # shifts epsilon -> epsilon + m omega. Their interpolation across k
        # therefore avoids compact-gauge spikes.
        z_grid=np.exp(-1j*raw_principal*period)
        z_path=periodic_interpolate(q,z_grid[ix],points)
        min_phase_modulus=float(np.min(np.abs(z_path)))
        if min_phase_modulus<0.10:
            raise al.CompareError('Floquet phase interpolation nearly vanishes '
                                  f'(|z|min={min_phase_modulus:.3g}); refine k mesh')
        principal=-np.angle(z_path)/period*al.HA_TO_EV
        if args.band_gauge=='principal':
            fl_energy=principal
        else:  # ks-connected, the default for overlaying KS bands
            nshift=np.rint((ks_assoc+reference_offset-principal)/omega)
            fl_energy=principal+nshift*omega-reference_offset

    if args.energy_zero_eV is None:
        zero=float(q['ks_vbm_eV'])
        zero_label='VBM'
    else:
        zero=float(args.energy_zero_eV)
        zero_label=f'{zero:.6g} eV'
    print(f'Floquet band gauge: {args.band_gauge}; omega={omega:.8f} eV')
    print(f'KS/Floquet phase-reference offset (median, modulo omega): {med_offset:+.8f} eV '
          f'(5th-95th width {spread:.3g} eV)')
    print(f'Applied global reference offset: {reference_offset:+.8f} eV '
          f'(--floquet-align={args.floquet_align}); reference zero: {zero_label}')
    if np.isfinite(min_phase_modulus):
        print(f'Minimum interpolated Floquet multiplier modulus: {min_phase_modulus:.6f}')
    if args.band_gauge=='ks-connected':
        print('Max |Floquet - KS| for represented modes on path:',
              float(np.max(np.abs(fl_energy-ks_assoc))),'eV')

    fig,ax=plt.subplots(figsize=(9.0,6.0))
    if args.ks_overlay:
        for j in range(ks_plot.shape[1]):
            ax.plot(dist,ks_plot[:,j]-zero,color='#d62728',lw=1.3,alpha=.9,
                    label=f'KS bands {args.ks_band_min}–{args.ks_band_max}' if j==0 else None)
    # Only the two propagated occupied modes (bands 3 and 4 for current IPA)
    # are drawn. Replicas change the quasienergy branch by integer hbar omega;
    # they are NOT additional independently computed eigenstates.
    for j in range(nst):
        for replica in range(-args.replicas,args.replicas+1):
            yy=fl_energy[:,j]+replica*omega-zero
            ax.plot(dist,yy,color='#1f77b4',ls='--',lw=1.05,
                    alpha=0.90 if replica==0 else .73,
                    label=(f'Floquet occupied modes + replicas (|m|≤{args.replicas})'
                           if j==0 and replica==-args.replicas else None))
    for tick in ticks:ax.axvline(tick,ls=':',lw=.7,color='0.45')
    ax.set_xticks(ticks,labels);ax.set_xlim(dist[0],dist[-1])
    ax.set_ylabel('Energy relative to VBM [eV]' if args.energy_zero_eV is None
                  else f'Energy relative to {zero_label} [eV]')
    ax.set_title('Floquet and KS bands along '+' → '.join(labels)
                 +f'\n' +rf'$\hbar\omega={omega:.6f}$ eV; represented occupied modes: '
                 +', '.join(map(str,active)))
    ax.legend(loc='lower left',fontsize=8)
    ax.grid(axis='y',alpha=.12)
    fig.tight_layout()
    fig.savefig(out/'fig-Floquet_bands_IBZ.png',dpi=240)
    plt.close(fig)

    # Save the exact values that were plotted, plus unshifted energies.
    al.atomic_npz(out/'bands_path.npz',points_red=points,distance_inv_ang=dist,
        energies_eV=fl_energy,tick_positions=np.array(ticks),
        tick_labels=np.array(labels),gauge=args.band_gauge,
        ks_bands=ksbands+1,ks_energies_eV=ks_plot,
        active_bands=active,ks_associated_energies_eV=ks_assoc,
        energy_reference_eV=zero,global_phase_reference_offset_eV=reference_offset,
        global_phase_reference_residual_median_eV=med_offset,
        global_phase_reference_spread_p05p95_eV=spread,
        minimum_multiplier_modulus=min_phase_modulus,
        floquet_energies_relative_to_reference_eV=fl_energy-zero,
        ks_energies_relative_to_reference_eV=ks_plot-zero)
    print('Band plot (with KS reference):',out/'fig-Floquet_bands_IBZ.png')

def bz_plot(q,f,labels,vertices,args,out):
    harm=f['harmonics'].astype(int);ix=np.flatnonzero(harm==args.eta)
    if len(ix)!=1:raise al.CompareError(f'eta={args.eta} is outside saved support {harm.tolist()}')
    basis=list(q['basis_bands'].astype(int));modes=np.arange(len(q['active_bands']))
    bands=np.arange(len(basis)) if args.bz_bands=='all' else np.array([basis.index(int(v)) for v in args.bz_bands.split(',')])
    if args.bz_modes!='all': modes=np.array([int(v)-1 for v in args.bz_modes.split(',')])
    if modes.min()<0 or modes.max()>=len(q['active_bands']): raise al.CompareError('Invalid 1-based Floquet mode selection')
    per_k_mode=None
    if args.bz_initial_band is not None:
        active=list(q['active_bands'].astype(int))
        if args.bz_initial_band not in active:raise al.CompareError('Requested BZ initial band is not propagated')
        if args.bz_modes!='all':raise al.CompareError('Choose --bz-initial-band or --bz-modes, not both')
        per_k_mode=q['mode_for_band'][:,active.index(args.bz_initial_band)].astype(int)
    resolved=per_k_mode is not None or len(bands)!=len(basis) or len(modes)!=len(q['active_bands'])
    if resolved and not args.assume_band_covariance:
        raise al.CompareError('A band/mode-resolved star map needs sewing information near degeneracies. '
            'Use complete sums (default), or explicitly pass --assume-band-covariance for a diagnostic map.')
    if resolved: print('WARNING: resolved BZ map assumes symmetry-covariant band/mode labels; degeneracies may violate this.')
    coefficients=f['Lambda'][:,:,int(ix[0]),:]
    selected=(coefficients[np.arange(len(coefficients)),per_k_mode][:,None,:]
              if per_k_mode is not None else coefficients[:,modes,:])
    values=np.sum(abs(np.take(selected,bands,axis=2))**2,axis=(1,2))
    if args.bz_quantity=='amplitude':values=np.sqrt(values)
    expanded=values[q['bz_to_ibz'].astype(int)]
    B=q['reciprocal_no2pi_bohr'];axes=q['varying_axes'].astype(int);dim=len(axes)
    bz,hull,frame=al.first_bz(B,axes)
    kcart=q['kpoints_full_red']@B*2*np.pi/al.BOHR_TO_ANG
    imagek,imagev=al.periodic_images(q['kpoints_full_red'],expanded,B,axes,span=2)
    pts=imagek@frame;tol=1e-8
    inside=np.all(pts@hull.equations[:,:dim].T+hull.equations[:,dim]<=tol,axis=1)
    kplot=pts[inside];vplot=imagev[inside]
    rounded=np.round(kplot,10);_,unique=np.unique(rounded,axis=0,return_index=True)
    kplot=kplot[unique];vplot=vplot[unique]
    label=('Σ|Λ|²' if args.bz_quantity=='weight' else '√Σ|Λ|²')+f', η={args.eta}, compact gauge'
    if dim==2:
        xmin,ymin=bz.min(axis=0);xmax,ymax=bz.max(axis=0)
        x=np.linspace(xmin,xmax,args.bz_grid);y=np.linspace(ymin,ymax,args.bz_grid)
        X,Y=np.meshgrid(x,y);target=np.c_[X.ravel(),Y.ravel()]
        mask=np.all(target@hull.equations[:,:2].T+hull.equations[:,2]<=tol,axis=1)
        interp=LinearNDInterpolator(pts,imagev)
        vals=np.full(len(target),np.nan);vals[mask]=interp(target[mask])
        if not np.isfinite(vals[mask]).all():raise al.CompareError('Periodic interpolation failed inside first BZ; no nearest-value extrapolation used')
        fig,ax=plt.subplots(figsize=(7,6));im=ax.pcolormesh(X,Y,vals.reshape(X.shape),shading='auto')
        closed=np.vstack([bz,bz[0]]);ax.plot(closed[:,0],closed[:,1],lw=1.)
        if args.show_kpoints:ax.scatter(kplot[:,0],kplot[:,1],s=3,alpha=.4)
        path=(vertices@B*2*np.pi/al.BOHR_TO_ANG)@frame
        ax.plot(path[:,0],path[:,1],lw=1);ax.scatter(path[:,0],path[:,1],s=16)
        for lab,xy in zip(labels,path):ax.annotate(lab,xy,xytext=(4,4),textcoords='offset points')
        ax.set_aspect('equal');ax.set(xlabel='k₁ in plotting plane [Å⁻¹]',ylabel='k₂ in plotting plane [Å⁻¹]',title='Full first Brillouin zone')
        fig.colorbar(im,ax=ax,label=label);fig.tight_layout();fig.savefig(out/'fig-FKS-proj_BZ.png',dpi=240);plt.close(fig)
    else:
        fig=plt.figure(figsize=(8,6));ax=fig.add_subplot(111,projection='3d')
        im=ax.scatter(kplot[:,0],kplot[:,1],kplot[:,2],c=vplot,s=9,alpha=.8)
        face_by_edge={}
        for face,tri in enumerate(hull.simplices):
            for i,j in ((0,1),(1,2),(2,0)):
                edge=tuple(sorted((int(tri[i]),int(tri[j]))));face_by_edge.setdefault(edge,[]).append(face)
        edges=[]
        for (i,j),faces in face_by_edge.items():
            if len(faces)==2 and np.allclose(hull.equations[faces[0],:3],hull.equations[faces[1],:3],atol=1e-7):continue
            line=bz[[i,j]];ax.plot(*line.T,lw=.7);edges.append(line)
        ax.set(xlabel='kx [Å⁻¹]',ylabel='ky [Å⁻¹]',zlabel='kz [Å⁻¹]',title='Full 3D first Brillouin zone')
        ax.set_box_aspect(np.ptp(bz,axis=0));fig.colorbar(im,ax=ax,shrink=.7,label=label)
        fig.tight_layout();fig.savefig(out/'fig-FKS-proj_BZ.png',dpi=240);plt.close(fig)
        if args.html:
            import plotly.graph_objects as go
            fig=go.Figure(go.Scatter3d(x=kplot[:,0],y=kplot[:,1],z=kplot[:,2],mode='markers',
                marker={'size':3,'color':vplot,'showscale':True,'colorbar':{'title':label}},name='F-KS weight'))
            for line in edges:fig.add_trace(go.Scatter3d(x=line[:,0],y=line[:,1],z=line[:,2],mode='lines',showlegend=False))
            fig.update_layout(title=label,scene={'aspectmode':'data','xaxis_title':'kx [Å⁻¹]','yaxis_title':'ky [Å⁻¹]','zaxis_title':'kz [Å⁻¹]'})
            fig.write_html(out/'fig-FKS-proj_BZ.html',include_plotlyjs=True)
    al.atomic_npz(out/'FKS_BZ.npz',kpoints_full_cart_inv_ang=kcart,values_full=expanded,
        bz_to_ibz=q['bz_to_ibz'],bz_vertices=bz,plot_frame=frame,quantity=label,
        selected_modes_1based=(per_k_mode[:,None]+1 if per_k_mode is not None else modes+1),
        selected_bands=q['basis_bands'][bands],symmetry_caveat='resolved-band covariance explicitly assumed' if resolved else 'complete saved band/mode sums; fixed compact Floquet gauge')
    print(f'BZ expansion: {len(values)} stored points -> {len(expanded)} full mesh points; first-BZ clipping applied.')
    print('This is a Floquet coefficient weight, not an exciton probability or an excitonic wavefunction.')


def spectrum_plots(path,component,out):
    """Plot a field contraction, or a tensor element isolated by symmetry."""
    z=al.load_npz(path);H=int(z['harmonic_order']);units=str(z['units'])
    selection=al.parse_tensor_component(component,H)
    catalog=json.loads(str(z['tensor_catalog_json'])) if 'tensor_catalog_json' in z else None
    if selection['kind']=='tensor' and catalog is None:
        raise al.CompareError('This harmonic NPZ predates tensor diagnostics; '
                              'rerun FLQ_HG.py with the updated four-script suite')
    label='chi_'+selection['label']+'_eff' if selection['kind']=='effective' else selection['label']
    slug=''.join(x for x in label if x.isalnum() or x=='_')

    def tex_symbol():
        if selection['kind']=='effective':
            comp='xyz'[selection['output_axis']]
            return rf'\chi^{{({H})}}_{{\mathrm{{eff}},{comp}}}'
        comp=''.join('xyz'[i] for i in selection['index'])
        return rf'\chi^{{({H})}}_{{{comp}}}'

    def ylabel(kind):
        chi=tex_symbol()
        unit=f' ({units})'
        if kind=='magnitude':
            return rf'$|{chi}|${unit}'
        if kind=='real':
            return rf'$\mathrm{{Re}}({chi})${unit}'
        if kind=='imaginary':
            return rf'$\mathrm{{Im}}({chi})${unit}'
        raise ValueError(kind)

    def extract(array):
        if catalog is None:
            return np.asarray(array)[:,selection['output_axis']]
        values,info=al.tensor_response(catalog,component,array)
        return values

    hasrt='rt_energy_eV' in z;hasfl='fl_energy_eV' in z
    out.mkdir(parents=True,exist_ok=True)
    if catalog is not None:
        print('Tensor symmetry source:',catalog['symmetry_source'])
        print('Applied field:',catalog['field_unit_xyz'])
        if selection['kind']=='tensor':
            print('Requested tensor component:',catalog['components'][selection['key']])

    fig,ax=plt.subplots(figsize=(7.5,4.8))
    for name,flag in [('rt',hasrt),('fl',hasfl)]:
        if flag:
            ax.plot(z[name+'_energy_eV'],np.abs(extract(z[name+'_chi'])),
                ls='-' if name=='rt' else '--',marker=None if name=='rt' else '.',
                label='RT' if name=='rt' else 'Floquet')
    ax.set(xlabel='Driving photon energy [eV]', ylabel=ylabel('magnitude'))
    ax.legend();ax.grid(alpha=.2);fig.tight_layout()
    fig.savefig(out/f'HG_order_{H}_{slug}_magnitude.png',dpi=240);plt.close(fig)

    fig,(axr,axi)=plt.subplots(2,1,figsize=(7.5,7.2),sharex=True)
    for name,flag in [('rt',hasrt),('fl',hasfl)]:
        if flag:
            energy=z[name+'_energy_eV']
            vals=extract(z[name+'_chi'])
            style=dict(ls='-' if name=='rt' else '--',marker=None if name=='rt' else '.',
                       label='RT' if name=='rt' else 'Floquet')
            axr.plot(energy,np.real(vals),**style)
            axi.plot(energy,np.imag(vals),**style)
    axr.set_ylabel(ylabel('real'))
    axi.set_ylabel(ylabel('imaginary'))
    axi.set_xlabel('Driving photon energy [eV]')
    axr.legend();axr.grid(alpha=.2);axi.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(out/f'HG_order_{H}_{slug}_complex.png',dpi=240);plt.close(fig)

    if hasrt and hasfl:
        method=str(z.get('rt_reference_method','exact'))
        if 'rt_matched_chi' in z:
            a_match=z['rt_matched_chi']
        else:
            a_match,_,_=al.reference_on_floquet_grid(
                z['rt_energy_eV'],z['rt_chi'],z['fl_energy_eV'],method=method)
        a=extract(a_match);b=extract(z['fl_chi'])
        metrics=al.comparison_metrics(z['fl_energy_eV'],a[:,None],z['fl_energy_eV'],b[:,None],0)
        metrics['component']=label
        metrics['reference_matching']='exact_RT_subset' if method=='exact' else 'PCHIP_interpolated_RT'
        peak=max(float(abs(extract(z['rt_chi'])).max()),1e-300)
        fig,ax=plt.subplots(figsize=(7.5,4.0))
        ax.plot(z['fl_energy_eV'],100*(abs(b)-abs(a))/peak,marker='o')
        ax.set(xlabel='Driving photon energy [eV]',ylabel='Magnitude difference / RT peak [%]')
        ax.grid(alpha=.2);fig.tight_layout()
        fig.savefig(out/f'HG_order_{H}_{slug}_difference.png',dpi=240);plt.close(fig)
        al.atomic_text(out/f'HG_order_{H}_{slug}_metrics.json',
                       json.dumps(al.jsafe(metrics),indent=2)+'\n')
        print(json.dumps(al.jsafe(metrics),indent=2))

def frequency_plan_plot(path,out):
    """Plot the adaptive RT-driven energy mesh without reopening Lumen databases."""
    j=json.loads(Path(path).read_text())
    E=np.asarray(j['source_energy_eV'],float)
    amp=np.asarray(j['source_magnitude'],float)
    F=np.asarray(j['energies_eV'],float)
    residual=np.asarray(j['pointwise_relative_interpolation_error'],float)
    fig,axes=plt.subplots(2,1,figsize=(9,6.2),sharex=True,
                          gridspec_kw={'height_ratios':[2.6,1]})
    ax,bx=axes
    ax.plot(E,amp,lw=1.7,label=f'RT complex-response magnitude ({j["component"]})')
    ax.scatter(F,np.interp(F,E,amp),s=22,zorder=3,label=f'{len(F)} selected FL frequencies')
    ax.set_ylabel('RT |χ| (response units)');ax.grid(alpha=.22);ax.legend()
    bx.plot(E,residual,lw=1.4,label='Linear interpolation mismatch / local scale')
    for energy in F:bx.axvline(float(energy),color='gray',alpha=.11,lw=.65)
    bx.set_xlabel('Driving photon energy [eV]');bx.set_ylabel('Estimated error')
    bx.grid(alpha=.22)
    if 'original_RT_energy_eV' in j and j.get('method')=='PCHIP_surrogate_off_grid':
        ax.scatter(j['original_RT_energy_eV'],j['original_RT_magnitude'],
                   s=8,marker='x',color='black',alpha=.5,label='Original RT nodes')
        ax.legend()
    fig.suptitle('Adaptive Floquet frequencies from RT response'
                 +(' (PCHIP proposal)' if j.get('method')=='PCHIP_surrogate_off_grid' else ''))
    fig.tight_layout()
    out.mkdir(parents=True,exist_ok=True)
    png=out/'FLQ_adaptive_frequency_grid.png'
    fig.savefig(png,dpi=220);plt.close(fig)
    print('Wrote',png)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workdir',default='.');p.add_argument('--out-name',default='PyFLQ')
    p.add_argument('--calc');p.add_argument('--calc-glob',default='FLQ_[0-9][0-9][0-9]')
    p.add_argument('--target-eV',type=float,default=4.1);p.add_argument('--ibz')
    p.add_argument('--kpt',type=int,default=6);p.add_argument('--initial-band',type=int)
    p.add_argument('--basis-band',type=int,default=5);p.add_argument('--projection-bands',default='all')
    p.add_argument('--eta',type=int,default=2);p.add_argument('--bz-bands',default='all');p.add_argument('--bz-modes',default='all')
    p.add_argument('--bz-initial-band',type=int,help='Per-k mode association for a chosen initial occupied band')
    p.add_argument('--bz-quantity',choices=['weight','amplitude'],default='weight')
    p.add_argument('--assume-band-covariance',action='store_true');p.add_argument('--show-kpoints',action='store_true')
    p.add_argument('--bz-grid',type=int,default=250);p.add_argument('--html',action='store_true')
    p.add_argument('--band-gauge',choices=['principal','compact','ks-connected'],default='ks-connected')
    p.add_argument('--replicas',type=int,default=2);p.add_argument('--path-points',type=int,default=30)
    p.add_argument('--ks-overlay',action='store_true',default=True)
    p.add_argument('--no-ks-overlay',action='store_false',dest='ks_overlay')
    p.add_argument('--floquet-align',choices=['global','none'],default='global',
                   help='single global Floquet/KS energy-reference alignment for comparison')
    p.add_argument('--ks-band-min',type=int,default=1)
    p.add_argument('--ks-band-max',type=int,default=8)
    p.add_argument('--energy-zero-eV',type=float,default=None,
                   help='override default KS VBM zero (e.g. 0 for absolute SAVE energies)')
    p.add_argument('--spectra',help='HG_order_N.npz; plot harmonic spectra instead of state plots')
    p.add_argument('--frequency-plan',help='FLQ_frequencies.json; preview adaptive frequency selection')
    p.add_argument('--component',default='y',help='x/y/z field-contracted response, or tensor indices like [1,1,2] for --spectra')
    args=p.parse_args(argv);work=Path(args.workdir).expanduser().resolve()
    if args.frequency_plan:
        plan=Path(args.frequency_plan);plan=plan if plan.is_absolute() else work/plan
        frequency_plan_plot(plan,work/args.out_name/'plots');return
    if args.spectra:
        path=Path(args.spectra);path=path if path.is_absolute() else work/path
        with al.run_log(work/args.out_name/'FLQPlot_spectra.log'):
            spectrum_plots(path,args.component,work/args.out_name/'plots')
        return
    if not args.ibz:p.error('--ibz is required for state/BZ/path plots')
    calc=select_calculation(work,args);out=calc/args.out_name/'plots';out.mkdir(parents=True,exist_ok=True)
    with al.run_log(calc/args.out_name/'FLQPlot.log'):
        q,f,meta=al.load_pair(calc,args.out_name)
        labels,vertices=csv_path(work/args.ibz)
        print(f'Calculation: {calc.name}; photon energy: {float(q["omega_au"])*al.HA_TO_EV:.8f} eV')
        local_plots(q,f,args,out);band_plot(q,f,labels,vertices,args,out);bz_plot(q,f,labels,vertices,args,out)
        al.atomic_text(out/'plot_settings.json',json.dumps(al.jsafe(vars(args)),indent=2)+'\n')
        print('Plots:',out)

if __name__=='__main__':main()
