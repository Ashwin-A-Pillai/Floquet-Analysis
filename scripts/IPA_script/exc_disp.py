"""
- Interpolate exciton band structure from a regular grid to a path.
- Obtain plot for both original points on the grid and interpolation
- Dump data in auxiliary files to be re-loaded for fast plotting

Inputs:
    - ns.db1 database
    - ndb.BS_diago_Q* databases
    - Path in BZ (crystal coordinates)
"""
import numpy as np
import matplotlib.pyplot as plt
import os.path
from yambopy import YamboLatticeDB
from yambopy.tools.skw import SkwInterpolator
import new_pathfix

def dump_exc_energies(ylat,BSE_path,out_name):
    """
    Dump exc energies to auxiliary file to avoid IO load
    due to opening many ndb.BS databases each time

    :: ylat -> YamboLatticeDB object with `Expand=False`
    """
    energies = []
    print("Reading exciton energies...")
    for iQ in range(ylat.nkpoints):
        iQ_fortran = iQ+1
        print(f'Q: {iQ_fortran}...')
        yexc = YamboExcitonDB.from_db_file(ylat,filename=BSE_path+f'/ndb.BS_diago_Q{iQ_fortran}',Load_WF=False)
        energies.append( yexc.eigenvalues.real )

    energies = np.array(energies)
    print(f'Nq   = {energies.shape[0]}')
    print(f'Nexc = {energies.shape[1]}')
    np.save(out_name,energies)
    return energies

### Input parameters by user ###

# whether to overwrite previously dumped interpolation data
Rerun_Interpolation = False

# Path points in crystal coordinates: use coordinates within (-1,1)
Z_red = [0.5, 0.5, 0.5]   # 
Gamma = [0.0, 0.0, 0.0]   # 
X_red = [0.0, -0.5,-0.5]  #  Any equivalent coordinates are ok, 
P_red = [0.25,-0.25,-0.25] # even if that copy of the special point
N_red = [0.0,0.0,-0.5]    #  is not directly sampled inside the grid
# Points labels and corresponding coordinates
klist = [[Z_red,'Z'],[Gamma,'$\Gamma$'],[X_red,'X'],[P_red,'P'],[N_red,'N'],[Gamma,'$\Gamma$']]
# Interpolation steps per symmetry line (Nlines = Npoints-1)
nsteps = 100
intervals =  [nsteps,nsteps,int(nsteps/2),int(nsteps/2),nsteps]

Nexc     = 30  # Number of exciton branches to interpolate
Nexc2plt = 25  # Number of exciton branches in plot
lpratio = 20   # shell ratio for SKW interpolation (parameter to be adjusted)

# netcdf databases
savepath = './SAVE' # Path to ns.db1 database (I usually download only this on local)
excpath  = './BSE'  # Path to ndb.BS_diago_Q* databases

# auxiliary files
BSE_energies = './exciton_energies.npy' # Path to exciton energies pre-extracted from DBs
out_data_file = f'interpolated_exc_dispersion.dat' # interpolated exciton energies

# output figure files
out_plot_pdf_file = f'ZGXPNG_disp.pdf'
out_plot_png_file = f'ZGXPNG_disp.png'

### End input parameters by user ###

# Load lattice info, IBZ pts & symmetries | do not expand
ylat   = YamboLatticeDB.from_db_file(filename=f'{savepath}/ns.db1',Expand=False)

# Get band path info
# Add `verbose= True` for more information on this step
# Add `debug= True` for even more information
path = new_pathfix.NEW_Path(klist,intervals,ylat,verbose=True)

# Get path of calculated / interpolated points within IBZ (NEW_Path attributes)
# Original points
#car_coarse_qpath     = path.sampled_kpath_cc # Points collinear to the path in original grid (cc)
qpath_coarse_indices = path.band_indices # IBZ indices of the collinear points
qx_coarse            = path.sampled_qx_axis # x-axis of band plot for coarse points
# For interpolation
#car_fine_qpath       = path.kpath_cc # Points in the path for interpolation (cc)
red_fine_qpath       = path.kpath # Points in the path for interpolation (rlu)
qx_fine              = path.qx_axis # x-axis of band plot for interpolated points

# Calculated energies
# If energy auxiliary file present, load it
if os.path.isfile(BSE_energies):
    excdisp = np.load(BSE_energies)[:,:Nexc]
# Otherwise, load data from databases and dump
else:
    # Note: `dump_exc_energies` can be run independently on the cluster in order to pre-produce 
    #       the auxiliary file there, while running this script on the local machine 
    #       (my preferred solution)
    excdisp = dump_exc_energies(ylat,excpath,BSE_energies)[:,:Nexc]
energies_coarse   = excdisp[qpath_coarse_indices] # y-axis of plot for coarse points

# Interpolated energies
# If interpolated data present, load them and go to plot
if os.path.isfile(out_data_file) and not Rerun_Interpolation: 
        energies = np.loadtxt(out_data_file)[:,1:] # y-axis of plot for interpolated points
# Otherwise, interpolation runs and files are dumped
else:
    # SKW interpolation
    energies = new_pathfix.SKW_interpolation(red_fine_qpath,excdisp,ylat,0,lpratio)
    np.savetxt(out_data_file,np.c_[qx_fine,energies],fmt='%.8f',header='# qpt # exciton energies')

#
# Figure
#
# Note: all of the below can be done with an independent script or plotting function
#       It requires as input
#       :: qx_fine, energies (already available in `out_data_file`)
#       :: qx_coarse, energies_coarse
#       :: NEW_Path object (for x-axis layout and labels, not strictly necessary)
#                          (note that an object *can* be dumped in .npy files and read with np.load)
#
# Plot layout options depend on the system and figure
ftsz = 12
ftsz_xticks = ftsz+2
yticks = [3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4]  # Energy range and ticks

fig, ax = plt.subplots(figsize=(8,6)) # Start plot
plt.gca().spines['top'].set_zorder(1) # axis frame in background
plt.gca().spines['right'].set_zorder(1)
plt.gca().spines['bottom'].set_zorder(1)
plt.gca().spines['left'].set_zorder(1)
ax.tick_params(axis='x', which='major', top=True,bottom=True,direction='in') #tick marks inside
ax.tick_params(axis='y', which='major', left=True,right=True,direction='in')

ax.set_ylim(yticks[0],yticks[-1]) # Y axis layout
ax.set_yticks(yticks)
ax.set_yticklabels([str(ytick) for ytick in yticks], fontsize=ftsz)
ax.set_ylabel('Energy (eV)',fontsize=ftsz)

ax.set_xlim(path.distances[0],path.distances[-1]) # X-axis layout using path object
ax.set_xticks([K for K in path.distances])
ax.set_xticklabels(path.klabels,fontsize=ftsz_xticks)

for K in path.distances: ax.axvline(K,color='black') # vertical lines at special points

ax.grid(linestyle='dashed',color='gray') # add grid

# Plot interpolated curves
for i in range(Nexc2plt): ax.plot(qx_fine,energies[:,i],'-',lw=1,c='red',label='BSE interpolation' if i==0 else None)

# Plot calculated points
for i in range(Nexc2plt): ax.plot(qx_coarse,energies_coarse[:,i],'ro',markersize=4.5,zorder=0,label='BSE calculation' if i==0 else None)

plt.legend(fontsize=10) # Final instructions
plt.savefig(out_plot_png_file,dpi=200)
print(f"Saved file at {out_plot_png_file}")
plt.show() # comment this if on remote machine without graphics
