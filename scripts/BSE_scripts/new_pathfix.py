"""
This is a second attempt to fix paths after tmp_pathfix
"""
from yambopy import YamboLatticeDB
from yambopy.lattice import red_car,car_red
from yambopy.kpoints import build_ktree,find_kpt
from yambopy.tools.skw import SkwInterpolator
import numpy as np
from copy import deepcopy

def latex(s: str) -> str:
    """
    Print to terminal Greek letters
    written in latex math mode
    """
    import re
    LATEX_GREEK = {
        r'\alpha': 'α', \
        r'\beta': 'β', \
        r'\gamma': 'γ', \
        r'\delta': 'δ',\
        r'\epsilon': 'ε',\
        r'\zeta': 'ζ',\
        r'\eta': 'η',\
        r'\theta': 'θ',\
        r'\iota': 'ι',\
        r'\kappa': 'κ',\
        r'\lambda': 'λ',\
        r'\mu': 'μ',\
        r'\nu': 'ν', r'\xi': 'ξ', r'\omicron': 'ο', r'\pi': 'π',
        r'\rho': 'ρ', r'\sigma': 'σ', r'\tau': 'τ', r'\upsilon': 'υ',
        r'\phi': 'φ', r'\chi': 'χ', r'\psi': 'ψ', r'\omega': 'ω',
        r'\Gamma': 'Γ', r'\Delta': 'Δ', r'\Theta': 'Θ', r'\Lambda': 'Λ',
        r'\Xi': 'Ξ', r'\Pi': 'Π', r'\Sigma': 'Σ', r'\Upsilon': 'Υ',
        r'\Phi': 'Φ', r'\Psi': 'Ψ', r'\Omega': 'Ω',
    }

    # remove dollars
    s = re.sub(r'\$(.*?)\$', r'\1', s)
    # replace latex commands
    for symbol, unicode in LATEX_GREEK.items():
        s = s.replace(symbol, unicode)

    return s

def kfmt(kpt,digits=6):
    """
    Format kpts
    (useful for rlu coords.)
    """
    fmt_kcoords = []
    for i in range(len(kpt)):
        fmt_kcoords.append(round(kpt[i], digits))
    return fmt_kcoords

def SKW_interpolation(qpts_fine_red,exc_disp,LAT,fermie,lpratio):
    """
    This is how the interpolator should be called in yambopy,
    with a single "agnostic" function.

    Working with symmetries in IBZ is *suggested* rather than
    supplying full-BZ values

    Input:
    :: fine q-path for interpolated energies in rlu (attribute of Path object)
    :: noninterpolated energies on original IBZ / BZ grid
    :: YamboLatticeDB object unexpanded (IBZ) / expanded (BZ)
    :: Fermi Energy (can be put as zero for excitons, semiconductors, etc)
    :: interpolation parameter (~number of Fourier components)

    Output:
    :: interpolated energies on fine q-path
    """
    na = np.newaxis

    # Inputs for SKW
    qpts_red      = LAT.red_kpoints
    exc2int       = exc_disp[na,:,:]
    symrel        = [ sym for sym,trev in zip(LAT.sym_rec_red,LAT.time_rev_list) if trev==False ]
    time_rev      = bool(LAT.time_rev)
    cell          = (LAT.lat,LAT.red_atomic_positions,LAT.atomic_numbers)

    # Actual interpolation: find interpolating function
    skw           = SkwInterpolator(lpratio,qpts_red,exc2int,fermie,LAT.nelectrons,cell,symrel,time_rev)
    # Obtain function values on fine grid points
    energies      = skw.interp_kpts(qpts_fine_red).eigens

    return energies[0]

def generate_G_shells(rlat,Nshells=3,unshifted=False):
    """
    Generate G vectors in Cartesian coordinates

    :: rlat -> reciprocal lattice vectors
    :: Nshells -> number of shells generated (in some rare cases, a BZ could border a second-order shell)
    :: unshifted -> if true., add [0.,0.,0.] to the Gvectors
    """
    car_G = []
    for i in range(-Nshells,Nshells+1):
        for j in range(-Nshells,Nshells+1):
            for k in range(-Nshells,Nshells+1):
                if (i, j, k) == (0,0,0): continue
                car_G.append( i*rlat[0]+j*rlat[1]+k*rlat[2] )
    car_G = np.array(car_G)

    if unshifted: car_G = np.insert(car_G,0,np.zeros(3),axis=0)
    return car_G

def point_is_on_border(car_k,rlat,Nshells=3,tol=1e-6):
    """
    Evaluate the condition    k.G = 0.5 |G|^2
    for each kpoint and a certain number of G shells.
    If true, then the kpoint in on the BZ border.

    :: car_k   : kpoints in cartesian coordinates
    :: rlat    : reciprocal lattice vectors
    :: Nshells : shells of G-vectors generated
                 (default 3: should not be increased unless pathological BZ)
    :: tol     : numeric tolerance (in principle depends on kpt density)
    """
    # First build Nshells shells of lattice vectors
    car_G = generate_G_shells(rlat,Nshells=Nshells)
    # Now evaluate Bragg's condition k.G = 0.5*|G|^2/2.
    bragg_l = car_k @ car_G.T  # list of (Nk,NG) scalar products
    bragg_r = 0.5 * np.linalg.norm(car_G,axis=1)**2. # (NG) norms
    bragg_satisfied = np.abs(bragg_l - bragg_r[np.newaxis,:]) < tol # (Nk,NG) bool values: for each k, NG tests
    # Find border points
    border_points_indx = []
    for ik in range(len(car_k)):
        is_border = np.count_nonzero(bragg_satisfied[ik])
        if is_border>0: border_points_indx.append(ik)
    return border_points_indx

class NEW_Path(object):
    """ 
    Class that defines a path in the Brillouin zone
    from user input of special points path.

    It can be used for:

    [1]
    Find a band structure path when our calculations 
    are only in the standard Monkhorst-Pack grid and
    cannot be easily rerun with a `bands` setting.
    This is the case for heavy GW or BSE finite-q runs.

    In this case, we have to find the collinear points
    in the existing grid connecting the supplied special
    points. We have to take care about actually getting 
    the *shortest* paths between special points.

    [2]
    Obtain a band structure path on an unsampled fine grid
    to be the used for (SKW-type) interpolation


    Input format:

        :: Path(klist, intervals, lattice)

        :: klist = [ [[ k1_x, k1_y, k1_z ], 'K1_label'],
                        ...                 ...
                     [[ kN_x, kN_y, kN_z ], 'KN_label'] ]

        :: intervals = [ n_steps_1, n_steps_2, ... , n_steps_N-1 ]

        NB: interval steps are the number of kpoints along
            each high-symmetry direction. If provided 
            consistently with the existing grid, we are
            in case [1]. If any (large) number is provided,
            this is used for interpolation in case [2].

        :: lattice -> YamboLatticeDB object

    Important attributes:

    * klist_opt -> new klist with closest special pts (may be different than user-supplied list)
    * kpoints -> user-supplied special points in rlu
    * klabels -> user-supplied labels
    * kpath -> list of reduced kpt coords. along path
    * kpath_cc -> same as above but Cart. coords.
    * intervals -> as input
    * kdict -> dictionary with information about special pts
    * distances -> cumulative lengths of band segments
    * band_indices -> indices along path in case [1]
    """

    def __init__(self,klist,intervals,lattice,ktree=None,verbose=False,debug=False):
    
        self.intervals = intervals
        self.missing_points = []
        self.warnings(klist,intervals)

        # User-supplied input
        self.kpoints, self.klabels = self.from_klist(klist)

        # Safe copy (deleted below)
        lat = deepcopy(lattice)
    
        # Expand if `lattice` unexpanded
        if lat.ibz_nkpoints == lat.nkpoints:
            lat.expand_kpoints(atol=1e-6,verbose=0)
    
        self.rlat         = lat.rlat
        self.red_kpts     = lat.red_kpoints
        self.car_kpts     = lat.car_kpoints
        self.ibz_red_kpts= lat.get_ibz_kpoints(units='red')
        self.kmap         = lat.BZ_to_IBZ_indexes
        self.inv_kmap     = lat.IBZ_to_BZ_indexes
    
        del lat
       
        if ktree is None: 
            ktree = build_ktree(self.red_kpts)
        self.ktree = ktree

        # [a.] Check special points and grid
        self.kdict = self.check_special_points()
        self.grid, self.dk_step = self.check_grid()

        # [b.] Calculate distances in band segments
        # [BUG] At the moment, if len(self.missing_points)<1 the function get_band_segments doesn't work
        self.distances, self.klist_opt = self.get_band_segments(verbose=debug)
        self.kpoints_opt = self.from_klist(self.klist_opt)[0]

        # [c.] Now generate points along the path
        self.kpath = self.get_kpath()[:,:3]
        self.kpath_cc = red_car(self.kpath,self.rlat)

        # [d.] Get qx axis for plot
        self.qx_axis = self.get_qaxis()
       
        self.sampled_intervals = None
        self.sampled_kpath = None
        self.sampled_kpat_cc = None
        self.sampled_qx_axis = None
        self.band_indices = None
        if len(self.missing_points)<1:
            # [e.] Find collinear points in existing grid
            self.sampled_intervals = self.find_coarse_grid()
            self.sampled_kpath = self.get_kpath(intervals=self.sampled_intervals)[:,:3]
            self.sampled_kpath_cc = red_car(self.sampled_kpath,self.rlat)
            # [f.] Get qx axis for coarse grid
            self.sampled_qx_axis = self.get_qaxis(intervals=self.sampled_intervals,special_pt_count_twice=True)

            # [g.] Get k-indices of calculated energies along the path (i.e., the IBZ equivalent indices)
            self.band_indices = self.find_band_indices(kpath=self.sampled_kpath,special_pt_count_twice=True)

        if verbose: self.print_info()

    def from_klist(self,klist):
        klabels = []
        kpoints = []
        for kpoint, klabel in klist:
            kpoints.append(kpoint)
            klabels.append(klabel)
        return kpoints, klabels

    def check_grid(self):
        """ We find the minimal k-steps and grid size
        """
        def ind_min_pos(a,tol=1e-5):
            return np.where(a>tol,a,np.inf).argmin()

        kx = self.red_kpts[:,0]
        ky = self.red_kpts[:,1]
        kz = self.red_kpts[:,2]
        ksteps = kfmt( [ kx[ind_min_pos(kx)],\
                         ky[ind_min_pos(ky)],\
                         kz[ind_min_pos(kz)] ] )
    
        ksteps = [1. if stp == 0. else stp for stp in ksteps] 
        Ngrid  = [ int(np.round(1./stp)) for stp in ksteps]
        min_dk_rlu = min(ksteps)
        return Ngrid, min_dk_rlu
    
    def find_coarse_grid(self,dk=None,tol=1e-9):
        """ We find the number of collinear sampled points
        """
        if dk is None: dk = self.dk_step

        sampled_intervals = []
        for ipt in range(1,len(self.kpoints_opt)):
            DIST = np.abs(self.kpoints_opt[ipt]-self.kpoints_opt[ipt-1])
            STEP = np.min(DIST[np.abs(DIST) > tol])
            interval = int(np.round(STEP/dk))
            sampled_intervals.append(interval)

        return sampled_intervals

    def get_kpath(self,kpoints=None,intervals=None):
        """
        Output in the format of quantum espresso == [ [kx, ky, kz, 1], ... ]
        """
        if kpoints is None:   
            # Default: get closest special points
            kpoints = self.from_klist(self.klist_opt)[0]
        if intervals is None: intervals = self.intervals

        kout  = np.zeros([sum(intervals)+1,4])
        kout[:,3] = 1
        io = 0
        for ik,interval in enumerate(intervals):
            for ip in range(interval):
                kout[io,:3] = kpoints[ik] + float(ip)/interval*(kpoints[ik+1] - kpoints[ik])
                io = io + 1
        kout[io,:3] = kpoints[ik] + float(ip+1)/interval*(kpoints[ik+1] - kpoints[ik])
        return kout

    def check_special_points(self):
        """
        Performs checks on given list of special points:

        (1) if sampled in expanded BZ grid
        (2) equivalent points
        (3) if it's on the border of BZ

        :: Inputs:

        ylat -> YamboLatticeDB object
        klist -> list of special points in crystal coords. and labels

        NB: points rlu coords. should be |k_i|<1 for all i

        Returns dictionary with special point info
        """

        # Find out which points are on the BZ borders/edges
        border_points = point_is_on_border(self.car_kpts,self.rlat)

        # Now check for the points
        ktree = self.ktree
        
        points_info = []
        for ipt,pt in enumerate(self.kpoints):
            
            # Create dictionary with info
            kdict = {}
            kdict['label'] = self.klabels[ipt]
            # Find point in grid
            try:
                bz_idx  = find_kpt(ktree,pt,tol=1e-5) 
                ibz_idx = self.kmap[bz_idx]
                kdict['bz_index'] = bz_idx
                kdict['bz_coords']= kfmt(self.red_kpts[bz_idx])
                kdict['ibz_index'] = ibz_idx
                kdict['ibz_coords']= self.ibz_red_kpts[ibz_idx]
                kdict['star'] = self.inv_kmap[ibz_idx]
                kdict['is_on_border'] = bz_idx in border_points
            except AssertionError:
                kdict['bz_index'] = None
                kdict['bz_coords']= None
                kdict['ibz_index'] = None
                kdict['ibz_coords']= None
                kdict['star'] = None
                kdict['is_on_border'] = None
                self.missing_points.append(self.klabels[ipt])
                
            points_info.append(kdict)

        return points_info

    def print_info(self):

        kdict = self.kdict
        print("== Analysis of special points supplied ==")
        for i in range(len(kdict)):
            print(f"Point name: {latex(kdict[i]['label'])}")
            print( "Coordinates supplied: ",np.array(self.kpoints[i]))
            print( "Coordinates found in IBZ:   ",kdict[i]['ibz_coords'])
            print(f"On BZ border (face/edge): {kdict[i]['is_on_border']}")
            print(f"Index: {kdict[i]['bz_index']} (BZ), {kdict[i]['ibz_index']} (IBZ)")
            print("===")

        print("== Analysis of collinear sampled points ==")
        print("Grid size: ", self.grid)
        if len(self.missing_points)<1:
            print("Numbers of sampled points collinear with path (intervals): ", self.sampled_intervals)
        else:
            print(f"No sampling possible since point(s) {[latex(pt) for pt in self.missing_points ]} are missing from the kgrid")
        print("==")

    def get_band_segments(self,add_zero_in_front=True,verbose=False):
        """
        Get lengths of band segments to plot. 
        Checks also Gshifted copies of endpoint.

        Returns: 
         - distances -> cumulative segment lengths.
         - new_klist -> klist with path-compliant special points
        """
        nrm = np.linalg.norm
        kdict = self.kdict
        rlat  = self.rlat
        Gshifts = generate_G_shells(self.rlat,unshifted=True,Nshells=1) 
        segments  = [] # length of each "band" segment
        distances = [] # cumulative length
        dist_sum  = 0. # accumulator for distances
        # CHOICE: We start from the IBZ coords of the first point
        pt_start  = red_car([kdict[0]['ibz_coords']],rlat)[0] 
        # klist with closest points to be given to Path obj.
        new_klist = [[car_red([pt_start],rlat)[0],\
                      kdict[0]['label']]]
        for ipt in range(1,len(kdict)):
            # Final point in the IBZ
            if self.kdict[ipt]['label'] in self.missing_points:
                print(f"[WARNING] Point {latex(self.kdict[ipt]['label'])} is not sampled: using supplied coordinates")
                pt_end = red_car([self.kpoints[ipt]],rlat)[0]
            else:   
                pt_end = red_car([kdict[ipt]['ibz_coords']],rlat)[0]
            dist = nrm(pt_end-pt_start)
            PT_end = pt_end
            # FP: need checks in the star of special PT?
            #for istar in kdict[ipt]['star']:
                #tmp_pt_end = car_kpts[istar]
            for G in Gshifts:
                tmp_pt_end = pt_end+G
                dist_G = nrm(tmp_pt_end-pt_start) 
                # If we find a closer special point:
                if dist_G < dist: 
                    #print(dist_G,dist)
                    dist=dist_G
                    PT_end = tmp_pt_end
            red_PT_end = car_red([PT_end],rlat)[0]
            if verbose:
                print('PT start: ',car_red([pt_start],rlat)[0])
                print('PT start (IBZ sampled): ',kdict[ipt-1]['ibz_coords'])
                print("--") 
                print('PT end (closest): ',red_PT_end), 
                print('PT end (IBZ sampled)',car_red([pt_end],rlat)[0])
                print("==")
            # We start again from the new point (equivalent)
            pt_start = PT_end
            dist_sum += dist 
            segments.append(dist)
            # These are the closest distances between special points
            distances.append(dist_sum)

            # Building the new klist to be used in Path.klist
            new_klist.append([red_PT_end,kdict[ipt]['label']])
        
        # Add initial 0 (useful for plotting)
        if add_zero_in_front: distances.insert(0,0.)

        return distances, new_klist

    def get_qaxis(self,intervals=None,special_pt_count_twice=False):
        """ In the case of the coarse grid [1], it can be 
            convenient that each segment restarts counting 
            the initial point (special_pt_count_twice=True).
            For interpolation [2], it shouldn't be used.
        """

        if intervals is None: intervals = self.intervals
        intervals = np.array(intervals,dtype=int)
        
        # Include start point as well
        intervals = intervals + np.ones(len(intervals),dtype=int)
        distances = self.distances
        if distances[0]==0.: distances = distances[1:]

        qx = np.zeros(np.sum(intervals))
        lines_j = 0.
        line_start = 0
        end_pts = []

        for j in range(len(intervals)):
            
            if j==0:
                dstep = distances[j]/(intervals[j]-1.)
            else:
                dstep = (distances[j]-distances[j-1])/(intervals[j]-1.)
            line_end = line_start+intervals[j]
            for i in range(intervals[j]):
                qx[i+line_start] = i*dstep + lines_j
            line_start = line_end
            end_pts.append(line_end)
            lines_j = distances[j]
        end_pts.pop()

        if not special_pt_count_twice: 
            qx = np.delete(qx,end_pts)

        return qx

    def find_band_indices(self,kpath=None,special_pt_count_twice=False):
        """ 
        For each kpoint along the band kpath, which may 
        not be sampled, find the equivalent IBZ index

        The band kpath is the one obtained from the 
        existing Monkhorst-Pack grid

        :: special_pt_count_twice -> True is often easier for plots
        """
        if kpath is None: kpath = self.kpath

        sp_pts_indx = []
        for i in range(len(self.kdict)):
            sp_pts_indx.append(self.kdict[i]['ibz_index'])
        
        indices = []
        for kpt in kpath:
            bz_index = find_kpt(self.ktree,kpt)
            ibz_index = self.kmap[bz_index]
            if special_pt_count_twice:
                if ibz_index in sp_pts_indx:
                    indices.append(ibz_index)
            indices.append(ibz_index)
        if special_pt_count_twice: indices = indices[1:-1]
        return indices

    def warnings(self,klist,intervals):

        if len(intervals)>len(klist)-1:
            print(f"[WARNING] number of intervals is {len(intervals)}, should be {len(klist)-1} (no. of path lines). Taking the first {len(klist)-1} intervals.")
            self.intervals = intervals[:len(klist)-1]
        if len(intervals)<len(klist)-1:
            Nsteps=20
            print(f"[WARNING] number of intervals (path lines) inconsistent with special points specified, using default of {Nsteps} steps per direction")
            self.intervals = [ Nsteps for ik in range(len(klist)-1) ]

if __name__ == "__main__":

    # Tetragonal 3D case
    Z_red = [0.5, 0.5, 0.5]   # Path points in crystal coordinates
    Gamma = [0.0, 0.0, 0.0]   # 0
    X_red = [0.0, 0.5,-0.5]   #
    P_red = [0.25,0.75,-0.25] # 
    N_red = [0.0,1.0,-0.5]    #
    N_red = [0.0,0.0,-0.5]    #
    klist = [[Z_red,'Z'],[Gamma,'$\Gamma$'],[X_red,'X'],[P_red,'P'],[N_red,'N'],[Gamma,'$\Gamma$']]
    intervals =  [20,20,10,10,20]
    # Get kgrid and rlat
    ylat = YamboLatticeDB.from_db_file(filename='./ns.db1')
   
    # Hexagonal 2D case
    #Gamma = [0.0, 0.0, 0.0]   # 0
    #M     = [0.0, 0.5, 0.0]
    #K     = [1./3.,1./3.,0.0]

    #klist = [[Gamma,'$\Gamma$'],[M,'M'],[K,'K'],[Gamma,'$\Gamma$']]
    #intervals = [100,100,100]

    #ylat = YamboLatticeDB.from_db_file(filename='DATABASES_MoS2/SAVE/ns.db1')
    
    path = NEW_Path(klist,intervals,ylat,verbose=True,debug=False)
    print(path.kdict[2]['star'])
    print(ylat.red_kpoints[path.kdict[2]['star']])
    print(ylat.car_kpoints[4095]*ylat.alat[0])
    print(path.distances)
    print(path.sampled_intervals)
    print(path.band_indices)
