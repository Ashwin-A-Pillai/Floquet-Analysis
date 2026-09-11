import numpy as np
import matplotlib.pyplot as plt
from yambopy.exciton_phonon.excph_luminescence import exc_ph_luminescence
from yambopy.exciton_phonon.excph_input_data import exc_ph_get_inputs

# This directory should play the role of "3D_hBN" in the tutorial:
path = '/home/pillai/my_codes/yambo/C-BN/pho/outs/cBN.save'

bsepath    = f'{path}/BSE'       # Lfull BSE (finite-Q)
bseBARpath = f'{path}/BSE-lbar'  # Lbar BSE (Q=0)
elphpath   = path                # where ndb.elph is
dipolespath = bsepath            # where ndb.dipoles is

# Adjust this if SAVE lives somewhere else (see section 1.1)
savepath   = f'{path}/SAVE'      # where ns.db1 is

bands_range   = [1, 7]   # Yambo bands 2-7 (0-based, upper bound excluded)
phonons_range = [0, 6]   # 6 phonon branches for c-BN

nexc_out = 8    # finite-Q excitons
nexc_in  = 16   # Q=0 excitons

T_ph  = 10
T_exc = 10

emin  = 10.3
emax  = 11.0
estep = 0.0002
broad = 0.005

# Load all inputs
input_data = exc_ph_get_inputs(
    savepath, elphpath, bsepath,
    bse_path2=bseBARpath, dipoles_path=dipolespath,
    nexc_in=nexc_in, nexc_out=nexc_out,
    bands_range=bands_range, phonons_range=phonons_range
)

ph_energies, exc_energies, exc_energies_in, G, exc_dipoles = input_data

w, PL = exc_ph_luminescence(
    T_ph, ph_energies, exc_energies, exc_dipoles, G,
    exc_energies_in=exc_energies_in, exc_temp=T_exc,
    nexc_out=nexc_out, nexc_in=nexc_in,
    emin=emin, emax=emax, estep=estep, broad=broad
)

data = np.column_stack((w, PL))
np.savetxt("cBN_luminescence_8x16.dat", data, fmt="%.8f")