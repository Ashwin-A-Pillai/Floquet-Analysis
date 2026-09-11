#!/usr/bin/env python3

import numpy as np
from yambopy import *
from yambopy.plot  import *

NLDB=YamboNLDB()

pol =NLDB.Polarization[0]
time=NLDB.IO_TIME_points
t_initial=NLDB.Efield[0]["initial_time"]
pol_damped=np.empty_like(pol)
for i_d in range(3):
   pol_damped[i_d,:]=damp_it(pol[i_d,:],time,t_initial,damp_type='LORENTZIAN',damp_factor=0.1/ha2ev)
Plot_Pol_or_Curr(time=time, pol=pol_damped, xlim=[0,55], save_file='polarization.pdf')   
Linear_Response(time=time,pol=pol_damped,efield=NLDB.Efield[0],plot=False,plot_file='eps.pdf')