#!/usr/bin/env python3

from yambopy import YamboVbandsDB, VbPP
import numpy as np

calc = "fl_sample"

mydb = YamboVbandsDB(calc=calc)
print(mydb)

if mydb.fl_order <= 0:
    raise RuntimeError(
        f"Invalid Floquet order in database: {mydb.fl_order}. "
        "The RT calculation did not generate a valid Floquet sampling database."
    )

if mydb.n_timesteps < 2 * mydb.fl_order + 2:
    raise RuntimeError(
        f"Insufficient Floquet samples: database contains "
        f"{mydb.n_timesteps} timesteps for FLOrder={mydb.fl_order}."
    )

kpt, band = 6, 4

vbs = VbPP(mydb, kpt, band)
print(vbs)

optimized_evecs = vbs.find_qe()

if (
    optimized_evecs.nr_it > vbs.max_iter
    or optimized_evecs.err > vbs.err_thrs
    or optimized_evecs.nr_acc > vbs.qe_thrs
    or not np.isfinite(optimized_evecs.FL_qe)
):
    raise RuntimeError(
        "Floquet quasienergy optimization failed:\n"
        f"  QE                = {optimized_evecs.FL_qe} eV\n"
        f"  iterations        = {optimized_evecs.nr_it}\n"
        f"  QE accuracy       = {optimized_evecs.nr_acc} eV\n"
        f"  periodicity error = {optimized_evecs.err}\n"
        f"  target QE tol     = {vbs.qe_thrs}\n"
        f"  target error tol  = {vbs.err_thrs}"
    )

print(
    f"QE= {optimized_evecs.FL_qe} eV "
    f"with accuracy {optimized_evecs.nr_acc} eV "
    f"after {optimized_evecs.nr_it} iterations "
    f"- error in periodicity = {optimized_evecs.err}"
)

# Projection onto physical band 5
plot_band = 5
plot_index = plot_band - mydb.basis_index[0]

optimized_evecs.plot_realtime(
    band_to_plot=plot_index,
    t_step=0.004,
)

optimized_evecs.plot_floquet()
optimized_evecs.output()