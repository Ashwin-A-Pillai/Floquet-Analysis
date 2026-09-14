from yambopy import *
NLDB=YamboNLDB()
pol =NLDB.Polarization[0]
time=NLDB.IO_TIME_points
SIN = Xn_from_sine(NLDB,X_order=5)
print(SIN)
OUT = SIN.perform_analysis()
SIN.output_analysis(OUT)