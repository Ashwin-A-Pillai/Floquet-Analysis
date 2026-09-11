module purge
module load oneapi_2021.2.0
module load yambo_5.3.0

#ypp -F trial.in -J GW0
#mv "o-GW0.bands_interpolated" "bands[KS].out"
#python3 band.py
#python3 band_plt.py "bands[KS].out" "bands-KS.png"

yambo -o b -k sex -y h -V qp -F BSE.in -J GW0
yambo -F BSE.in -J GW0
mv "o-GW0.eps_q1_haydock_bse" "BSE[abs].out"
python3 abs.py