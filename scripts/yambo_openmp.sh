#/bin/csh -f

# Variables d'environnement SGE
#$ -S /bin/csh
#$ -cwd
#$ -j y
#$ -pe qlogic 96
#$ -q ib8.q

setenv OMP_NUM_THREADS 1
setenv NUM_MPI `expr ${NSLOTS} / ${OMP_NUM_THREADS}`
setenv UCX_TLS ud,sm,self 

echo ""
echo "Total number of cores $NSLOTS"
echo "Run with $OMP_NUM_THREADS threads/process"
echo "Run with $NUM_MPI mpi processes"
echo ""

module purge
module load oneapi_2021.2.0
module load yambo_5.3.0


awk '{print $1}' $PE_HOSTFILE > filehosts_${JOB_ID}

echo "Hosts: "
cat filehosts_${JOB_ID}

# newx to run 
mpirun -genv OMP_NUM_THREADS ${OMP_NUM_THREADS} -f filehosts_${JOB_ID} -perhost 1 -genv I_MPI_PIN_DOMAIN omp:scatter -np $NUM_MPI yambo -F GW.in  -J GW0 > output
ypp -s b -V qp -F interplt.in -J GW0
ypp -F interplt.in -J GW0
mv "o-GW0.bands_interpolated" "bands[KS].out"
python3 band_plt.py "bands[KS].out" "bands-KS.png"
ypp -F interplt.in -J GW0
mv "o-GW0.bands_interpolated" "bands[GW].out"
python3 band_plt.py "bands[GW].out" "bands-GW.png"
python3 compr_bnd.py "bands[GW].out" "bands[KS].out" "bandstructure.png"

yambo -o b -k sex -y h -V qp -F BSE.in -J GW0
yambo -F BSE.in -J GW0
mv "o-GW0.eps_q1_haydock_bse" "BSE[abs].out"
python3 abs.py