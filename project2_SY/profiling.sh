#!/bin/bash
#
#
#SBATCH --job-name=CosmoBoost Profiling
#SBATCH --output=profiling.txt
#
#SBATCH --ntasks=5
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=50GB
#SBATCH --constraint=IB
#SBATCH --mail-user=maamari@usc.edu
#SBATCH --mail-type=ALL
#
cd /home/rcf-proj3/ep4/yasini/software/montepython_3/
source /usr/usc/openmpi/default/setup.sh

#source /usr/usc/mpich2/default/setup.sh
#export PYTHONPATH=/home/rcf-proj3/ep4/yasini/miniconda2/lib/python2.7:$PYTHONPATH
#source /home/rcf-proj3/ep4/yasini/software/plc-2.0/bin/clik_profile.sh

stepsize=3.
Nstep=100000
updatestep=50
#echo clik_sourced

#directories
params=input/cmb_fullsky_r.param
output=chains/cmb_fullsky/ell_2-2000_fsky_1/01/r/
#restart=chains/masked_plancknoise_sims/4/rest/f38/ell_2-2500_fsky_65/2018-08-08_100000__1.txt
covmat=covmat/base2015TTTEEE.covmat
#bestfit=bestfit/base2015.bestfit
#export OMP_NUM_THREADS=2




srun --ntasks=5 --mpi=pmi2 python montepython/MontePython.py run \
    --conf default.conf -j fast\
	-p $params -o $output -c $covmat --superupdate 20\
	 -N $Nstep --update $updatestep --silent


source activate cosmoboost
pytho profiling.py
