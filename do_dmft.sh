#!/bin/sh
##SBATCH --job-name=dmft
#SBATCH --error=do_dmft.err
#SBATCH --output=do_dmft.out
#SBATCH --nodelist=pv[15-20]
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --mem-per-cpu=2G

module purge
module load NRGLjubljana/master-foss-2023b

module list

export PATH=$HOME/NRGRUN/TEST/scripts/tools:$PATH
echo $PATH

export EASYBUILD_MODULES_TOOL=EnvironmentModules
export EASYBUILD_MODULE_SYNTAX=Tcl

export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

./DMFT
#python ../get_Gloc.py -wstep=0.02
#python ../get_rhodc.py #figure out T from param.loop and model based on whether Phi.dat exists
