#!/bin/sh
##SBATCH --job-name=dmft
#SBATCH --error=do_dmft.err
#SBATCH --output=do_dmft.out
##SBATCH --nodelist=pv04
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
##SBATCH --mem-per-cpu=2G

module purge
module load NRGLjubljana/master-foss-2023b

module list

export PATH=$HOME/NRGRUN/TEST/scripts/tools:$PATH
echo $PATH

export EASYBUILD_MODULES_TOOL=EnvironmentModules
export EASYBUILD_MODULE_SYNTAX=Tcl

export MKL_NUM_THREADS=8 

./DMFT
