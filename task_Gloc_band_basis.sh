#!/bin/sh
#SBATCH --job-name=get_Gloc_band_basis
#SBATCH --error=get_Gloc_band_basis.err
#SBATCH --output=get_Gloc_band_basis.out
#SBATCH --nodes=3
#SBATCH --ntasks=134
#SBATCH --cpus-per-task=1
##SBATCH --mem-per-cpu=2G

mpirun -- python ../get_Gloc_band_basis.py -wstep=0.2
