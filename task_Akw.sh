#!/bin/sh
#SBATCH --job-name=Akw
#SBATCH --error=Akw.err
#SBATCH --output=Akw.out
##SBATCH --nodelist=pv04
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
##SBATCH --mem-per-cpu=2G

export LC_NUMERIC="en_US.UTF-8"

src="test"
echo $src
for mu in 0.0 1.0 2.0 3.0 3.5 4.0 4.25 4.5 5.0 5.5
do
    dest=$(printf "test.mu%.3f" $mu)
    echo $dest
    pushd $dest
    mpirun -- python ../get_Akw_band_basis.py -wstep=0.02 -nk=128 -skip=[0,1]
    popd
    #src=$dest
done
