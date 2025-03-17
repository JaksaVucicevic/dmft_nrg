#!/bin/sh
#SBATCH --job-name=Gloc
#SBATCH --error=Gloc.err
#SBATCH --output=Gloc.out
##SBATCH --nodelist=pv04
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
##SBATCH --mem-per-cpu=2G

export LC_NUMERIC="en_US.UTF-8"

src="test"
echo $src
for mu in 4.3 4.4 4.6 4.7 4.8 4.9 5.1 5.2 5.3 0 1 2 3 3.5 4.0 4.25 4.5 5.0 5.5 #0.0 1.0 2.0 3.0 3.5 4.0 4.25 4.5 5.0 5.5
do
    dest=$(printf "test.emery.mu%.3f" $mu)
    echo $dest
    pushd $dest
    mpirun -- python ../get_Gloc.py -wstep=0.02
    #mpirun -- python ../get_Gloc_band_basis.py -wstep=0.02 -skip=[1]
    popd
    #src=$dest
done
