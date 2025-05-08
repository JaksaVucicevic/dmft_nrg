#!/bin/sh
#SBATCH --job-name=Gloc
#SBATCH --error=Gloc.err
#SBATCH --output=Gloc.out
##SBATCH --nodelist=pv04
#SBATCH --nodes=2
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
##SBATCH --mem-per-cpu=2G

export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

export LC_NUMERIC="en_US.UTF-8"

comp=1

declare -a Ts
Ts=(0.1 0.07 0.05 0.03 0.02 0.01 0.007 0.005 0.003)
echo these are Ts "${Ts[@]}"

declare -a mus
#mus=(3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0) # 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8)
#mus=(3.3 3.4 3.5 3.6 3.7 3.8)
mus=(2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2)

declare -a Us
Us=(7.0 6.0)

echo this is mus "${mus[@]}"

for U in "${Us[@]}"
do
    for mu in "${mus[@]}"
    do
        for T in "${Ts[@]}" 
        do    
            #dest=$(printf "comp$comp.emery.mu%.3f.T%.3f" $mu $T)
            dest=$(printf "comp$comp.emery.mu%.3f.T%.3f.U%.3f" $mu $T $U)
            echo this is dest $dest
            pushd $dest
            pwd
            mpirun -- python ../get_Gloc.py -wstep=0.002
            #mpirun -- python ../get_Gloc_band_basis.py -wstep=0.02 -skip=[1]
            popd
        done    
    done
done
