#!/bin/sh
#SBATCH --job-name=dmft
#SBATCH --error=dmft.err
#SBATCH --output=dmft.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
##SBATCH --mem-per-cpu=2G

export LC_NUMERIC="en_US.UTF-8"

#src="test.emery.seed"
src="test.hubbard.seed"
echo $src
#for mu in 0 1 2 3 3.5 4.0 4.25 4.5 5.0 5.5
for mu in 3.25 3.5 3.75 4.0 4.25 4.5 
do
    dest=$(printf "test.hubbard.mu%.3f" $mu)
    echo $dest
    python prepare_run.py -source=$src -dest=$dest -keys=["'mu'"] -vals=[$mu] -clean_up -from_scratch
    pushd $dest
    sbatch ../do_dmft.sh -J $dest
    popd
    #src=$dest
done
