#!/bin/sh
#SBATCH --job-name=dmft
#SBATCH --error=dmft.err
#SBATCH --output=dmft.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --mem-per-cpu=2G

export LC_NUMERIC="en_US.UTF-8"

declare -a Ts
Ts=(0.1 0.07 0.05 0.03 0.02 0.01 0.007 0.005)
echo these are Ts "${Ts[@]}"

declare -a mus

for model in emery hubbard
do    
    echo this is the model $model
    if [ model=="emery" ]; then
        mus=(4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5)
    fi
    if [ model=="hubbard" ]; then
        mus=(5.7 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8 6.9 7.0)
    fi   
    echo this is mus "${mus[@]}"
    src="test.$model.seed"
    echo src=$src
    for mu in "${mus[@]}"
    do
        for T in "${Ts[@]}" 
        do    
            dest=$(printf "comp11a.$model.mu%.3f.T%.3f" $mu $T)
            echo dest=$dest
            python prepare_run.py -source=$src -dest=$dest -keys=["'mu'","'T'"] -vals=[$mu,$T] -clean_up -from_scratch
            pushd $dest
            sbatch -J $dest ../do_dmft.sh 
            popd            
        done
    done
done
