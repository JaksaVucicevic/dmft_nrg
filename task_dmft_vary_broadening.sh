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
Ts=(0.1 0.07 0.05 0.03 0.02 0.01 0.007 0.005 0.003)
echo these are Ts "${Ts[@]}"

declare -a mus

comp=1
emery=emery
hubbard=hubbard

for model in $emery $hubbard
do    
    echo this is the model $model
    if [ "$model" = "$emery" ]; then
        echo indeed $emery
        #mus=(3.3 3.4 3.5 3.6 3.7 3.8)
        #mus=(3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8)
        mus=(3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0)
    fi
    if [ "$model" = "$hubbard" ]; then
        echo indeed $hubbard
        #mus=(4.8 4.9 5.0 5.1 5.2 5.3)
        #mus=(5.4 5.5 5.6 5.7 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8 6.9 7.0 7.1 7.2 7.3)
        mus=(4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5)
    fi   
    echo this is mus "${mus[@]}"

    for ba in 0.01 0.05 0.1
    do
        for mu in "${mus[@]}"
        do
            for T in "${Ts[@]}" 
            do    
                src=$(printf "comp$comp.$model.mu%.3f.T%.3f" $mu $T)
                echo src=$src
                dest=$(printf "comp$comp.$model.mu%.3f.T%.3f.ba%.3f" $mu $T $ba)
                echo dest=$dest
                python prepare_run.py -source=$src -dest=$dest -keys=["'broaden_alpha'"] -vals=[$ba] -clean_up
                pushd $dest
                sbatch -J $dest ../do_dmft.sh 
                popd            
            done
        done
    done
done
