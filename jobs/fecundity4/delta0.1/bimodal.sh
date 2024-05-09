#!/bin/bash

# Define arrays of parameters
fecundity=4
delta=0.1
alpha1s=(0.6 0.9)
p1s=(0.5 0.6 0.7 0.8 0.9 1)

# Loop over each combination of parameters
for alpha1 in "${alpha1s[@]}"; do
    for p1 in "${p1s[@]}"; do

        beta2=$alpha1
        beta1=$(echo "1 - $alpha1" | bc -l)
        alpha2=$beta1

        save_dir="bimodal/alpha1_${alpha1}__p1_${p1}"
        
        command="julia --project=. src/minimize_extinction_probability_multipatch.jl --fecundity $fecundity --delta $delta --alpha1 $alpha1 --beta1 $beta1 --alpha2 $alpha2 --beta2 $beta2 --p1 $p1 --save_dir $save_dir"

        echo $command
        
        eval $command
    done
done
