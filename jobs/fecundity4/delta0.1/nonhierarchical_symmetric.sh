#!/bin/bash

# Define arrays of parameters
fecundity=4
delta=0.1
alpha1s=(0.5 1 3 5 15)
p1=1

# Loop over each combination of parameters
for alpha1 in "${alpha1s[@]}"; do

    beta1=$alpha1

    alpha2=$beta1
    beta2=$alpha1

    save_dir="fecundity4/delta0.1/nonhierarchical_symmetric/alpha1_${alpha1}"
    
    command="julia --project=. src/minimize_extinction_probability_multipatch.jl brute-force-minimize --fecundity $fecundity --delta $delta --alpha1 $alpha1 --beta1 $beta1 --alpha2 $alpha2 --beta2 $beta2 --p1 $p1 --save-dir $save_dir --use-educated-guess"

    echo $command
    
    eval $command
done

