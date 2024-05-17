#!/bin/bash

# Define arrays of parameters
fecundity=4
delta=0.1
pairs=("4,2" "7,3" "20,10" "23,7")
p1s=(0.5 0.6 0.7 0.8 0.9 1)

# Loop over each combination of parameters
for pair in "${pairs[@]}"; do
    for p1 in "${p1s[@]}"; do
        IFS=',' read -r alpha1 beta1 <<< "$pair"

        alpha2=$beta1
        beta2=$alpha1

        save_dir="fecundity4/delta0.1/central_mode__varying_width/alpha1_${alpha1}__beta1_${beta1}__p1_${p1}"
        
        command="julia --project=. src/minimize_extinction_probability_multipatch.jl brute-force-minimize --fecundity $fecundity --delta $delta --alpha1 $alpha1 --beta1 $beta1 --alpha2 $alpha2 --beta2 $beta2 --p1 $p1 --save-dir $save_dir --use-educated-guess"

        echo $command
        
        eval $command
    done
done