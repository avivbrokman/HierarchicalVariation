include("src/extinction_probability_multipatch_functions.jl")
using .ExtinctionMultipatch

fecundity = 4
delta = 0.1
alpha1s = [1.2, 1.8, 2, 5]
p1s = [0.5, 0.6, 0.7, 0.8, 0.9, 1]

for alpha1 in alpha1s
    for p1 in p1s
        beta1 = 1
        
        alpha2 = beta1
        beta2 = alpha1

        save_dir = "hierarchical/fecundity$(fecundity)/delta$(delta)/extreme_mode__varying_slope/alpha1_$(alpha1)__p1_$(p1)"

        minimize_extinction_probability(fecundity, delta, alpha1, beta1, alpha2, beta2, p1, save_dir, 100, true, true, 1e-8, 0.5, 75, 100)
    end
end

