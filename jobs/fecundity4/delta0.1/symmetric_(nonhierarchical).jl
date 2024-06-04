include("src/extinction_probability_multipatch_functions.jl")
using .ExtinctionMultipatch

fecundity = 4
delta = 0.1
alpha1s = [0.5, 1, 3, 5, 15]
p1 = 1

for alpha1 in alpha1s
    beta1 = alpha1
    
    alpha2 = alpha1
    beta2 = beta1

    save_dir = "hierarchical/fecundity$(fecundity)/delta$(delta)/symmetric_(nonhierarchical)/alpha1_$(alpha1)"

    minimize_extinction_probability(fecundity, delta, alpha1, beta1, alpha2, beta2, p1, save_dir, 100, true, true, 1e-8, 0.5, 75, 100)
end

