include("src/extinction_probability_multipatch_functions.jl")
using .ExtinctionMultipatch

fecundity = 4
delta = 0.1
alpha1s = [0.6, 0.9]
p1s = [0.5, 0.6, 0.7, 0.8, 0.9, 1]

for alpha1 in alpha1s
    for p1 in p1s
        print("alpha1: $(alpha1), p1: $(p1)")
        beta1 = 1 - alpha1
        
        alpha2 = beta1
        beta2 = alpha1

        save_dir = "hierarchical/fecundity$(fecundity)/delta$(delta)/bimodal/alpha1_$(alpha1)__p1_$(p1)"

        minimize_extinction_probability(fecundity, delta, alpha1, beta1, alpha2, beta2, p1, save_dir, 100, true, true, 1e-8, 0.5, 75, 100)
    end
end

