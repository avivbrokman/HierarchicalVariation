include("src/extinction_probability_multipatch_functions.jl")
using .ExtinctionMultipatch

fecundity = 4
delta = 0.1
pairs = [(4,2), (7,3), (20,10), (23,7)] 
p1s = [0.5, 0.6, 0.7, 0.8, 0.9, 1]

for pair in pairs
    for p1 in p1s
        alpha1 = pair[1]
        beta1 = pair[2]
        
        alpha2 = beta1
        beta2 = alpha1

        save_dir = "hierarchical/fecundity$(fecundity)/delta$(delta)/central_mode__varying_width/alpha1_$(alpha1)__beta1_$(beta1)__p1_$(p1)"

        minimize_extinction_probability(fecundity, delta, alpha1, beta1, alpha2, beta2, p1, save_dir, 100, true, true, 1e-8, 0.5, 75, 100)
    end
end

