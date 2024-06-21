include("src/extinction_probability_multipatch_functions.jl")
using .ExtinctionMultipatch



fecundity = 4
delta = 0.1
p1s = 0.5:1/40:1

for i in 0:4
    for p1 in p1s
        alpha1 = 20
        beta1 = 10
        
        alpha2 = beta1
        beta2 = alpha1

        save_dir = "repeated1/hierarchical/fecundity$(fecundity)/delta$(delta)/central_mode__varying_width/alpha1_20__beta1_10__p1_$(p1)_seed_$(i)"

        minimize_extinction_probability(fecundity, delta, alpha1, beta1, alpha2, beta2, p1, save_dir, 100, true, true, 1e-4, i, 1000, 0.05, 0.5, 75, 100)

    end
end



