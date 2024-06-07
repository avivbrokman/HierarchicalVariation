include("src/extinction_probability_multipatch_functions.jl")
using .ExtinctionMultipatch

fecundity = 4
delta = 0.1
alpha1s = [0.6, 0.9]
p1s = 0.5:1/40:1

for alpha1 in alpha1s
    for p1 in p1s
        print("alpha1: $(alpha1), p1: $(p1)")
        beta1 = 1 - alpha1
        
        alpha2 = beta1
        beta2 = alpha1

        save_dir = "finer2/hierarchical/fecundity$(fecundity)/delta$(delta)/bimodal/alpha1_$(alpha1)__p1_$(p1)"

        minimize_extinction_probability(fecundity, delta, alpha1, beta1, alpha2, beta2, p1, save_dir, 100, true, true, 1e-8, 0, 0.5, 75, 100)
    end
end









fecundity = 4
delta = 0.1
pairs_ = [(4,2), (7,3), (20,10), (23,7)] 
p1s = 0.5:1/40:1

for pair in pairs_
    for p1 in p1s
        alpha1 = pair[1]
        beta1 = pair[2]
        
        alpha2 = beta1
        beta2 = alpha1

        save_dir = "finer2/hierarchical/fecundity$(fecundity)/delta$(delta)/central_mode__varying_width/alpha1_$(alpha1)__beta1_$(beta1)__p1_$(p1)"

        minimize_extinction_probability(fecundity, delta, alpha1, beta1, alpha2, beta2, p1, save_dir, 100, true, true, 1e-8, 0, 0.5, 75, 100)
    end
end





fecundity = 4
delta = 0.1
beta1s = [0.05, 0.95]
p1s = 0.5:1/40:1

for beta1 in beta1s
    for p1 in p1s
        alpha1 = 1
        
        alpha2 = beta1
        beta2 = alpha1

        save_dir = "finer2/hierarchical/fecundity$(fecundity)/delta$(delta)/extreme_mode__varying_mass_everywhere_else/beta1_$(beta1)__p1_$(p1)"

        minimize_extinction_probability(fecundity, delta, alpha1, beta1, alpha2, beta2, p1, save_dir, 100, true, true, 1e-8, 0, 0.5, 75, 100)
    end
end







fecundity = 4
delta = 0.1
# alpha1s = [1.2, 1.8, 2, 5]
alpha1_p1_combs = Dict(1.2 => [], 1.8 => 0.7:1/40:1, 2 => 0.7:1/40:1, 3 => 0.6:1/40:1, 4 => 0.5:1/40:1, 5 => 0.5:1/40:1, 7.5 => 0.5:1/40:1, 10 => 0.5:1/40:1)

for (alpha1, p1s) in pairs(alpha1_p1_combs)
    for p1 in p1s
        beta1 = 1
        
        alpha2 = beta1
        beta2 = alpha1

        save_dir = "finer2/hierarchical/fecundity$(fecundity)/delta$(delta)/extreme_mode__varying_slope/alpha1_$(alpha1)__p1_$(p1)"

        minimize_extinction_probability(fecundity, delta, alpha1, beta1, alpha2, beta2, p1, save_dir, 100, true, true, 1e-8, 0, 0.5, 75, 100)
    end
end



