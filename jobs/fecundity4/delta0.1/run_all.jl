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





fecundity = 4
delta = 0.1
beta1s = [0.05, 0.95]
p1s = [0.5, 0.6, 0.7, 0.8, 0.9, 1]

for beta1 in beta1s
    for p1 in p1s
        alpha1 = 1
        
        alpha2 = beta1
        beta2 = alpha1

        save_dir = "hierarchical/fecundity$(fecundity)/delta$(delta)/extreme_mode__varying_mass_everywhere_else/beta1_$(beta1)__p1_$(p1)"

        minimize_extinction_probability(fecundity, delta, alpha1, beta1, alpha2, beta2, p1, save_dir, 100, true, true, 1e-8, 0.5, 75, 100)
    end
end







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

