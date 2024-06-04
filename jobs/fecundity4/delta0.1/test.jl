get_extinction_prob([0.1,0.3,0.5,0.7], [[1,2,3,4]], 0.1, 4, 2, 2, 4, 0.5, 4)
get_extinction_prob([0.125,0.375,0.625,0.875], [[1,2,3,4]], 0.1, 4, 2, 2, 4, 0.5, 4)
get_extinction_prob([0.125,0.375,0.625,0.875], [[1],[2],[3],[4]], 0.13, 4, 2, 2, 4, 0.5, 4)
get_extinction_prob([0.125,0.375,0.625,0.875], [[1],[2],[3],[4]], 0.13, 1, 1, 1, 1, 0.5, 4)
get_extinction_prob([0.125,0.375,0.625,0.875], [[1,2,3,4]], 0.13, 1, 1, 1, 1, 0.5, 4)
get_extinction_prob([0.125,0.375,0.625,0.875], [[1],[2],[3],[4]], 0.13, 1, 1, 1, 1, 0.5, 4)
get_extinction_prob([0.65,0.65,0.65,0.65], [[1,2],[3,4]], 0.1, 4, 2, 2, 4, 0.6, 4)
get_extinction_prob([0.7567992219390307,0.5593650801878238,0.7594179019715713,0.5567825159332134], [[1],[2],[3],[4]], 0.1, 4, 2, 2, 4, 0.7, 4)

e0 = get_extinction_prob([0.1000000000797513,0.1,0.9,0.9], [[1,2],[3],[4]], 0.1, 0.6, 0.4, 0.4, 0.6, 0.7, 4)
e1 = get_extinction_prob([0.1,0.1,0.9,0.9], [[1,2],[3],[4]], 0.1, 0.6, 0.4, 0.4, 0.6, 0.7, 4)
e2 = get_extinction_prob([0.1,0.1,0.9,0.9], [[1,3],[2],[4]], 0.1, 0.6, 0.4, 0.4, 0.6, 0.7, 4)
e3 = get_extinction_prob([0.1,0.1,0.9,0.9], [[1,3],[2,4]], 0.1, 0.6, 0.4, 0.4, 0.6, 0.7, 4)
e4 = get_extinction_prob([0.1,0.1,0.9,0.9], [[1], [2], [3], [4]], 0.1, 0.6, 0.4, 0.4, 0.6, 0.7, 4)
e5 = get_extinction_prob([0.1,0.9,0.9,0.9], [[1], [2], [3], [4]], 0.1, 0.6, 0.4, 0.4, 0.6, 0.7, 4)
e6 = get_extinction_prob([0.1,0.9,0.9,0.9], [[1,2], [3], [4]], 0.1, 0.6, 0.4, 0.4, 0.6, 0.7, 4)
e1 < e0
e2 < e1
e3 < e2
e4 < e3
e5 < e3
e6 < e3


alpha1 = beta1 = alpha2 = beta2 = 0.5

alpha1 = beta1 = alpha2 = beta2 = 1.0

alpha1 = beta1 = alpha2 = beta2 = 2.0

alpha1 = beta2 = 0.6
alpha2 = beta1 = 0.4

alpha1 = beta2 = 1
alpha2 = beta1 = 0.4

julia --project=. src/minimize_extinction_probability_multipatch.jl brute-force-minimize --fecundity 4 --delta 0.1 --alpha1 1.0 --beta1 0.05 --alpha2 0.05 --beta2 1.0 --p1 0.6 --save-dir fecundity4/delta0.1/extreme_mode__varying_mass_everywhere_else/beta1_0.05__p1_0.6 --use-educated-guess --population-size 100

julia --project=. src/minimize_extinction_probability_multipatch.jl brute-force-minimize --fecundity 4 --delta 0.1 --alpha1 1.0 --beta1 0.05 --alpha2 0.05 --beta2 1.0 --p1 0.6 --save-dir fecundity4/delta0.1/extreme_mode__varying_mass_everywhere_else/beta1_0.05__p1_0.6 --use-educated-guess

julia --compile=all --project=. src/minimize_extinction_probability_multipatch.jl brute-force-minimize --fecundity 4 --delta 0.1 --alpha1 1.0 --beta1 0.05 --alpha2 0.05 --beta2 1.0 --p1 0.6 --save-dir test8 --use-educated-guess --q0 0.5 --num-generations 75 --num-runs 100

julia --compiled-modules=yes --project=. src/minimize_extinction_probability_multipatch.jl brute-force-minimize --fecundity 4 --delta 0.1 --alpha1 1.0 --beta1 0.05 --alpha2 0.05 --beta2 1.0 --p1 0.6 --save-dir test8 --use-educated-guess --q0 0.5 --num-generations 75 --num-runs 100

julia --project=. src/minimize_extinction_probability_multipatch.jl brute-force-minimize --fecundity 4 --delta 0.1 --alpha1 1.0 --beta1 0.05 --alpha2 0.05 --beta2 1.0 --p1 0.6 --save-dir test8 --use-educated-guess --q0 0.5 --num-generations 75 --num-runs 100

inner(a,b) = a + b
inner(a,b,c,d) = a * b * c * d

function outer(a,b,args...)
    inner(a,b,args...)
end


include("src/extinction_probability_multipatch_functions.jl")
using .ExtinctionMultipatch

fecundity = 4
delta = 0.1
alpha1s = [0.6, 0.9]
p1s = [0.5, 0.6, 0.7, 0.8, 0.9, 1]

alpha1 = 0.6
p1 = 1.0

beta1 = 1 - alpha1
        
alpha2 = beta1
beta2 = alpha1

save_dir = "hierarchical/fecundity$(fecundity)/delta$(delta)/bimodal/alpha1_$(alpha1)__p1_$(p1)"

out = minimize_extinction_probability(fecundity, delta, alpha1, beta1, alpha2, beta2, p1, save_dir, population_size, true, true, 1e-8, 0.5, 75, 100)



fecundity = 4
delta = 0.1
alpha1 = beta2 = 0.6
alpha2 = beta1 = 0.4
p1 = 0.5
save_dir = "test"
use_educated_guess = true
analyze_centers = true
population_size = 100
num_generations = 75
num_runs = 100
tol = 1e-8
q0 = 0.5
num_generations = 75
num_runs = 100

partition = [[1,2,3,4]]

out = minimize_partition_extinction_probability(partition, fecundity, delta, alpha1, beta1, alpha2, beta2, p1, population_size, use_educated_guess, analyze_centers, tol, q0, num_generations, num_runs)

out = minimize_extinction_probability(fecundity, delta, alpha1, beta1, alpha2, beta2, p1, save_dir, population_size, use_educated_guess, analyze_centers, tol)


objective_function = make_objective_function(partition, fecundity, delta, alpha1, beta1, alpha2, beta2, p1, q0, num_generations, num_runs)

objective_function([0.1, 0.1, 0.1, 0.1])


centers = [0.1, 0.3, 0.7, 0.9]

get_extinction_probability(centers, partition, delta, alpha1, beta1, alpha2, beta2, p1, fecundity, q0, num_generations, num_runs)

using Distributions

approximate_extinction_probability2(probabilities_by_environment, q0, p1, num_generations)

probabilities_by_environment = [[0.12194767121906397, 0.8780523287809361, 0.5, 0.0, 0.0]]