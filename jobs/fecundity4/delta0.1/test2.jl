using Revise
using Optimization, Optim
using OptimizationOptimJL
using ForwardDiff

include("/Users/avivbrokman/Documents/Kentucky/Grad School/ms_project/branching1/src/extinction_probability_multipatch_functions.jl")
using .ExtinctionMultipatch


out = minimize_extinction_probability(4, 0.1, 7, 3, 3, 7, 0.8, "test1", 100, true, true, 1e-4, 0, 1000, 0.05, 0.5, 75, 100)

print(3)

fecundity = 4
delta = 0.1
alpha1 = beta2 = 20
alpha2 = beta1 = 10
p1 = 0.5



partition1 = [[1,2,3,4]]
partition2 = [[1,2,3], [4]]

probs1 = get_probabilities_by_environment(centers1, partition1, delta, alpha1, beta1, alpha2, beta2, fecundity)

probs2 = get_probabilities_by_environment(centers2, partition2, delta, alpha1, beta1, alpha2, beta2, fecundity)

centers1_05 = [0.41282216352693996,0.7579786080936614,0.23362895555371976,0.5810405164870959]
centers2_05 = [0.4712654753170071,0.2717734431260377,0.6732609981480132,0.6830986482734478]

centers1_075 = [0.2958713367765998,0.7308261162716772,0.4804577342334033,0.62435873015993]
centers2_075 = [0.3109673841788403,0.7057630937627754,0.5058301984214542,0.6514044883680641]

centers1 = centers1_075
centers2 = centers2_075

ext1_1 = get_extinction_probability(centers1, partition1, delta, alpha1, beta1, alpha1, beta1, p1, fecundity)
ext1_2 = get_extinction_probability(centers1, partition1, delta, alpha2, beta2, alpha2, beta2, p1, fecundity)

ext2_1 = get_extinction_probability(centers2, partition2, delta, alpha1, beta1, alpha1, beta1, p1, fecundity)
ext2_2 = get_extinction_probability(centers2, partition2, delta, alpha2, beta2, alpha2, beta2, p1, fecundity)

sign(ext1_1 - ext2_1) == sign(ext1_2 - ext2_2)



fecundity = 4
delta = 0.1
alpha1 = beta2 = 20
alpha2 = beta1 = 10
p1 = 0.5

q0 = 0.5
num_generations = 75
num_runs = 500

seed = 0

partition = partition1

centers1_08 = [0.6509420062656308,0.5038538534526574,0.7729138823822734,0.2767975398899015]

objective = centers -> get_extinction_probability(centers, partition, delta, alpha1, beta1, alpha2, beta2, p1, fecundity, q0, num_generations, num_runs, seed)
guess = centers1_08
solution = optimize(objective, fill(delta, fecundity), fill(1 - delta, fecundity), guess, Fminbox(LBFGS()), Optim.Options(show_trace=false))
best = solution.minimizer
