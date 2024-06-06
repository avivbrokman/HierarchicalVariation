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

out = minimize_extinction_probability(4, 0.1, 0.9, 0.1, 0.1, 0.9, 0.825, "test", 100, true, true, 1e-8, 1, 0.5, 75, 100)

using Random
Random.seed!(0)

partitions = get_partitions(4)

centers = [0.1, 0.9, 0.9, 0.9]
delta = 0.1
alpha1 = beta2 = 0.9
alpha2 = beta1 = 0.1
p1 = 0.806
fecundity = 4
q0 = 0.5
num_generations = 75
num_runs = 100
partition4 = [[1], [2], [3], [4]]
partition3 = [[1,2], [3], [4]]

threes = []
fours = []
for el in 1:100
    out3 = get_extinction_probability(centers, partition3, delta, alpha1, beta1, alpha2, beta2, p1, fecundity, q0, num_generations, num_runs)
    out4 = get_extinction_probability(centers, partition4, delta, alpha1, beta1, alpha2, beta2, p1, fecundity, q0, num_generations, num_runs)
    push!(threes, out3)
    push!(fours, out4)
end

# figure out all random seed stuff
# figure out if you really want to run this with the averaging, or if you want to get raw run output from the get_extinction_probability, so you can compare those, and reduce the number of times this function must be run.
 
function col_argmax(matrix::Matrix)
    cartesian_argmax = argmax(matrix, dims = 1)
    return [el[1] for el in cartesian_argmax]
end

function col_argmax(data::Vector{Vector})
    matrix = hcat(data...)
    matrix = transpose(matrix)
    return col_argmax(matrix)
end



function get_extinction_probability_runs(centers, partitions, delta, alpha1, beta1, alpha2, beta2, p1, fecundity, q0, num_generations, num_runs)
    results = Dict()
    for el in partitions
        Random.seed!(1)
        results[el] = Vector{Float64}()
        for _ in 1:100
            push!(results[el], get_extinction_probability(centers, el, delta, alpha1, beta1, alpha2, beta2, p1, fecundity, q0, num_generations, num_runs))
        end
    end
    
    return results
end

function check_optimal_partition(centers, partitions, delta, alpha1, beta1, alpha2, beta2, p1, fecundity, q0, num_generations, num_runs)

    results = get_extinction_probability_runs(centers, partitions, delta, alpha1, beta1, alpha2, beta2, p1, fecundity, q0, num_generations, num_runs)

    best_index_by_run = col_argmax(values(results))
    best_index = mode(best_index_by_run)
    best_partition = keys(results)[best_index]

    return best_partition
end

d = Dict("a" => [1,2,3], "b" => [4,5,6])









function approximate_extinction_probability_med(probabilities_by_environment, q0, p1, num_generations, num_runs)
    # print("approximate 33333")

    extinction_probability_runs = [approximate_extinction_probability(probabilities_by_environment, q0, p1, num_generations) for el in 1:num_runs]

    extinction_probability = median(extinction_probability_runs)

    return extinction_probability
end


function get_extinction_probability_med(centers, partition, delta, alpha1, beta1, alpha2, beta2, p1, fecundity, q0, num_generations, num_runs) 
    # get Distributions
    dist1 = Beta(alpha1, beta1)
    dist2 = Beta(alpha2, beta2)
    
    # get intervals
    intervals_ = get_survival_intervals(centers, delta)

    patch_probabilities_given_H1_by_patch = Vector{Dict{Int64, Float64}}()
    patch_probabilities_given_H2_by_patch = Vector{Dict{Int64, Float64}}()

    for (i, members) in enumerate(partition)
        patch_intervals = [intervals_[i] for i in members]
        segments, segment_survival_counts = get_segment_survival_counts(patch_intervals)
        
        segment_probabilities_given_H1 = get_segment_probabilities(segments, dist1)
        segment_probabilities_given_H2 = get_segment_probabilities(segments, dist2)

        push!(patch_probabilities_given_H1_by_patch, get_probabilities_in_patch_given_H(segment_survival_counts, segment_probabilities_given_H1, fecundity))
        push!(patch_probabilities_given_H2_by_patch, get_probabilities_in_patch_given_H(segment_survival_counts, segment_probabilities_given_H2, fecundity))
    end

    probabilities_given_H1 = get_probabilities_given_H(patch_probabilities_given_H1_by_patch, fecundity)
    probabilities_given_H2 = get_probabilities_given_H(patch_probabilities_given_H2_by_patch, fecundity)

    probabilities_by_environment = [probabilities_given_H1, probabilities_given_H2]

    extinction_probability = approximate_extinction_probability_med(probabilities_by_environment, q0, p1, num_generations, num_runs)

    return extinction_probability
end

Random.seed!(0)
medians = [get_extinction_probability_med(centers, partition3, delta, alpha1, beta1, alpha2, beta2, p1, fecundity, q0, num_generations, num_runs) for _ in 1:100]


Random.seed!(0)
means = [get_extinction_probability(centers, partition3, delta, alpha1, beta1, alpha2, beta2, p1, fecundity, q0, num_generations, num_runs) for el in 1:100]

std(means)
std(medians)


function reproduce(probabilities::Vector{Float64})
    return StatsBase.sample([0,1,2,3,4], Weights(probabilities))
end

function reproduce(probabilities::Vector{Vector{Float64}}, environment::Int)
    probabilities = probabilities[environment]
    return reproduce(probabilities)
end

function reproduce(population_size::Int, probabilities::Vector{Vector{Float64}}, p1::Number)
    environment = StatsBase.sample([1,2], Weights([p1, 1-p1]))
    
    return sum(reproduce(probabilities, environment) for _ in 1:population_size)
end

function goes_extinct(probabilities::Vector{Vector{Float64}}, p1::Number, max_generations::Int, max_population_size::Int)
    population_size = 1
    generation = 0
    while population_size < max_population_size && generation < max_generations
        population_size = reproduce(population_size, probabilities, p1)
        generation += 1
        if population_size == 0
            return true
        end
    end
    return false
end

function estimate_extinction_probability(probabilities::Vector{Vector{Float64}}, p1::Number, max_generations::Int, max_population_size::Int, num_runs::Int)
    results = [goes_extinct(probabilities, p1, max_generations, max_population_size) for _ in 1:num_runs]
    return mean(results)
end

function estimate_extinction_probability(centers, partition, delta, alpha1, beta1, alpha2, beta2, p1::Number, max_generations::Int, max_population_size::Int, num_runs::Int)
    # get Distributions
    dist1 = Beta(alpha1, beta1)
    dist2 = Beta(alpha2, beta2)
    
    # get intervals
    intervals_ = get_survival_intervals(centers, delta)

    patch_probabilities_given_H1_by_patch = Vector{Dict{Int64, Float64}}()
    patch_probabilities_given_H2_by_patch = Vector{Dict{Int64, Float64}}()

    for (i, members) in enumerate(partition)
        patch_intervals = [intervals_[i] for i in members]
        segments, segment_survival_counts = get_segment_survival_counts(patch_intervals)
        
        segment_probabilities_given_H1 = get_segment_probabilities(segments, dist1)
        segment_probabilities_given_H2 = get_segment_probabilities(segments, dist2)

        push!(patch_probabilities_given_H1_by_patch, get_probabilities_in_patch_given_H(segment_survival_counts, segment_probabilities_given_H1, fecundity))
        push!(patch_probabilities_given_H2_by_patch, get_probabilities_in_patch_given_H(segment_survival_counts, segment_probabilities_given_H2, fecundity))
    end

    probabilities_given_H1 = get_probabilities_given_H(patch_probabilities_given_H1_by_patch, fecundity)
    probabilities_given_H2 = get_probabilities_given_H(patch_probabilities_given_H2_by_patch, fecundity)

    probabilities_by_environment = [probabilities_given_H1, probabilities_given_H2]

    extinction_probability = estimate_extinction_probability(probabilities_by_environment, p1, max_generations, max_population_size, num_runs)
    
    return extinction_probability

end

max_generations = 100
max_population_size = 100
num_runs = 1000000

out = estimate_extinction_probability(centers, partition3, delta, alpha1, beta1, alpha2, beta2, p1, max_generations, max_population_size, num_runs)