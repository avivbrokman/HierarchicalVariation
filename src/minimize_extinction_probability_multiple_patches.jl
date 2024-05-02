# packages
using Intervals
using Distributions
using Polynomials
using Evolutionary
using LinearAlgebra 
using IterTools
using Parameters
using Combinatorics
using Random



##
# Function to print a partition for clarity
# function print_partition(partition)
#     println([collect(subset) for subset in partition])
# end

## carries global parameter values (from command line?)
# struct GlobalParams
#     delta::Float64
#     alpha1::Float64
#     beta1::Float64
#     alpha2::Float64
#     beta2::Float64
#     p1::Float64
#     fecundity::Int8
#     population_size::Int8
#     p_init_same_patch::Float64
#     bit_flipping_prob::Float64
# end

@with_kw struct GlobalParams
    delta::Float64
    alpha1::Float64
    beta1::Float64
    alpha2::Float64
    beta2::Float64
    p1::Float64
    fecundity::Int8
    population_size::Int8
    p_init_same_patch::Float64
    bit_flipping_prob::Float64
end

## utilities
better_binomial(n) = Int(round(binomial(n,2)))

## stuff to construct adjacency2partition
# Recursive function to generate all partitions
function partition_set(set, index, ans)
    if index == length(set) + 1
        # If we have considered all elements in the set, print the partition
        # print_partition(ans)
        return ans
    end

    # For each subset in the partition, add the current element to it and recall
    for i in 1:length(ans)
        push!(ans[i], set[index])
        partition_set(set, index + 1, ans)
        pop!(ans[i])
    end

    # Add the current element as a singleton subset and recall
    push!(ans, [set[index]])
    partition_set(set, index + 1, ans)
    pop!(ans)
end

# Function to initiate the partition process
function all_partitions(set)
    ans = []
    partition_set(set, 1, ans)
end

# Example usage with set [1, 2, 3, 4]
# a = all_partitions([1, 2, 3, 4])


# Recursive function to generate and collect all partitions
function partition_set(set, index, current_partition, collected_partitions)
    if index == length(set) + 1
        # If we have considered all elements in the set, add the partition to collected_partitions
        push!(collected_partitions, deepcopy(current_partition))
        return
    end

    # For each subset in the partition, add the current element to it and recall
    for subset in current_partition
        push!(subset, set[index])
        partition_set(set, index + 1, current_partition, collected_partitions)
        pop!(subset)
    end

    # Add the current element as a singleton subset and recall
    push!(current_partition, [set[index]])
    partition_set(set, index + 1, current_partition, collected_partitions)
    pop!(current_partition)
end

# Function to initiate the partition process and collect results
function all_partitions(set)
    collected_partitions = []
    partition_set(set, 1, [], collected_partitions)
    return collected_partitions
end

function partition2adjacency_single(partition, fecundity)
    
    # Initialize an adjacency matrix of size total_elements x total_elements
    adjacency = zeros(Int, fecundity, fecundity)

    # For each subset in the partition, mark edges between all pairs of vertices within the subset
    for subset in partition
        for i in 1:length(subset)
            for j in i+1:length(subset)
                adjacency[subset[i], subset[j]] = 1
                adjacency[subset[j], subset[i]] = 1  # Ensure the matrix is symmetric
            end
        end
    end

    return adjacency
end

function get_adjacency2partition(partitions, fecundity)
    return Dict(partition2adjacency_single(el, fecundity) => el for el in partitions)
end

## segmenting (0,1) by how many offspring survive in the patch for a given environmental value
function survival_interval(center, delta)
    left = max(center - delta, 0)
    right = min(center + delta, 1)
    return Interval{Open,Open}(left, right)
end

function get_survival_intervals(centers, delta)
    return [survival_interval(el, delta) for el in centers]
end


function segment_unit_interval_by_survival_overlap(intervals)
    endpoint_pairs = [[el.first, el.last] for el in intervals]
    endpoints = vcat(endpoint_pairs...)
    endpoints = [0; endpoints; 1]
    unique_endpoints = unique(endpoints)
    sorted_endpoints = sort(unique_endpoints)

    # for i in 1:length(sorted_endpoints) - 1
    #     print(i)
    #     print(sorted_endpoints[i])
    # end

    segments = [Interval{Open,Open}(sorted_endpoints[i], sorted_endpoints[i + 1]) for i in 1:length(sorted_endpoints) - 1]

    return segments
end


function get_segment_survival_counts(survival_intervals)

    segments = segment_unit_interval_by_survival_overlap(survival_intervals)
    segment_survival_by_individual = Intervals.find_intersections(segments, survival_intervals)
    segment_survival_counts = [length(el) for el in segment_survival_by_individual]
    return segments, segment_survival_counts
end

## Getting probabilities of survival within a patch for a given hierchical variable value
function get_segment_probability(segment, distribution)
    return cdf(distribution, segment.last) - cdf(distribution, segment.first)
end

function get_segment_probabilities(segments, distribution)
    return [get_segment_probability(el, distribution) for el in segments]
end

function get_probabilities_in_patch_given_H(segment_survival_counts, segment_probabilities, fecundity)
    
    # probabilities = zeros(Float64, fecundity + 1)
    
    probabilities = Dict(i => 0.0 for i in 0:fecundity)

    # Iterate over the keys and values simultaneously
    for (key, value) in zip(segment_survival_counts, segment_probabilities)
        probabilities[key] += value
    end
    
    return probabilities
end

## get probabilities of survival across patches and hierarchical variable values 
function get_probabilities_given_H(probabilities_vec, fecundity)
    
    # probabilities = Dict(i => 0.0 for i in 0:fecundity)
    probabilities = fill(0.0, fecundity)

    survival_iterator = (keys(el) for el in probabilities_vec)
    probability_iterator = (values(el) for el in probabilties_vec)

    cartesian_survival = product(survival_iterator)
    cartesian_probability = product(probability_iterator)

    cartesian_survival_sums = (sum(el) for el in cartesian_survival)
    cartesian_probability_products = (prod(el) for el in cartesian_probability)

    for (el_key, el_value) in zip(cartesian_survival_sums, cartesian_probability_products)
        probabilities[el_key] += el_value
    end

    return probabilities
end

function get_probabilities(probabilities1, probabilities2, p1)
    return p1 .* probabilities1 + (1-p1) .* probabilities2
end

## converting overall survival probabilities into extinction probability
function probabilities2extinction_coefficients!(probabilities)
    probabilities[2] -= 1
end 

function get_extinction_prob_from_coefficients(coefficients)
    # make polynomial
    coefficients_vector = [coefficients[i] for i in 0:maximum(keys(coefficients))]

    p = Polynomial(coefficients_vector)

    # get roots
    roots_p = roots(p)

    # filter out complex roots
    real_roots = filter(root -> isreal(root), roots_p)
    real_roots = [real(el) for el in real_roots]

    # Filter out negative roots
    positive_roots = filter(root -> root >= 0, real_roots)

    # get lower root
    extinction_prob = minimum(positive_roots)

    return extinction_prob
end

function warshall_algorithm(A)
    n = size(A, 1)
    T = copy(A)  # Create a copy of the adjacency matrix

    for k in 1:n
        for i in 1:n
            for j in 1:n
                T[i, j] = T[i, j] | (T[i, k] & T[k, j])
            end
        end
    end
    return T
end

function get_adjacency_lower_triangle2closure(fecundity)
    all_vectors = collect(product(0:1, fecundity - 1))
    adjacencies = [reconstruct_adjacency(el) for el in all_vectors]
    closures = [warshall_algorithm(el) for el in adjacencies]
    flattened_closures = [flatten_adjacency(el) for el in closures]

    
    
    
    

## workhorse function
function get_extinction_prob(centers, adjacency, delta, alpha1, beta1, alpha2, beta2, p1, fecundity, adjacency2partition)
    # get Distributions
    dist1 = Beta(alpha1, beta1)
    dist2 = Beta(alpha2, beta2)
    
    # get intervals
    intervals = get_survival_intervals(centers, delta)

    # adjacency closure
    # adjacency_closure = sum(adjacency^i for i in 1:fecundity)
    adjacency_closure = warshall_algorithm(adjacency)
    
    # get clusters from closure
    patch_membership = adjacency2partition[adjacency_closure]

    patch_probabilities_given_H1_by_patch = fill(0,length(patch_membership))
    patch_probabilities_given_H2_by_patch = fill(0,length(patch_membership))

    for (i, members) in enumerate(patch_membership)
        patch_intervals = [intervals[i] for i in members]
        segments, segment_survival_counts = get_segment_survival_counts(patch_intervals)
        
        segment_probabilities_given_H1 = get_segment_probabilities(segments, dist1)
        segment_probabilities_given_H2 = get_segment_probabilities(segments, dist2)

        patch_probabilities_given_H1_by_patch[i] = get_probabilities_in_patch_given_H(segment_survival_counts, segment_probabilities_given_H1, fecundity)
        patch_probabilities_given_H2_by_patch[i] = get_probabilities_in_patch_given_H(segment_survival_counts, segment_probabilities_given_H2, fecundity)
    end

    probabilities_given_H1 = get_probabilities_given_H(patch_probabilities_given_H1_by_patch, fecundity)
    probabilities_given_H2 = get_probabilities_given_H(patch_probabilities_given_H2_by_patch, fecundity)

    probabilities = get_probabilities(probabilities_given_H1, probabilities_given_H2, p1)
    
    probabilities2extinction_coefficients!(probabilities)

    extinction_prob = get_extinction_prob_from_coefficients(coefficients)

    return extinction_prob
end


function get_extinction_prob2(centers, partition, delta, alpha1, beta1, alpha2, beta2, p1, fecundity)
    # get Distributions
    dist1 = Beta(alpha1, beta1)
    dist2 = Beta(alpha2, beta2)
    
    # get intervals
    intervals = get_survival_intervals(centers, delta)

    patch_probabilities_given_H1_by_patch = fill(0,length(partition))
    patch_probabilities_given_H2_by_patch = fill(0,length(partition))

    for (i, members) in enumerate(partition)
        patch_intervals = [intervals[i] for i in members]
        segments, segment_survival_counts = get_segment_survival_counts(patch_intervals)
        
        segment_probabilities_given_H1 = get_segment_probabilities(segments, dist1)
        segment_probabilities_given_H2 = get_segment_probabilities(segments, dist2)

        patch_probabilities_given_H1_by_patch[i] = get_probabilities_in_patch_given_H(segment_survival_counts, segment_probabilities_given_H1, fecundity)
        patch_probabilities_given_H2_by_patch[i] = get_probabilities_in_patch_given_H(segment_survival_counts, segment_probabilities_given_H2, fecundity)
    end

    probabilities_given_H1 = get_probabilities_given_H(patch_probabilities_given_H1_by_patch, fecundity)
    probabilities_given_H2 = get_probabilities_given_H(patch_probabilities_given_H2_by_patch, fecundity)

    probabilities = get_probabilities(probabilities_given_H1, probabilities_given_H2, p1)
    
    probabilities2extinction_coefficients!(probabilities)

    extinction_prob = get_extinction_prob_from_coefficients(coefficients)

    return extinction_prob
end

## using a closure as a workaround because DE's objective function must take a single vector as input

function flatten_adjacency(matrix)
    nrows, ncols = size(matrix)
    lower_triangle = []
    for row in 2:nrows
        for col in 1:row-1
            push!(lower_triangle, matrix[row, col])
        end
    end
    return lower_triangle
end

function reconstruct_adjacency(vector)
    # Calculate the number of rows in the original matrix
    # n * (n - 1) / 2 = length of vector
    n = Int(round((1 + sqrt(1 + 8*length(vector))) / 2))

    # Initialize an nxn matrix of zeros
    matrix = zeros(Int, n, n)

    # Fill the lower and upper triangles from the vector
    index = 1
    for row in 2:n
        for col in 1:(row - 1)
            matrix[row, col] = vector[index]
            matrix[col, row] = vector[index]
            index += 1
        end
    end
    return matrix
end

function make_objective_function(params::GlobalParams, adjacency2partition)
    function objective_function(vars)
        centers = vars[1:params.fecundity]
        adjacency_vars = vars[params.fecundity + 1:end]
        adjacency = reconstruct_adjacency(adjacency_vars)
        extinction_prob = get_extinction_prob(centers, adjacency, params.delta, params.alpha1, params.beta1, params.alpha2, params.beta2, params.p1, params.fecundity, adjacency2partition)
        return extinction_prob
    end

    return objective_function
end

## Customizing optimization algorithm to handle a mix of continuous and binary parameters
function initialize_population(fecundity, population_size, p_init_same_patch)
    # Prepare initial population
    num_patch_indices = better_binomial(fecundity)
    # patch_indices = fecundity .+ 1:num_patch_indices

    standard_uniform_dist = Uniform(0,1)
    bernoulli_dist = Bernoulli(p_init_same_patch)

    initial_population = [[rand(standard_uniform_dist, fecundity); rand(bernoulli_dist, num_patch_indices)] for _ in 1:population_size]

    return initial_population
end

function initialize_population2(fecundity, population_size)
    # Prepare initial population
    all_partitions = collect(partitions(1:fecundity))

    standard_uniform_dist = Uniform(0,1)

    initial_population = [[rand(standard_uniform_dist, fecundity); rand(all_partitions)] for _ in 1:population_size]

    return initial_population
end


# Custom mutation for DE
function createCustomDifferentiation(fecundity, bit_flipping_prob)
    function customDifferentiation(recombinant::T, mutators::AbstractVector{T}; F::Real = 1.0) where {T <: AbstractVector}
       
        m = length(mutators)
        @assert m % 2 == 0 "Must be an even number of target mutators"
        
        for i in 1:2:m
            recombinant[1:fecundity] .+= F .* (mutators[i][1:fecundity] .- mutators[i + 1][1:fecundity])
        end

        # Handle binary parameters
        for i in binary_indices
            recombinant[i] = rand() < bit_flipping_prob ? 1 - recombinant[i] : recombinant[i]
        end

        return recombinant
    end
    return customDifferentiation
end

function createCustomDifferentiation2(fecundity, bit_flipping_prob)
    function customDifferentiation2(recombinant::T, mutators::AbstractVector{T}; F::Real = 1.0) where {T <: AbstractVector}
       
        m = length(mutators)
        @assert m % 2 == 0 "Must be an even number of target mutators"
        
        for i in 1:2:m
            recombinant[1:fecundity] .+= F .* (mutators[i][1:fecundity] .- mutators[i + 1][1:fecundity])
        end

        # Handle partition stuff
        all_partitions = collect(Combinatorics.partitions(elements))

        recombinant[end] = rand() < partition_mutation_rate ? rand(all_partitions) : recombinant[end]
        end

        return recombinant
    end
    return customDifferentiation2
end
## Run the DE algorithm
# result = Evolutionary.optimize(fitness, initial_population, DE(customMutation; cr=0.8), options)




# running
params = GlobalParams(delta = 0.1, alpha1 = 0.7, beta1 = 0.3, alpha2 = 0.3, beta2 = 0.7, p1 = 0.5, fecundity = 4, population_size = 100, p_init_same_patch = 0.2, bit_flipping_prob = 0.2)

partitions = all_partitions(1:params.fecundity)
adjacency2partition = get_adjacency2partition(partitions, params.fecundity)

objective_function = make_objective_function(params, adjacency2partition)

initial_population = initialize_population(params.fecundity, params.population_size, params.p_init_same_patch)
mixed_mutation = createCustomDifferentiation(params.fecundity, params.bit_flipping_prob)

num_pars = params.fecundity + better_binomial(params.fecundity)
lower_constraint = fill(0.0, num_pars)
upper_constraint = fill(1.0, num_pars)
constraints = BoxConstraints(lower_constraint, upper_constraint)
de_algorithm = DE(populationSize = params.population_size, recombination = mixed_mutation)
result = Evolutionary.optimize(objective_function, initial_population, de_algorithm, Evolutionary.Options())

# output = Evolutionary.optimize(f = objective_function, population = initial_population, constraints = constraints, method = de_algorithm)
    
print(output)











# #
# A = [0 0 0 0; 1 0 0 0; 0 1 0 0; 0 0 1 0]  
# B = [0 0 0 0; 1 0 0 0; 1 0 0 0; 1 0 0 0]  
# C = [0 0 0 0; 1 0 0 0; 0 0 0 0; 0 0 1 0]  

# A = A + transpose(A)
# B = B + transpose(B)
# C = C + transpose(C)

# n = size(A, 1)
# function compute_transitive_closure!(A)
#     for k in 1:n
#         for i in (k+1):n  # Start from k+1 to skip diagonal and upper triangle
#             for j in 1:(i-1)  # Iterate up to i-1 to stay within the lower triangle
#                 if A[i, k] == 1 && A[k, j] == 1
#                     A[i, j] = 1  # Only update the lower triangle entry
#                 end
#             end
#         end
#     end
#     return A
# end

# function warshall(A)
#     n = size(A,1)
#     for k in 1:n
#         for i in 1:n
#             for j in 1:n
#                 A[i, j] = Bool(A[i, j]) || (Bool(A[i, k]) && Bool(A[k, j]))
#             end
#         end
#     end
#     return A
# end

# compute_transitive_closure!(A)
# compute_transitive_closure!(B)

# C = B + transpose(B)
# warshall(A)
# warshall(B)
# warshall(C)