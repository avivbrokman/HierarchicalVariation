# packages
using Intervals
using Distributions
using Polynomials
using CustomEvolutionary1
using LinearAlgebra 
using IterTools
using Parameters
using Combinatorics
using Random
using ArgParse


##
# Function to print a partition for clarity
# function print_partition(partition)
#     println([collect(subset) for subset in partition])
# end




## carries global parameter values (from command line?)
@with_kw struct GlobalParams
    delta::Float64
    alpha1::Float64
    beta1::Float64
    alpha2::Float64
    beta2::Float64
    p1::Float64
    fecundity::Int64
    population_size::Int64
    partition_mutation_rate::Float64
end

## Individual
# struct Individual
#     centers::Vector{Float64}
#     partition::Vector{Vector{Int}}
# end

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
    probabilities = fill(0.0, fecundity + 1)

    survival_iterator = collect(keys(el) for el in probabilities_vec)
    probability_iterator = collect(values(el) for el in probabilities_vec)

    cartesian_survival = collect(product(survival_iterator...))
    cartesian_probability = collect(product(probability_iterator...))

    cartesian_survival_sums = collect(sum(el) for el in cartesian_survival)
    cartesian_probability_products = collect(prod(el) for el in cartesian_probability)

    for (el_key, el_value) in zip(cartesian_survival_sums, cartesian_probability_products)
        if el_value > 0.0
            probabilities[el_key + 1] += el_value
        end
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
    # coefficients_vector = [coefficients[i] for i in 0:maximum(keys(coefficients))]

    # p = Polynomial(coefficients_vector)
    p = Polynomial(coefficients)

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


## workhorse function
function get_extinction_prob(centers, partition, delta, alpha1, beta1, alpha2, beta2, p1, fecundity)
    # get Distributions
    dist1 = Beta(alpha1, beta1)
    dist2 = Beta(alpha2, beta2)
    
    # get intervals
    intervals = get_survival_intervals(centers, delta)

    patch_probabilities_given_H1_by_patch = Vector{Dict{Int64, Float64}}()
    patch_probabilities_given_H2_by_patch = Vector{Dict{Int64, Float64}}()

    for (i, members) in enumerate(partition)
        patch_intervals = [intervals[i] for i in members]
        segments, segment_survival_counts = get_segment_survival_counts(patch_intervals)
        
        segment_probabilities_given_H1 = get_segment_probabilities(segments, dist1)
        segment_probabilities_given_H2 = get_segment_probabilities(segments, dist2)

        push!(patch_probabilities_given_H1_by_patch, get_probabilities_in_patch_given_H(segment_survival_counts, segment_probabilities_given_H1, fecundity))
        push!(patch_probabilities_given_H2_by_patch, get_probabilities_in_patch_given_H(segment_survival_counts, segment_probabilities_given_H2, fecundity))
    end

    probabilities_given_H1 = get_probabilities_given_H(patch_probabilities_given_H1_by_patch, fecundity)
    probabilities_given_H2 = get_probabilities_given_H(patch_probabilities_given_H2_by_patch, fecundity)

    probabilities = get_probabilities(probabilities_given_H1, probabilities_given_H2, p1)
    
    probabilities2extinction_coefficients!(probabilities)

    extinction_prob = get_extinction_prob_from_coefficients(probabilities)

    return extinction_prob
end

## using a closure as a workaround because DE's objective function must take a single vector as input
function make_objective_function(params::GlobalParams, idx2partition)
    function objective_function(vars)
        centers = vars[1:params.fecundity]
        partition_idx = vars[end]
        partition = idx2partition[partition_idx]
        extinction_prob = get_extinction_prob(centers, partition, params.delta, params.alpha1, params.beta1, params.alpha2, params.beta2, params.p1, params.fecundity)
        return extinction_prob
    end

    return objective_function
end

function get_idx2partition(fecundity)
    all_partitions = collect(partitions(1:fecundity))
    idx2partition = Dict(float(i) => el for (i, el) in enumerate(all_partitions))
    return idx2partition
end


## Customizing optimization algorithm to handle a mix of continuous and binary parameters
function initialize_population(population_size, fecundity, idx2partition)
    # Prepare initial population
    # all_partitions = collect(partitions(1:fecundity))

    standard_uniform_dist = Uniform(0,1)

    # initial_population = [[rand(standard_uniform_dist, fecundity); [rand(all_partitions)]] for _ in 1:population_size]
    initial_population = [[rand(standard_uniform_dist, fecundity); rand(keys(idx2partition))] for _ in 1:population_size]

    return initial_population
end


# function create_custom_differentiation(fecundity, partition_mutation_rate, idx2partition)
#     # function customDifferentiation(recombinant::T, mutators::AbstractVector{T}; F::Real = 1.0) where {T <: AbstractVector}
#     function custom_differentiation(recombinant::T, mutators::AbstractVector{T}; F = 1.0, rng = default_rng()) where {T <: AbstractVector}

       
#         m = length(mutators)
#         @assert m % 2 == 0 "Must be an even number of target mutators"
        
#         for i in 1:2:m
#             recombinant[1:fecundity] .+= F .* (mutators[i][1:fecundity] .- mutators[i + 1][1:fecundity])
#         end

#         # Handle partition stuff
#         # all_partitions = collect(Combinatorics.partitions(elements))

#         recombinant[end] = rand() < partition_mutation_rate ? rand(keys(idx2partition)) : recombinant[end]
#         return recombinant
#     end
#     return custom_differentiation
# end

# function original_update_state!(args...; kwargs...)
#     Evolutionary.update_state!(args...; kwargs...)
# end


# struct CustomDEWrapper
#     de: Evolutionary.DE  # The original DE method
#     differentiation: Function  # Your custom differentiation function
# end




# function Evolutionary.custom_update_state!(objfun, constraints, state, population::AbstractVector{IT}, method::CustomDEWrapper, options, itr) where {IT}

#     # setup
#     Np = method.de.populationSize
#     n = method.de.n
#     F = method.de.F
#     rng = options.rng

#     offspring = Array{IT}(undef, Np)

#     # select base vectors
#     bases = method.de.selection(state.fitness, Np)

#     # select target vectors
#     for (i,b) in enumerate(bases)
#         # mutation
#         base = population[b]
#         offspring[i] = copy(base)
#         # println("$i => base:", offspring[i])

#         targets = randexcl(rng, 1:Np, [i], 2*n)
#         offspring[i] = method.differentiation(offspring[i], @view population[targets]; F=F)
#         # println("$i => mutated:", offspring[i], ", targets:", targets)

#         # recombination
#         offspring[i], _ = method.de.recombination(offspring[i], base, rng=rng)
#         # println("$i => recombined:", offspring[i])
#     end

#     # Create new generation
#     fitidx = 0
#     minfit = Inf
#     for i in 1:Np
#         o = apply!(constraints, offspring[i])
#         v = value(objfun, o) + penalty(constraints, o)
#         if (v <= state.fitness[i])
#             population[i] = o
#             state.fitness[i] = v
#             if v < minfit
#                 minfit = v
#                 fitidx = i
#             end
#         end
#     end

#     # set best individual
#     if fitidx > 0
#         state.fittest = population[fitidx]
#     end

#     return false
# end




# function create_custom_constraints(fecundity)
#     function custom_constraint(x)
#         # Apply scalar constraints to the first four parameters
#         for i in 1:fecundity
#             x[i] = clamp(x[i], 0, 1)  
#         end
#         return x
#     end
#     return custom_constraint
# end

# struct CustomConstraints <: Evolutionary.AbstractConstraints
#     fecundity::Int
#     lower::Float64
#     upper::Float64
# end


## Run the DE algorithm
# result = Evolutionary.optimize(fitness, initial_population, DE(customMutation; cr=0.8), options)


function main()
    s = ArgParseSettings()  # Create a settings object

    @add_arg_table s begin
        "input"
            help = "input file"
            required = true
            arg_type = String
        "--output", "-o"
            help = "output file"
            default = "out.txt"  # Default value if not specified
            arg_type = String
        "--verbose", "-v"
            help = "print more verbose output"
            action = :store_true  # This makes it a boolean flag
    end

    parsed_args = parse_args(s)  # This function parses the command line arguments based on the specified settings
    println("Input file: ", parsed_args["input"])
    println("Output file: ", parsed_args["output"])
    if parsed_args["verbose"]
        println("Verbose mode is on.")
    else
        println("Verbose mode is off.")
    end
end

# This conditional ensures that the script runs main only if it is not being included as a module
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end



# running
params = GlobalParams(delta = 0.1, alpha1 = 0.7, beta1 = 0.3, alpha2 = 0.3, beta2 = 0.7, p1 = 0.5, fecundity = 4, population_size = 100, partition_mutation_rate = 0.2)

idx2partition = get_idx2partition(params.fecundity)

objective_function = make_objective_function(params, idx2partition)

initial_population = initialize_population(params.population_size, params.fecundity, idx2partition)
custom_differentiation = create_custom_differentiation(params.fecundity, params.partition_mutation_rate, idx2partition)
# update_state! = create_custom_update_state(mixed_mutation)



lower_constraint = fill(0.0, params.fecundity)
upper_constraint = fill(1.0, params.fecundity)

lower_constraint = [lower_constraint; float(minimum(keys(idx2partition)))]
upper_constraint = [upper_constraint; float(maximum(keys(idx2partition)))]

constraints = BoxConstraints(lower_constraint, upper_constraint)
# constraints = CustomConstraints(4, 0.0, 1.0)


# de_algorithm = CustomDE(populationSize = params.population_size)
# de_algorithm = CustomDEWrapper(populationSize = params.population_size)
de_algorithm = DE(populationSize = params.population_size, differentiation = custom_differentiation)
# update_state!(objfun, constraints, state, population, de_algorithm, options, itr) = custom_update_state!(objfun, constraints, state, population, de_algorithm, options, itr)

# result = Evolutionary.optimize(objective_function, initial_population, de_algorithm, Evolutionary.Options())
# result = Evolutionary.optimize(objective_function, constraints, de_algorithm, initial_population, Evolutionary.Options())
result = CustomEvolutionary1.optimize(objective_function, constraints, de_algorithm, initial_population, CustomEvolutionary1.Options())

# output = Evolutionary.optimize(f = objective_function, population = initial_population, constraints = constraints, method = de_algorithm)
    
print(result)





