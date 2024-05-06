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
using JSON



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
function make_objective_function(fecundity, delta, alpha1, beta1, alpha2, beta2, p1, idx2partition)
    function objective_function(vars)
        centers = vars[1:fecundity]
        partition_idx = vars[end]
        partition = idx2partition[partition_idx]
        extinction_prob = get_extinction_prob(centers, partition, delta, alpha1, beta1, alpha2, beta2, p1, fecundity)
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
    

function main()
    s = ArgParseSettings()  # Create a settings object

    @add_arg_table s begin
        "--fecundity"
            arg_type = Int64
            required = true
        "--delta"
            arg_type = Float64
            required = true
        "--alpha1"
            arg_type = Float64
            required = true
        "--beta1"
            arg_type = Float64
            required = true
        "--alpha2"
            arg_type = Float64
            required = true
        "--beta2"
            arg_type = Float64
            required = true
        "--p1"
            arg_type = Float64
            required = true
        "--save_dir"
            arg_type = String
            required = true
        "--population_size"
            arg_type = Int64
            required = false
            default = 100
        "--partition_mutation_rate"
            arg_type = Float64
            required = false
            default = 0.2
    end

    args = parse_args(s)  # This function parses the command line arguments based on the specified settings


    idx2partition = get_idx2partition(args["fecundity"])

    objective_function = make_objective_function(args["fecundity"], args["delta"], args["alpha1"], args["beta1"], args["alpha2"], args["beta2"], args["p1"], idx2partition)

    initial_population = initialize_population(args["population_size"], args["fecundity"], idx2partition)
    custom_differentiation = create_custom_differentiation(args["fecundity"], args["partition_mutation_rate"], idx2partition)

    lower_constraint = fill(0.0, args["fecundity"])
    upper_constraint = fill(1.0, args["fecundity"])

    lower_constraint = [lower_constraint; float(minimum(keys(idx2partition)))]
    upper_constraint = [upper_constraint; float(maximum(keys(idx2partition)))]

    constraints = BoxConstraints(lower_constraint, upper_constraint)

    de_algorithm = DE(populationSize = args["population_size"], differentiation = custom_differentiation)

    results = CustomEvolutionary1.optimize(objective_function, constraints, de_algorithm, initial_population, CustomEvolutionary1.Options())

    # output
    centers = results.minimizer[1:args["fecundity"]]
    partition = idx2partition[results.minimizer[end]]
    extinction_probability = results.minimum

    output = Dict("fecundity" => args["fecundity"],
                  "delta" => args["delta"],
                  "alpha1" => args["alpha1"],
                  "beta1" => args["beta1"],
                  "alpha2" => args["alpha2"],
                  "beta2" => args["beta2"],
                  "p1" => args["p1"],
                  "centers" => centers,
                  "partition" => partition,
                  "extinction_probability" => extinction_probability
                  )
    
    # saving
    # Convert the dictionary to JSON format
    json_data = JSON.json(output)

    # Save the JSON string to a file
    save_dir = "output/" * args["save_dir"]
    mkpath(save_dir)

    open(save_dir * "/" * "output.json", "w") do file
        write(file, json_data)
    end
    
    return output
end

# This conditional ensures that the script runs main only if it is not being included as a module
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end








