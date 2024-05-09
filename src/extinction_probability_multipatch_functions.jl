module ExtinctionMultipatch

    using Intervals
    using Distributions
    using Polynomials
    using CustomEvolutionary1
    using CustomEvolutionary1: optimize
    using LinearAlgebra 
    using IterTools
    using Parameters
    using Combinatorics
    using Random
    using JSON
    using JSON:json
    using Optimization
    using OptimizationOptimJL
    using ForwardDiff


    export initialize_population, get_idx2partition, make_objective_function, create_custom_differentiation, BoxConstraints, DE, optimize, json, survival_prob, maximize_survival_prob, educated_guess, minimize_extinction_probability

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

    function survival_prob(center, delta, alpha, beta)
        dist = Beta(alpha, beta)
        return cdf(dist, center + delta) - cdf(dist, center - delta)
    end

    function maximize_survival_prob(delta, alpha, beta)
        
        objective = center -> -survival_prob(center, delta, alpha, beta)
        
        if alpha > 0 && beta > 0
            if alpha == beta
                best = 0.5
            else
                dist = Beta(alpha, beta)
                guess = mode(dist)
                guess = clamp(guess, delta/2, 1-delta/2)
                # constraints = [delta/2,1 - delta/2]
                optf = OptimizationFunction(objective, Optimization.AutoForwardDiff())
                problem = OptimizationProblem(optf, guess, lb = delta/2, ub = 1 - delta/2)
                solution = solve(problem, LBFGS())
                best = solution.u
            end
        elseif alpha > beta
            best = 1 - delta/2
        elseif beta < alpha
            best = 1 + delta/2
        elseif alpha == beta < 1
            @warn "two best: delta/2 and 1 - delta/2"
            best = nothing
        elseif alpha == beta == 1
            @warn "std uniform dist"
            best = nothing
        end
    return best
    end        
        
    function educated_guess(optimized, delta, alpha1, beta1, alpha2, beta2, p1)
        best1 = maximize_survival_prob(delta, alpha1, beta1)
        best2 = maximize_survival_prob(delta, alpha2, beta2)

        if isnothing(best)
            return nothing
        else
            distance1 = abs(optimized - best1)
            distance2 = abs(optimized - best2)

            best_idx = argmin([distance1, distance2])
            guess = [best1, best2][best_idx]

            return guess
        end
    end

    function minimize_extinction_probability(fecundity, delta, alpha1, beta1, alpha2, beta2, p1, save_dir, population_size, partition_mutation_rate)

        idx2partition = get_idx2partition(fecundity)

        objective_function = make_objective_function(fecundity, delta, alpha1, beta1, alpha2, beta2, p1, idx2partition)

        initial_population = initialize_population(population_size, fecundity, idx2partition)
        custom_differentiation = create_custom_differentiation(fecundity, partition_mutation_rate, idx2partition)

        lower_constraint = fill(0.0, fecundity)
        upper_constraint = fill(1.0, fecundity)

        lower_constraint = [lower_constraint; float(minimum(keys(idx2partition)))]
        upper_constraint = [upper_constraint; float(maximum(keys(idx2partition)))]

        constraints = BoxConstraints(lower_constraint, upper_constraint)

        de_algorithm = DE(populationSize = population_size, differentiation = custom_differentiation)

        results = optimize(objective_function, constraints, de_algorithm, initial_population)

        # output
        centers = results.minimizer[1:fecundity]
        partition = idx2partition[results.minimizer[end]]
        extinction_probability = results.minimum

        output = Dict("fecundity" => fecundity,
                    "delta" => delta,
                    "alpha1" => alpha1,
                    "beta1" => beta1,
                    "alpha2" => alpha2,
                    "beta2" => beta2,
                    "p1" => p1,
                    "centers" => centers,
                    "partition" => partition,
                    "extinction_probability" => extinction_probability
                    )
        
        # saving
        # Convert the dictionary to JSON format
        json_data = json(output)

        # Save the JSON string to a file
        save_dir = "output/" * save_dir
        mkpath(save_dir)

        open(save_dir * "/" * "output.json", "w") do file
            write(file, json_data)
        end

        return output
    end
    
end