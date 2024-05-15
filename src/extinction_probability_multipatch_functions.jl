module ExtinctionMultipatch

    using Intervals
    using Distributions
    using Polynomials
    using Evolutionary
    using CustomEvolutionary1
    using LinearAlgebra
    using IterTools
    using Parameters
    using Combinatorics
    using Random
    using JSON
    using JSON:json
    using Optimization, Optim
    using OptimizationOptimJL
    using ForwardDiff

    export survival_interval, get_survival_intervals, segment_unit_interval_by_survival_overlap, get_segment_survival_counts, get_segment_probability, get_segment_probabilities, get_probabilities_in_patch_given_H, get_probabilities_given_H, get_probabilities, probabilities2extinction_coefficients!, get_extinction_prob_from_coefficients, get_extinction_prob, make_objective_function, make_objective_function, get_idx2partition, initialize_population, survival_prob, maximize_survival_prob, near_mode_for_beta_mixture, diverse_push!, educated_guess, save_results, minimize_extinction_probability, make_objective_function, minimize_partition_extinction_probability


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
    function make_objective_function(fecundity::Int64, delta::Float64, alpha1::Float64, beta1::Float64, alpha2::Float64, beta2::Float64, p1::Float64, idx2partition::Dict)
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

    function survival_prob(center, delta, alpha1, beta1, alpha2, beta2, p1)
        dist = MixtureModel(Beta, [(alpha1, beta1), (alpha2, beta2)], [p1, 1-p1])

        return cdf(dist, center + delta) - cdf(dist, center - delta)
    end

    function maximize_survival_prob(delta, alpha, beta)
        
        objective = center -> -survival_prob(center[1], delta, alpha, beta)
        
        if alpha > 0 && beta > 0
            if alpha == beta
                best = 0.5
                desc = "best"
            else
                dist = Beta(alpha, beta)
                guess = [mode(dist)]
                guess = clamp.(guess, delta/2, 1-delta/2)
                # constraints = [delta/2,1 - delta/2]
                solution = optimize(objective, [delta/2], [1 - delta/2], guess, Fminbox(LBFGS()), Optim.Options(show_trace=false))
                best = solution.minimizer[1]

                desc = "best"
            end
        elseif 1 > alpha > beta
            best = 1 - delta/2, delta/2
            desc = "best", "good"
        elseif 1 > beta > alpha
            best = delta/2, 1 -  delta/2
            desc = "best", "good"
        elseif alpha > beta
            best = 1 - delta/2
            desc = "best"
        elseif beta < alpha
            best = delta/2
            desc = "best"
        elseif alpha == beta < 1
            @warn "two best: delta/2 and 1 - delta/2"
            best = delta/2, 1 - delta/2
            desc = "best", "best"
        elseif alpha == beta == 1
            @warn "std uniform dist"
            best = nothing
            desc = nothing
        end
    return best, desc
    end
    
    function near_mode_for_beta_mixture(alpha1, beta1, alpha2, beta2, p1)
        dist = MixtureModel(Beta, [(alpha1, beta1), (alpha2, beta2)], [p1, 1-p1])

        sample = [pdf(dist, el) for el in 0:0.01:1]

        return sample[argmax(sample)]
    end

    function maximize_survival_prob(delta, alpha1, beta1, alpha2, beta2, p1)
        
        objective = center -> -survival_prob(center[1], delta, alpha1, beta1, alpha2, beta2, p1)

        if alpha1 == beta1 == alpha2 == beta2 == 1 || (alpha1 == beta1 == 1 && p1 == 1) || (alpha1 == beta2 == 2 && alpha2 == beta1 == 0.5) # uniform dist
            # all centers equal
            @warn "Std Uniform distribution"
            best = nothing
        else 
            # dist = MixtureModel(Beta, [(alpha1, beta1), (alpha2, beta2)], [p1, 1-p1])
            guess = [near_mode_for_beta_mixture(alpha1, beta1, alpha2, beta2, p1)]
            guess = clamp.(guess, delta/2, 1-delta/2)
            # optf = OptimizationFunction(objective,Optimization.FiniteDiff())
            # problem = OptimizationProblem(optf, guess, lb = [delta/2], ub = [1 - delta/2])
            # solution = solve(problem, LBFGS())
            solution = optimize(objective, [delta/2], [1 - delta/2], guess, Fminbox(LBFGS()), Optim.Options(show_trace=false))

            # best = solution.u
            best = solution.minimizer[1]
        end
    return best
    end     

    ###### educated guessing
    function diverse_push!(vector, value)
        if value isa Float64 || value isa String
            push!(vector, value)
        else
            push!(vector, value[1])
            push!(vector, value[2])
        end
    end
    
    function get_educated_guess_candidates(delta, alpha1, beta1, alpha2, beta2, p1)
        candidates = Float64[]
        descriptions = String[]
        
        # maximize surv porb in single distributions
        candidates1, descriptions1 = maximize_survival_prob(delta, alpha1, beta1)
        candidates2, descriptions2 = maximize_survival_prob(delta, alpha2, beta2)
        
        diverse_push!(candidates, candidates1)
        diverse_push!(descriptions, descriptions1)
        diverse_push!(candidates, candidates2)
        diverse_push!(descriptions, descriptions2)

        # maximize surv prob in mixture dist
        candidate3 = maximize_survival_prob(delta, alpha1, beta1, alpha2, beta2, p1)
        description3 = "mean_max"
        push!(candidates, candidate3)
        push!(descriptions, description3)
        
        # might have eq surv prob to candidate3 in mixture dist or may be local max
        candidate4 = 1 - candidate3
        description4 = "mean_max_complement"
        push!(candidates, candidate4)
        push!(descriptions, description4)

        return candidates, descriptions
    end

    # function educated_guess(optimized, delta, alpha1, beta1, alpha2, beta2, p1)
    #     candidates, descriptions = get_educated_guess_candidates(delta, alpha1, beta1, alpha2, beta2, p1)

    #     distances = Float64[]
    #     guess_types = String[]
    #     for (el_cand, el_desc) in zip(candidates, descriptions)
    #         push!(distances, abs(optimized - el_cand))
    #         push!(guess_types, el_desc)
    #     end
    #     best_idx = argmin(distances)
    #     guess = candidates[best_idx]
    #     guess_type = guess_types[best_idx]

    #     return guess, guess_type
    # end

    function educated_guess(optimized, candidates, descriptions)
        # candidates, descriptions = get_educated_guess_candidates(delta, alpha1, beta1, alpha2, beta2, p1)

        distances = Float64[]
        guess_types = String[]
        for (el_cand, el_desc) in zip(candidates, descriptions)
            push!(distances, abs(optimized - el_cand))
            push!(guess_types, el_desc)
        end
        best_idx = argmin(distances)
        guess = candidates[best_idx]
        guess_type = guess_types[best_idx]

        return guess, guess_type
    end

    # function try_educated_guess(minimizer, extinction_probability, fecundity, delta, alpha1, beta1, alpha2, beta2, p1, objective_function)
    #     candidates, descriptions = get_educated_guess_candidates(delta, alpha1, beta1, alpha2, beta2, p1)
    #     replaced_by = fill("N/A", fecundity)
    #         current_pars = copy(minimizer)
    #         for (i, el) in enumerate(current_pars[1:fecundity])
    #             new_pars = copy(current_pars)
    #             new_pars[i], candidate_type = educated_guess(el, candidates, descriptions)
    #             new_extinction_probability = objective_function(new_pars)
    #             if new_extinction_probability < extinction_probability
    #                 current_pars = new_pars
    #                 extinction_probability = new_extinction_probability
    #                 replaced_by[i] = candidate_type
    #             end
    #         end
    #         centers = current_pars[1:fecundity]
    #     return centers, replaced_by, extinction_probability
    # end

    # function try_educated_guess(minimizer, extinction_probability, candidates, descriptions, objective_function, fecundity)
    #     replaced_by = fill("N/A", fecundity)
    #         current_pars = copy(minimizer)
    #         for (i, el) in enumerate(current_pars[1:fecundity])
    #             new_pars = copy(current_pars)
    #             new_pars[i], candidate_type = educated_guess(el, candidates, descriptions)
    #             new_extinction_probability = objective_function(new_pars)
    #             if new_extinction_probability < extinction_probability
    #                 current_pars = new_pars
    #                 extinction_probability = new_extinction_probability
    #                 replaced_by[i] = candidate_type
    #             end
    #         end
    #         centers = current_pars[1:fecundity]
    #     return centers, replaced_by, extinction_probability
    # end

    function try_educated_guess(centers, partition, extinction_probability, candidates, descriptions, delta, alpha1, beta1, alpha2, beta2, p1)
        fecundity = length(centers)
        replaced_by = fill("N/A", fecundity)
            current_centers = copy(centers)
            for (i, el) in enumerate(centers)
                new_centers = copy(centers)
                new_centers[i], candidate_type = educated_guess(el, candidates, descriptions)
                new_extinction_probability = get_extinction_prob(new_centers, partition, delta, alpha1, beta1, alpha2, beta2, p1, fecundity)
                if new_extinction_probability < extinction_probability
                    current_centers = new_centers
                    extinction_probability = new_extinction_probability
                    replaced_by[i] = candidate_type
                end
            end
            centers = current_centers
        return centers, replaced_by, extinction_probability
    end
     
    ###### helper functions
    function save_results(output, save_dir)

        # saving
        # Convert the dictionary to JSON format
        json_data = json(output)

        # Save the JSON string to a file
        save_dir = "output/" * save_dir
        mkpath(save_dir)

        open(save_dir * "/" * "output.json", "w") do file
            write(file, json_data)
        end

    end

    function get_constraints(delta, fecundity, idx2partition)
        lower_constraint = fill(delta/2, fecundity)
        upper_constraint = fill(1-delta/2, fecundity)

        lower_constraint = [lower_constraint; float(minimum(keys(idx2partition)))]
        upper_constraint = [upper_constraint; float(maximum(keys(idx2partition)))]

        constraints = CustomEvolutionary1.BoxConstraints(lower_constraint, upper_constraint)

        return constraints
    end

    function get_constraints(delta, fecundity)
        lower_constraint = fill(delta/2, fecundity)
        upper_constraint = fill(1-delta/2, fecundity)

        constraints = CustomEvolutionary1.BoxConstraints(lower_constraint, upper_constraint)

        return constraints
    end

    ##### running everything
    function minimize_extinction_probability(fecundity::Int64, delta::Float64, alpha1::Float64, beta1::Float64, alpha2::Float64, beta2::Float64, p1::Float64, save_dir::String, population_size::Int64, partition_mutation_rate::Float64, use_educated_guess::Bool)
        
        idx2partition = get_idx2partition(fecundity)

        objective_function = make_objective_function(fecundity, delta, alpha1, beta1, alpha2, beta2, p1, idx2partition)

        initial_population = initialize_population(population_size, fecundity, idx2partition)
        custom_differentiation = create_custom_differentiation(fecundity, partition_mutation_rate, idx2partition)

        lower_constraint = fill(delta/2, fecundity)
        upper_constraint = fill(1-delta/2, fecundity)
        lower_constraint = [lower_constraint; float(minimum(keys(idx2partition)))]
        upper_constraint = [upper_constraint; float(maximum(keys(idx2partition)))]
        constraints = CustomEvolutionary1.BoxConstraints(lower_constraint, upper_constraint)

        de_algorithm = CustomEvolutionary1.DE(populationSize = population_size, differentiation = custom_differentiation)

        results = CustomEvolutionary1.optimize(objective_function, constraints, de_algorithm, initial_population)

        # exctract results
        centers = results.minimizer[1:fecundity]
        partition = idx2partition[results.minimizer[end]]
        extinction_probability = results.minimum

        # educated guess true optimum
        if use_educated_guess
            candidates, descriptions = get_educated_guess_candidates(delta, alpha1, beta1, alpha2, beta2, p1)
            centers, replaced_by, extinction_probability = try_educated_guess(centers, partition, extinction_probability, candidates, descriptions, delta, alpha1, beta1, alpha2, beta2, p1)
        end

        # gets mean fitness maximizer 
        # mean_maximizer = maximize_survival_prob(delta, alpha1, beta1, alpha2, beta2, p1)

        # output
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

        if use_educated_guess
            output["replaced_by"] = replaced_by
            output["special_centers"] = hcat(candidates, descriptions)
        end
        
        save_results(output, save_dir)

        return output
    end

    # brute force functions
    function make_objective_function(partition::Vector{Vector{Int64}}, fecundity::Int64, delta::Float64, alpha1::Float64, beta1::Float64, alpha2::Float64, beta2::Float64, p1::Float64)
        function objective_function(vars)
            centers = vars[1:fecundity]
            extinction_prob = get_extinction_prob(centers, partition, delta, alpha1, beta1, alpha2, beta2, p1, fecundity)
            return extinction_prob
        end

        return objective_function
    end


    function minimize_partition_extinction_probability(partition::Vector{Vector{Int64}}, fecundity::Int64, delta::Float64, alpha1::Float64, beta1::Float64, alpha2::Float64, beta2::Float64, p1::Float64, population_size::Int64, use_educated_guess::Bool)

        # idx2partition = get_idx2partition(fecundity)

        objective_function = make_objective_function(partition, fecundity, delta, alpha1, beta1, alpha2, beta2, p1)

        # initial_population = initialize_population(population_size, fecundity, idx2partition)
        # custom_differentiation = create_custom_differentiation(fecundity, partition_mutation_rate, idx2partition)

        lower_constraint = fill(delta/2, fecundity)
        upper_constraint = fill(1-delta/2, fecundity)

        # lower_constraint = [lower_constraint; float(minimum(keys(idx2partition)))]
        # upper_constraint = [upper_constraint; float(maximum(keys(idx2partition)))]

        constraints = Evolutionary.BoxConstraints(lower_constraint, upper_constraint)

        de_algorithm = Evolutionary.DE(populationSize = population_size)

        results = Evolutionary.optimize(objective_function, constraints, de_algorithm)

        # exctract results
        centers = results.minimizer[1:fecundity]
        extinction_probability = results.minimum

        # educated guess true optimum
        if use_educated_guess
            replaced_by = fill("N/A", fecundity)
            current_pars = copy(results.minimizer)
            for (i, el) in enumerate(current_pars[1:fecundity])
                new_pars = copy(current_pars)
                new_pars[i], candidate_type = educated_guess(el, delta, alpha1, beta1, alpha2, beta2, p1)
                new_extinction_probability = objective_function(new_pars)
                if new_extinction_probability < extinction_probability
                    current_pars = new_pars
                    extinction_probability = new_extinction_probability
                    replaced_by[i] = candidate_type
                end
            end
            centers = current_pars[1:fecundity]
        end

        # gets mean fitness maximizer 
        mean_maximizer = maximize_survival_prob(delta, alpha1, beta1, alpha2, beta2, p1)

        # output
        output = Dict("fecundity" => fecundity,
                    "delta" => delta,
                    "alpha1" => alpha1,
                    "beta1" => beta1,
                    "alpha2" => alpha2,
                    "beta2" => beta2,
                    "p1" => p1,
                    "centers" => centers,
                    "partition" => partition,
                    "extinction_probability" => extinction_probability,
                    "mean_maximizer" => mean_maximizer
                    )

        if use_educated_guess
            output["replaced_by"] = replaced_by
        end
        
        return output
    end

    function minimize_extinction_probability(fecundity::Int64, delta::Float64, alpha1::Float64, beta1::Float64, alpha2::Float64, beta2::Float64, p1::Float64, save_dir::String, population_size::Int64, use_educated_guess::Bool)

        extinction_probability = 10
        output = nothing

        all_partitions = collect(partitions(1:fecundity))
        
        for el in all_partitions
            current_output = minimize_partition_extinction_probability(el, fecundity, delta, alpha1, beta1, alpha2, beta2, p1, population_size, use_educated_guess)
            
            if current_output["extinction_probability"] < extinction_probability
                extinction_probability = current_output["extinction_probability"]
                output = current_output
            end
        end

        save_results(output, save_dir)
        
        return output
    end
end