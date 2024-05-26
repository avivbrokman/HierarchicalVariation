# For mean-maximizer strategy, if arrived at through optimization, check whether any of the single distribution best strategies are better 
# Figure out why the replacement thingy isn't working.
# removed duplicates for "good" and "best" and "mean-maximizer"

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
    using StatsBase

    export survival_interval, get_survival_intervals, segment_unit_interval_by_survival_overlap, get_segment_survival_counts, get_segment_probability, get_segment_probabilities, get_probabilities_in_patch_given_H, get_probabilities_given_H, get_probabilities, probabilities2extinction_coefficients!, get_extinction_prob_from_coefficients, get_extinction_prob, make_objective_function, make_objective_function, get_idx2partition, initialize_population, survival_prob, maximize_survival_prob, near_mode_for_beta_mixture, candidate_push!, educated_guess, save_results, minimize_extinction_probability, make_objective_function, minimize_partition_extinction_probability, get_desc2priority, get_educated_guess_candidates_from_single_distributions, try_educated_guess, get_educated_guess_candidates, is_close, close_replace


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
        
        if alpha > 1 && beta > 1
            if alpha == beta
                best = 0.5
                desc = "best"
            else
                dist = Beta(alpha, beta)
                guess = [mode(dist)]
                guess = clamp.(guess, delta, 1-delta)
                solution = optimize(objective, [delta], [1 - delta], guess, Fminbox(LBFGS()), Optim.Options(show_trace=false))
                best = solution.minimizer[1]

                desc = "best"
            end
        elseif 1 > alpha > beta
            best = [1 - delta, delta]
            desc = ["best", "good"]
        elseif 1 > beta > alpha
            best = [delta, 1 -  delta]
            desc = ["best", "good"]
        elseif alpha > beta
            best = 1 - delta
            desc = "best"
        elseif beta < alpha
            best = delta
            desc = "best"
        elseif alpha == beta < 1
            @warn "two best: delta and 1 - delta"
            best = [delta, 1 - delta]
            desc = ["best", "best"]
        elseif alpha == beta == 1
            @warn "std uniform dist"
            best = nothing
            desc = nothing
        elseif 1 == alpha > beta
            best = 1 - delta
            desc = "best"
        elseif 1 == beta > alpha
            best = delta
            desc = "best"
        end
    return best, desc
    end

    function maximize_survival_prob(delta, alpha, beta)
        
        objective = center -> -survival_prob(center[1], delta, alpha, beta)
        

        if alpha == beta == 1
            @warn "std uniform dist"
            best = nothing
            desc = nothing
        elseif alpha == beta < 1
            @warn "two best: delta and 1 - delta"
            best = [delta, 1 - delta]
            desc = ["best", "best"]
        elseif alpha == beta > 1
            alpha == beta
            best = 0.5
            desc = "best"
        elseif alpha > 1 && beta > 1 && alpha != beta
            dist = Beta(alpha, beta)
            guess = [mode(dist)]
            guess = clamp.(guess, delta, 1-delta)
            solution = optimize(objective, [delta], [1 - delta], guess, Fminbox(LBFGS()), Optim.Options(show_trace=false))
            best = solution.minimizer[1]
            desc = "best"
        elseif 1 <= alpha > beta <= 1
            best = 1 - delta
            desc = "best"
        elseif 1 <= beta > alpha <= 1
            best = delta
            desc = "best"
        elseif 1 > alpha > beta
            best = [1 - delta, delta]
            desc = ["best", "good"]
        elseif 1 > beta > alpha
            best = [delta, 1 -  delta]
            desc = ["best", "good"]        
        end
    return best, desc
    end
    
    function near_mode_for_beta_mixture(alpha1, beta1, alpha2, beta2, p1)
        dist = MixtureModel(Beta, [(alpha1, beta1), (alpha2, beta2)], [p1, 1-p1])

        sample = [pdf(dist, el) for el in 0:0.01:1]

        return sample[argmax(sample)]
    end

    function is_close(val1, val2, tol = 1e-8)
        return abs(val1 - val2) < tol
    end

    function close_replace(candidate, desc, optimized, tol = 1e-8)
        if is_close(val1, val2, tol)
            return candidate, desc
        else
            return optimized, nothing
        end
    end


    function try_educated_guess(center::Float64, center_survival_prob::Float64, candidates::Vector{Float64}, descriptions::Vector{String}, delta::Float64, alpha1::Float64, beta1::Float64, alpha2::Float64, beta2::Float64, p1::Float64)
        desc = "nothing"
        for (el_cand, el_desc) in zip(candidates, descriptions)
            if is_close(el_cand, center)
                new_survival_prob = survival_prob(el_cand, delta, alpha1, beta1, alpha2, beta2, p1)
                center = el_cand
                center_survival_prob = new_survival_prob
                desc = el_desc
            else
                new_survival_prob = survival_prob(el_cand, delta, alpha1, beta1, alpha2, beta2, p1)
                if new_survival_prob >= center_survival_prob
                    center = el_cand
                    center_survival_prob = new_survival_prob
                    desc = el_desc
                end
            end
        end
        return center, desc
    end
    
    function maximize_survival_prob(delta, alpha1, beta1, alpha2, beta2, p1)
        
        objective = center -> -survival_prob(center[1], delta, alpha1, beta1, alpha2, beta2, p1)

        if alpha1 == beta1 == alpha2 == beta2 == 1 || (alpha1 == beta1 == 1 && p1 == 1) || (alpha1 == beta2 == 2 && alpha2 == beta1 == 1) # uniform dist
            # all centers equal
            @warn "Std Uniform distribution"
            return nothing, nothing
        else 
            # dist = MixtureModel(Beta, [(alpha1, beta1), (alpha2, beta2)], [p1, 1-p1])
            guess = [near_mode_for_beta_mixture(alpha1, beta1, alpha2, beta2, p1)]
            guess = clamp.(guess, delta, 1-delta)
            solution = optimize(objective, [delta], [1 - delta], guess, Fminbox(LBFGS()), Optim.Options(show_trace=false))

            # best = solution.u
            best = solution.minimizer[1]
            value = -solution.minimum

            # educated guessing
            candidates, descriptions = get_educated_guess_candidates_from_single_distributions(delta, alpha1, beta1, alpha2, beta2)

            best, desc = try_educated_guess(best, value, candidates, descriptions, delta, alpha1, beta1, alpha2, beta2, p1)

            return best, desc
        end
     end     

    ###### educated guessing

    # function diverse_push!(vector, value)
    #     if value isa Float64 || value isa String
    #         if value ∉ vector
    #             push!(vector, value)
    #         end
    #     else
    #         if value[1] ∉ vector
    #             push!(vector, value[1])
    #         end
    #         if if value[2] ∉ vector
    #             push!(vector, value[2])
    #         end
    #     end
    # end


    # function diverse_push!(vector::Vector, value::Float64)
    #     if value ∉ vector
    #         push!(vector, value)
    #     end
    # end

    # function diverse_push!(vector::Vector, value::Vector)
    #     diverse_push!(vector, value[1])
    #     diverse_push!(vector, value[2])
    # end

    function get_desc2priority()
        priority_order = ["best", "mean_max", "good", "mean_max_complement", "N/A"]
        desc2priority = Dict(zip(priority_order, 1:length(priority_order)))

        return desc2priority
    end
    #Vector{Union{Float64, Nothing}}()
    function candidate_push!(candidate_vector::Vector, desc_vector::Vector, candidate::Union{Float64, Nothing}, desc::Union{String, Nothing}, desc2priority::Dict)
        if isnothing(candidate)
            nothing 
        elseif candidate ∉ candidate_vector 
            push!(candidate_vector, candidate)
            push!(desc_vector, desc)
        else
            # index = findfirst(isequal(candidate), candidate_vector)
            index = findfirst(x -> is_close(x, candidate), candidate_vector)
            old_desc = desc_vector[index]
            old_priority = desc2priority[old_desc]
            new_priority = desc2priority[desc]
            if new_priority < old_priority
                deleteat!(candidate_vector, index)
                deleteat!(desc_vector, index)
                push!(candidate_vector, candidate)
                push!(desc_vector, desc)
            end
        end
    end

    function candidate_push!(candidate_vector::Vector{Float64}, desc_vector::Vector{String}, candidates::Vector, descs::Vector, desc2priority::Dict)
        for (el_cand, el_desc) in zip(candidates, descs)
            candidate_push!(candidate_vector, desc_vector, el_cand, el_desc, desc2priority)
        end
    end

    function get_educated_guess_candidates_from_single_distributions(delta, alpha1, beta1, alpha2, beta2)
        
        desc2priority = get_desc2priority()
        
        # candidates = Vector{Union{Float64, Nothing}}()
        # descriptions = Vector{Union{String, Nothing}}()

        candidates = Float64[]
        descriptions = String[]
        
        # maximize surv porb in single distributions
        candidates1, descriptions1 = maximize_survival_prob(delta, alpha1, beta1)
        candidates2, descriptions2 = maximize_survival_prob(delta, alpha2, beta2)

        candidate_push!(candidates, descriptions, candidates1, descriptions1, desc2priority)
        candidate_push!(candidates, descriptions, candidates2, descriptions2, desc2priority)

        return candidates, descriptions
    end

    function get_educated_guess_candidates(delta, alpha1, beta1, alpha2, beta2, p1)
        desc2priority = get_desc2priority()

        candidates, descriptions = get_educated_guess_candidates_from_single_distributions(delta, alpha1, beta1, alpha2, beta2)

        # maximize surv prob in mixture dist
        candidate3, description3 = maximize_survival_prob(delta, alpha1, beta1, alpha2, beta2, p1)
        # if isnothing(description3)
        #     description3 = "mean_max"
        # end
        # push!(candidates, candidate3)
        # push!(descriptions, description3)
        candidate_push!(candidates, descriptions, candidate3, description3, desc2priority)

        
        # might have eq surv prob to candidate3 in mixture dist or may be local max
        candidate4 = 1 - candidate3
        description4 = "mean_max_complement"
        # push!(candidates, candidate4)
        # push!(descriptions, description4)
        candidate_push!(candidates, descriptions, candidate4, description4, desc2priority)


        return candidates, descriptions
    end

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

    function try_educated_guess(centers::Vector{Float64}, partition::Vector{Vector{Int64}}, extinction_probability::Float64, candidates::Vector{Float64}, descriptions::Vector{String}, delta::Float64, alpha1::Float64, beta1::Float64, alpha2::Float64, beta2::Float64, p1::Float64)
        fecundity = length(centers)
        replaced_by = fill("N/A", fecundity)
        is_altered = true
        while is_altered
            is_altered = false
            for i in 1:fecundity
                new_centers = copy(centers)
                for (el_cand, el_desc) in zip(candidates, descriptions)
                    if is_close(el_cand, centers[i])
                        centers[i] = el_cand
                        replaced_by[i] = el_desc
                        is_altered = True
                    else
                        new_centers[i] = el_cand
                        new_extinction_probability = get_extinction_prob(new_centers, partition, delta, alpha1, beta1, alpha2, beta2, p1, fecundity)
                        if new_extinction_probability <= extinction_probability
                            centers = new_centers
                            extinction_probability = new_extinction_probability
                            replaced_by[i] = el_desc
                            is_altered = True
                        end
                    end
                end
            end
        end
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

    ##### running everything
    function minimize_extinction_probability(fecundity::Int64, delta::Float64, alpha1::Float64, beta1::Float64, alpha2::Float64, beta2::Float64, p1::Float64, save_dir::String, population_size::Int64, partition_mutation_rate::Float64, use_educated_guess::Bool)
        
        idx2partition = get_idx2partition(fecundity)

        objective_function = make_objective_function(fecundity, delta, alpha1, beta1, alpha2, beta2, p1, idx2partition)

        initial_population = initialize_population(population_size, fecundity, idx2partition)
        custom_differentiation = create_custom_differentiation(fecundity, partition_mutation_rate, idx2partition)

        lower_constraint = fill(delta, fecundity)
        upper_constraint = fill(1-delta, fecundity)
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

        lower_constraint = fill(delta, fecundity)
        upper_constraint = fill(1-delta, fecundity)

        constraints = Evolutionary.BoxConstraints(lower_constraint, upper_constraint)

        de_algorithm = Evolutionary.DE(populationSize = population_size)

        results = Evolutionary.optimize(objective_function, constraints, de_algorithm)

        # exctract results
        centers = results.minimizer[1:fecundity]
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

    # def f_eps(s, epsilon, p):
    #     conditional_gen_fun_coef = p[epsilon]
    #     val = sum(el * s**i for i, el in enumerate(conditional_gen_fun_coef))
    #     return val

    # survival = s, environment = ϵ, generating_fun_coef_by_environment = p
    function f_eps(survival, environment, generating_fun_coef_by_environment)
        conditional_generating_fun_coef = generating_fun_coef_by_environment[environment]
        val = sum(el * survival^i for (i, el) in enumerate(conditional_generating_fun_coef))
        return val
    end
    
    function generate_environment_sequence(p1, num_generations)
        return StatsBase.sample([1,2], Weights([p1, 1 - p1]), num_generations, replace = true)
    end

    # function generate_environment_sequences(p1, num_generations, num_runs)
    #     return [generate_environment_sequence(p1, num_generations) for el in 1:num_runs]
    # end

    # generating_fun_coef_by_environment = p, 
    # function approximate_extinction_probability(generating_fun_coef_by_environment, q0, environments::Vector{int64})
        
    #     q_seq = [q0]
    #     for el in environments
    #         q_new = f_eps(q_seq[end], el, generating_fun_coef_by_environment)
    #         q_seq.append(q_new)
    #     end
            
    #     return q_seq
    # end

    function approximate_extinction_probability(generating_fun_coef_by_environment, q0, environments::Vector{int64})
        
        q = q0
        for el in environments
            q = f_eps(q, el, generating_fun_coef_by_environment)
        end
            
        return q
    end

    function approximate_extinction_probability(generating_fun_coef_by_environment, q0, p1, num_generations)
        generate_environment_sequence(p1, num_generations)
        q = q0
        for el in environments
            q = f_eps(q, el, generating_fun_coef_by_environment)
        end
            
        return q
    end

    function approximate_extinction_probability(generating_fun_coef_by_environment, q0, p1, num_generations, num_runs)
        
        extinction_prob_runs = [approximate_extinction_probability(generating_fun_coef_by_environment, q0, p1, num_generations) for el in 1:num_runs]

        extinction_prob = mean(extinction_prob_runs)

        return extinction_prob
    end

    # function single_run_wrapper(p, q0, p0, num_iters = 500)
    #     environments = generate_environments(p0, num_iters)
    #     return single_run(p, q0, environments)
    # end

    # function multi_run(p, p0, num_iters = 500)
    #     environments = generate_environments(p0, num_iters)
    #     q0s = arange(0.01, 1, 0.01)
    #     runs = [single_run(p, el, environments) for el in q0s]
    #     # return np.array(runs)
    # end
end