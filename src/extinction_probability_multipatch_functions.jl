
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
    using Profile

    export survival_interval, get_survival_intervals, segment_unit_interval_by_survival_overlap, get_segment_survival_counts, get_segment_probability, get_segment_probabilities, get_probabilities_in_patch_given_H, get_probabilities_given_H, get_probabilities, probabilities2extinction_coefficients!, get_extinction_probability_from_coefficients, get_extinction_probability
    
    export f_eps, generate_environment_sequence, approximate_extinction_probability

    export make_objective_function, get_idx2partition, initialize_population
    
    export survival_probability, maximize_survival_probability, near_mode_for_beta_mixture, special_center_push!, educated_guess, get_special_centers_from_single_distributions, get_educated_guesses, try_educated_guesses, is_close, modify_description!, get_special_centers, find_matches
    
    export save_results, minimize_extinction_probability, minimize_partition_extinction_probability, get_partitions, partition2signature

    # export Beta, Interval, cdf, Polynomial, roots

    ############### Basic calculations
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

    function get_extinction_probability_from_coefficients(coefficients::Vector{Float64})
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
        extinction_probability = minimum(positive_roots)

        return extinction_probability
    end

    # survival = s, environment = ϵ, generating_fun_coef_by_environment = p
    function f_eps(survival, environment, generating_fun_coef_by_environment)
        # print("f_eps")
        conditional_generating_fun_coef = generating_fun_coef_by_environment[environment]
        val = sum(el * survival^(i-1) for (i, el) in enumerate(conditional_generating_fun_coef))
        return val
    end
    
    function generate_environment_sequence(p1, num_generations)
        # print("generate environment sequence")
        return StatsBase.sample([1,2], Weights([p1, 1 - p1]), num_generations, replace = true)
    end

    function approximate_extinction_probability(probabilities_by_environment, q0, environments::Vector{Int64})
        # print("approximate_extinction_probability 11111")

        q = q0
        for el in environments
            q = f_eps(q, el, probabilities_by_environment)
        end
            
        return q
    end

    function approximate_extinction_probability(probabilities_by_environment, q0, p1, num_generations)
        # print("approximate_extinction_probability 22222")
        environments = generate_environment_sequence(p1, num_generations)
        q = q0
        for el in environments
            q = f_eps(q, el, probabilities_by_environment)
        end
            
        return q
    end

    function approximate_extinction_probability(probabilities_by_environment, q0, p1, num_generations, num_runs)
        # print("approximate 33333")

        extinction_probability_runs = [approximate_extinction_probability(probabilities_by_environment, q0, p1, num_generations) for el in 1:num_runs]

        extinction_probability = mean(extinction_probability_runs)

        return extinction_probability
    end

    ## workhorse function
    function get_extinction_probability(centers, partition, delta, alpha1, beta1, alpha2, beta2, p1, fecundity) 
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

        extinction_probability = get_extinction_probability_from_coefficients(probabilities)

        return extinction_probability
    end

    function get_extinction_probability(centers, partition, delta, alpha1, beta1, alpha2, beta2, p1, fecundity, q0, num_generations, num_runs) 
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

        extinction_probability = approximate_extinction_probability(probabilities_by_environment, q0, p1, num_generations, num_runs)

        return extinction_probability
    end

    ################### Necessary for running DE
    ## using a closure as a workaround because DE's objective function must take a single vector as input
    function make_objective_function(fecundity::Int64, delta::Float64, alpha1::Number, beta1::Number, alpha2::Number, beta2::Number, p1::Number, idx2partition::Dict, args...)
        function objective_function(vars)
            centers = vars[1:fecundity]
            partition_idx = vars[end]
            partition = idx2partition[partition_idx]
            extinction_probability = get_extinction_probability(centers, partition, delta, alpha1, beta1, alpha2, beta2, p1, fecundity, args...)
            return extinction_probability
        end

        return objective_function
    end

    function get_idx2partition(fecundity)
        all_partitions = collect(partitions(1:fecundity))
        idx2partition = Dict(float(i) => el for (i, el) in enumerate(all_partitions))
        return idx2partition
    end

    # brute force functions
    function make_objective_function(partition::Vector{Vector{Int64}}, fecundity::Int64, delta::Float64, alpha1::Number, beta1::Number, alpha2::Number, beta2::Number, p1::Number, args...)
        function objective_function(vars)
            centers = vars[1:fecundity]
            extinction_probability = get_extinction_probability(centers, partition, delta, alpha1, beta1, alpha2, beta2, p1, fecundity, args...)
            return extinction_probability
        end

        return objective_function
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

    #################### Educated guessing

    function survival_probability(center, delta, alpha, beta)
        dist = Beta(alpha, beta)
        return cdf(dist, center + delta) - cdf(dist, center - delta)
    end

    function survival_probability(center, delta, alpha1, beta1, alpha2, beta2, p1)
        dist = MixtureModel(Beta, [(alpha1, beta1), (alpha2, beta2)], [p1, 1-p1])

        return cdf(dist, center + delta) - cdf(dist, center - delta)
    end

    function maximize_survival_probability(delta, alpha, beta)
        
        objective = center -> -survival_probability(center[1], delta, alpha, beta)
        
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

    function modify_description!(description::String, modifier)
        description *= modifier
    end

    function modify_description!(description::Nothing, modifier)
    end

    function modify_description!(description::Vector, modifier)
        for el in description
            modify_description!(el, modifier)
        end
    end

    function maximize_survival_probability(delta, alpha, beta, description_modifier = "")
        
        objective = center -> -survival_probability(center[1], delta, alpha, beta)
        

        if alpha == beta == 1
            @warn "std uniform dist"
            best = nothing
            description = nothing
        elseif alpha == beta < 1
            @warn "two best: delta and 1 - delta"
            best = [delta, 1 - delta]
            description = ["best", "best"]
        elseif alpha == beta > 1
            alpha == beta
            best = 0.5
            description = "best"
        elseif alpha > 1 && beta > 1 && alpha != beta
            dist = Beta(alpha, beta)
            guess = [mode(dist)]
            guess = clamp.(guess, delta, 1-delta)
            solution = optimize(objective, [delta], [1 - delta], guess, Fminbox(LBFGS()), Optim.Options(show_trace=false))
            best = solution.minimizer[1]
            description = "best"
        elseif 1 <= alpha > beta <= 1
            best = 1 - delta
            description = "best"
        elseif 1 <= beta > alpha <= 1
            best = delta
            description = "best"
        elseif 1 > alpha > beta
            best = [1 - delta, delta]
            description = ["best", "good"]
        elseif 1 > beta > alpha
            best = [delta, 1 -  delta]
            description = ["best", "good"]        
        end

        modify_description!(description, description_modifier)
        
        return best, description
    end
    
    function near_mode_for_beta_mixture(alpha1, beta1, alpha2, beta2, p1)
        dist = MixtureModel(Beta, [(alpha1, beta1), (alpha2, beta2)], [p1, 1-p1])

        sample_ = [pdf(dist, el) for el in 0:0.001:1]

        return sample_[argmax(sample_)]
    end

    function special_center_push!(special_centers::Vector, descriptions::Vector, new_center::Float64, new_description::String)
        if (new_center, new_description) ∉ zip(special_centers, descriptions)
            push!(special_centers, new_center)
            push!(descriptions, new_description)
        end
    end

    function special_center_push!(special_centers::Vector, descriptions::Vector, new_center::Nothing, new_description::Nothing)
    end
    function special_center_push!(special_centers::Vector, descriptions::Vector, new_center::Nothing, new_description::String)
    end

    function special_center_push!(special_centers::Vector, descriptions::Vector, new_centers::Vector, new_descriptions::Vector)
        for (el_cent, el_desc) in zip(new_centers, new_descriptions)
            special_center_push!(special_centers, descriptions, el_cent, el_desc)
        end
    end

    function get_special_centers_from_single_distributions(delta, alpha1, beta1, alpha2, beta2, p1)
        
        special_centers = Float64[]
        descriptions = String[]
        
        # maximize surv porb in single distributions
        if (alpha1, beta1) == (alpha2, beta2) || p1 == 1
            special_center, description = maximize_survival_probability(delta, alpha1, beta1)
            special_center_push!(special_centers, descriptions, special_center, description)
            return special_centers, descriptions
        elseif p1 == 0.5
            special_centers1, descriptions1 = maximize_survival_probability(delta, alpha1, beta1)
            special_centers2, descriptions2 = maximize_survival_probability(delta, alpha2, beta2)
        else
            special_centers1, descriptions1 = maximize_survival_probability(delta, alpha1, beta1, "major_")
            special_centers2, descriptions2 = maximize_survival_probability(delta, alpha2, beta2, "minor_")
        end

        special_center_push!(special_centers, descriptions, special_centers1, descriptions1)
        special_center_push!(special_centers, descriptions, special_centers2, descriptions2)

        return special_centers, descriptions
    end 

    function maximize_survival_probability(delta, alpha1, beta1, alpha2, beta2, p1)
        
        objective = center -> -survival_probability(center[1], delta, alpha1, beta1, alpha2, beta2, p1)

        if alpha1 == beta1 == alpha2 == beta2 == 1 || (alpha1 == beta1 == 1 && p1 == 1) || (alpha1 == beta2 == 2 && alpha2 == beta1 == 1) # uniform dist
            # all centers equal
            @warn "Std Uniform distribution"
            return nothing
        else 
            guess = [near_mode_for_beta_mixture(alpha1, beta1, alpha2, beta2, p1)]
            guess = clamp.(guess, delta, 1-delta)
            solution = optimize(objective, [delta], [1 - delta], guess, Fminbox(LBFGS()), Optim.Options(show_trace=false))

            # best = solution.u
            best = solution.minimizer[1]
            # value = -solution.minimum

            return best
        end
    end 

    function get_special_centers(delta, alpha1, beta1, alpha2, beta2, p1)

        special_centers, descriptions = get_special_centers_from_single_distributions(delta, alpha1, beta1, alpha2, beta2, p1)

        # maximize surv prob in mixture dist
        special_center3 = maximize_survival_probability(delta, alpha1, beta1, alpha2, beta2, p1)
        description3 = "mean-maximizer"
        
        special_center_push!(special_centers, descriptions, special_center3, description3)

        # might have eq surv prob to candidate3 in mixture dist or may be local max
        
        if (alpha1, beta1) != (alpha2, beta2) && p1 != 1 && !isnothing(special_center3)
            special_center4 = 1 - special_center3
            description4 = "mean-maximizer-complement"

            special_center_push!(special_centers, descriptions, special_center4, description4)
        end

        return special_centers, descriptions
    end

    function get_educated_guesses(special_centers)
        return unique(special_centers)
    end

    function try_educated_guesses(centers::Vector{Float64}, partition::Vector{Vector{Int64}}, extinction_probability::Float64, candidates::Vector{Float64}, delta::Float64, alpha1::Number, beta1::Number, alpha2::Number, beta2::Number, p1::Number, args...)
        fecundity = length(centers)
        is_altered = true
        while is_altered
            is_altered = false
            for i in 1:fecundity
                new_centers = copy(centers)
                for el in candidates
                    new_centers[i] = el
                    new_extinction_probability = get_extinction_probability(new_centers, partition, delta, alpha1, beta1, alpha2, beta2, p1, fecundity, args...)
                    if new_extinction_probability < extinction_probability
                        centers = new_centers
                        extinction_probability = new_extinction_probability
                        is_altered = true
                    end
                end
            end
        end
        return centers, extinction_probability
    end

    ################ post-hoc analysis
    function is_close(val1, val2, tol = 1e-8)
        return abs(val1 - val2) < tol
    end

    function find_matches(center::Float64, special_centers, descriptions, tol)
        return [(el_cent, el_desc) for (el_cent, el_desc) in zip(special_centers, descriptions) if is_close(center, el_cent, tol)]
    end

    function find_matches(centers::Vector{Float64}, special_centers, descriptions, tol)
        return [find_matches(el, special_centers, descriptions, tol) for el in centers]
    end
    
     
    #################### helper functions
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

    ############### running everything
    # function minimize_extinction_probability(fecundity::Int64, delta::Float64, alpha1::Number, beta1::Number, alpha2::Number, beta2::Number, p1::Number, save_dir::String, population_size::Int64, partition_mutation_rate::Float64, use_educated_guess::Bool, analyze_centers::Bool, tol::Float64, args...)
        
    #     Random.seed!(0)

    #     idx2partition = get_idx2partition(fecundity)

    #     objective_function = make_objective_function(fecundity, delta, alpha1, beta1, alpha2, beta2, p1, idx2partition, args...)

    #     initial_population = initialize_population(population_size, fecundity, idx2partition)
    #     custom_differentiation = create_custom_differentiation(fecundity, partition_mutation_rate, idx2partition)

    #     lower_constraint = fill(delta, fecundity)
    #     upper_constraint = fill(1-delta, fecundity)
    #     lower_constraint = [lower_constraint; float(minimum(keys(idx2partition)))]
    #     upper_constraint = [upper_constraint; float(maximum(keys(idx2partition)))]
    #     constraints = CustomEvolutionary1.BoxConstraints(lower_constraint, upper_constraint)

    #     de_algorithm = CustomEvolutionary1.DE(populationSize = population_size, differentiation = custom_differentiation)

    #     results = CustomEvolutionary1.optimize(objective_function, constraints, de_algorithm, initial_population)

    #     # exctract results
    #     centers = results.minimizer[1:fecundity]
    #     partition = idx2partition[results.minimizer[end]]
    #     extinction_probability = results.minimum

    #     # educated guess true optimum
    #     if use_educated_guess
    #         candidates, descriptions = get_educated_guesses(delta, alpha1, beta1, alpha2, beta2, p1)
    #         centers, replaced_by, extinction_probability = try_educated_guesses(centers, partition, extinction_probability, candidates, descriptions, delta, alpha1, beta1, alpha2, beta2, p1, args...)
    #     end

    #     # gets mean fitness maximizer 
    #     # mean_maximizer = maximize_survival_prob(delta, alpha1, beta1, alpha2, beta2, p1)

    #     # output
    #     output = Dict("fecundity" => fecundity,
    #                 "delta" => delta,
    #                 "alpha1" => alpha1,
    #                 "beta1" => beta1,
    #                 "alpha2" => alpha2,
    #                 "beta2" => beta2,
    #                 "p1" => p1,
    #                 "centers" => centers,
    #                 "partition" => partition,
    #                 "extinction_probability" => extinction_probability
    #                 )

    #     if use_educated_guess
    #         output["replaced_by"] = replaced_by
    #         output["special_centers"] = hcat(candidates, descriptions)
    #     end
        
    #     save_results(output, save_dir)

    #     return output
    # end

    function minimize_partition_extinction_probability(partition::Vector{Vector{Int64}}, fecundity::Int64, delta::Float64, alpha1::Number, beta1::Number, alpha2::Number, beta2::Number, p1::Number, population_size::Int64, use_educated_guess::Bool, analyze_centers::Bool, tol::Float64, args...)

        # idx2partition = get_idx2partition(fecundity)

        objective_function = make_objective_function(partition, fecundity, delta, alpha1, beta1, alpha2, beta2, p1, args...)

        lower_constraint = fill(delta, fecundity)
        upper_constraint = fill(1-delta, fecundity)

        constraints = Evolutionary.BoxConstraints(lower_constraint, upper_constraint)

        de_algorithm = Evolutionary.DE(populationSize = population_size)

        # results = Evolutionary.optimize(objective_function, constraints, de_algorithm, Evolutionary.Options(iterations = 1))
        results = Evolutionary.optimize(objective_function, constraints, de_algorithm)

        # exctract results
        centers = results.minimizer[1:fecundity]
        extinction_probability = results.minimum

        if use_educated_guess || analyze_centers
            special_centers, descriptions = get_special_centers(delta, alpha1, beta1, alpha2, beta2, p1)
        end

        # educated guess true optimum
        if use_educated_guess
            # print("start educated guess", "\n")
            candidates = get_educated_guesses(special_centers)
            centers, extinction_probability = try_educated_guesses(centers, partition, extinction_probability, candidates, delta, alpha1, beta1, alpha2, beta2, p1, args...)
            # print("end educated guess", "\n")
        end

        if analyze_centers
            if "mean-maximizer" ∈ descriptions 
                index = findfirst(isequal("mean-maximizer"), descriptions)
                mean_maximizer = special_centers[index]
                mean_maximizer_matches = find_matches(mean_maximizer, special_centers, descriptions, tol)
            end
            if "mean-maximizer-complement" ∈ descriptions 
                index = findfirst(isequal("mean-maximizer-complement"), descriptions)
                mean_maximizer_complement = special_centers[index]
                mean_maximizer_complement_matches = find_matches(mean_maximizer_complement, special_centers, descriptions, tol)
            end
            center_matches = find_matches(centers, special_centers, descriptions, tol)
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

        if analyze_centers
            if "mean-maximizer" ∈ descriptions
                output["mean_maximizer_matches"] = mean_maximizer_matches
            end
            if "mean-maximizer-complement" ∈ descriptions
                output["mean_maximizer_complement_matches"] = mean_maximizer_complement_matches
            end
            output["center_matches"] = center_matches
        end

        return output
    end


    function partition2signature(partition)
        patch_sizes = [length(el) for el in partition]
        return sort(patch_sizes)
    end

    function get_partitions(fecundity)
        all_partitions = collect(partitions(1:fecundity))
        
        unique_partitions = []
        signatures = []
        for el in all_partitions
            signature =  partition2signature(el)
            if signature ∉ signatures
                push!(unique_partitions, el)
                push!(signatures, signature)
            end
        end
        return unique_partitions
    end

    function minimize_extinction_probability(fecundity::Int64, delta::Float64, alpha1::Number, beta1::Number, alpha2::Number, beta2::Number, p1::Number, save_dir::String, population_size::Int64, use_educated_guess::Bool, analyze_centers::Bool, tol::Float64, args...)

        Random.seed!(0)

        partitions = get_partitions(fecundity)

        extinction_probability = 10
        output = nothing
        
        for el in partitions
            # print("starting", "\n", el, "\n")
            current_output = minimize_partition_extinction_probability(el, fecundity, delta, alpha1, beta1, alpha2, beta2, p1, population_size, use_educated_guess, analyze_centers, tol, args...)
            
            if current_output["extinction_probability"] < extinction_probability
                extinction_probability = current_output["extinction_probability"]
                output = current_output
            end
            # print("finished", "\n", el, "\n")
        end

        # print("saving", "\n")
        save_results(output, save_dir)
        # print("saved", "\n")
        return output
    end

end
