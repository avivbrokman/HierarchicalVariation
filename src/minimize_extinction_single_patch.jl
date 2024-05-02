# packages
using Intervals
using Distributions
using Polynomials
using Evolutionary
using LinearAlgebra #
using Graphs #


# 
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


function get_segment_probability(segment, distribution)
    return cdf(distribution, segment.last) - cdf(distribution, segment.first)
end

function get_segment_probabilities(segments, distribution)
    return [get_segment_probability(el, distribution) for el in segments]
end

function get_generating_coefficients(segment_survival_counts, segment_probabilities, fecundity)
    
    # probabilities = zeros(Float64, fecundity + 1)
    
    probabilities = Dict(i => 0.0 for i in 0:fecundity)

    # Iterate over the keys and values simultaneously
    for (key, value) in zip(segment_survival_counts, segment_probabilities)
        probabilities[key] += value
    end
    
    return probabilities
end



function get_extinction_coefficients_from_segments(segments, segment_survival_counts, alpha, beta, fecundity)
    dist1 = Beta(alpha, beta)

    segment_probabilities = get_segment_probabilities(segments, dist)
    
    generating_coefficients = get_generating_coefficients(segment_survival_counts, segment_probabilities,fecundity)

    generating_coefficients[1] -= 1

    return generating_coefficients
end

function get_extinction_coefficients(centers, delta, alpha, beta)
    fecundity = length(centers)
    
    intervals = [survival_interval(el, delta) for el in centers]

    segments, segment_survival = get_segment_survival_counts(intervals)

    coefficients = get_extinction_coefficients_from_segments(segments, segment_survival, alpha, beta, fecundity)

    return coefficients
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

function get_extinction_prob(centers, delta, alpha, beta)
    coefficients = get_extinction_coefficients(centers, delta, alpha, beta)
    extinction_prob = get_extinction_prob_from_coefficients(coefficients)
    return extinction_prob
end


# Workaround for the fact that Evolutionary.optimize doesn't allow for fixed parameters
struct GlobalParams
    delta::Float64
    alpha::Float64
    beta::Float64
    fecundity::Int8
end

function make_objective_function(params::GlobalParams)
    return centers -> get_extinction_prob(centers, params.delta, params.alpha, params.beta)
end


params = GlobalParams(0.1, 0.7, 0.3, 4)
objective_function = make_objective_function(params)

# optimization

# my mutation
# function gaussian(recombinant::AbstractVector, s::IsotropicStrategy;
#     rng::AbstractRNG = default_rng())
#     vals = randn(rng, length(recombinant)) * s.σ
#     recombinant += vals
#     return recombinant
# end


output = Evolutionary.optimize(objective_function, BoxConstraints([0,0,0,0], [1,1,1,1]), DE(populationSize = 100, F = 1.0))

# output = Evolutionary.optimize(objective_function, initial_centers, GA(populationSize = 100))
# output = Evolutionary.optimize(objective_function, BoxConstraints([0,0,0,0], [1,1,1,1]), GA(populationSize = 100))

# mutationRate: Probability of chromosome to be mutated
# ɛ/epsilon: Positive integer specifies how many individuals in the current generation are guaranteed to survive to the next generation. Floating number specifies fraction of population.
# selection: Selection function (default: tournament)
# crossover: Crossover function (default: genop)
# mutation: Mutation function (default: genop)
    
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