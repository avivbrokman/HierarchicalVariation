# packages
using Intervals
using Distributions
using Polynomials
using Evolutionary


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



function get_extinction_coefficients_from_segments(segments, segment_survival_counts, alpha1, beta1, alpha2, beta2, p1, fecundity)
    dist1 = Beta(alpha1, beta1)
    dist2 = Beta(alpha2, beta2)

    segment_probabilities1 = get_segment_probabilities(segments, dist1)
    segment_probabilities2 = get_segment_probabilities(segments, dist2)

    segment_probabilities =  p1 .* segment_probabilities1 + (1 - p1) .* segment_probabilities2
    
    generating_coefficients = get_generating_coefficients(segment_survival_counts, segment_probabilities,fecundity)

    generating_coefficients[1] -= 1

    return generating_coefficients
end

function get_extinction_coefficients(centers, delta, alpha1, beta1, alpha2, beta2, p1)
    fecundity = length(centers)
    
    intervals = [survival_interval(el, delta) for el in centers]

    segments, segment_survival = get_segment_survival_counts(intervals)

    coefficients = get_extinction_coefficients_from_segments(segments, segment_survival, alpha1, beta1, alpha2, beta2, p1, fecundity)

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

function get_extinction_prob(centers, delta, alpha1, beta1, alpha2, beta2, p1)
    coefficients = get_extinction_coefficients(centers, delta, alpha1, beta1, alpha2, beta2, p1)
    extinction_prob = get_extinction_prob_from_coefficients(coefficients)
    return extinction_prob
end


# Workaround for the fact that Evolutionary.optimize doesn't allow for fixed parameters
struct GlobalParams
    delta::Float64
    alpha1::Float64
    beta1::Float64
    alpha2::Float64
    beta2::Float64
    p1::Float64
    fecundity::Int8
end

function make_objective_function(params::GlobalParams)
    return centers -> get_extinction_prob(centers, params.delta, params.alpha1, params.beta1, params.alpha2, params.beta2, params.p1)
end


params = GlobalParams(0.1, 0.7, 0.1, 0.1, 0.7, 0.5, 4)
objective_function = make_objective_function(params)

# optimization
uniform = Uniform()
initial_centers = rand(uniform, params.fecundity)
# initial_centers = [0.15, 0.85, 0.4, 0.7]
# initial_centers = [0.5, 0.5, 0.5, 0.5]

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
