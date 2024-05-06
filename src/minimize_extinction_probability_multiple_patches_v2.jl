# packages
using ArgParse
include("extinction_probability_multipatch_functions.jl")
using .ExtinctionMultipatch#: initialize_population, get_idx2partition, make_objective_function, create_custom_differentiation
    
println(names(ExtinctionMultipatch))

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

    results = optimize(objective_function, constraints, de_algorithm, initial_population)

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
    json_data = json(output)

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








