# packages
using ArgParse
include("extinction_probability_multipatch_functions.jl")
using .ExtinctionMultipatch#: initialize_population, get_idx2partition, make_objective_function, create_custom_differentiation
    
# println(names(ExtinctionMultipatch))
function setup_minimize_parser()
    settings = ArgParseSettings()
    @add_arg_table settings begin
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
        "--save-dir"
            arg_type = String
            required = true
        "--population-size"
            arg_type = Int64
            required = false
            default = 100
        "--partition-mutation-rate"
            arg_type = Float64
            required = false
            default = 0.2
        "--use-educated-guess"
            # arg_type = Bool
            action = :store_true
            # required = false
            # default = false

    end
    return settings
end

function setup_brute_force_minimize_parser()
    settings = ArgParseSettings()
    @add_arg_table settings begin
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
        "--save-dir"
            arg_type = String
            required = true
        "--population-size"
            arg_type = Int64
            required = false
            default = 100
        "--use-educated-guess"
            # arg_type = Bool
            action = :store_true
            # required = false
            # default = false
            

    end
    return settings
end

# function setup_extinction_probability_parser()
#     settings = ArgParseSettings()
#     @add_arg_table settings begin
#         "--centers"
#             arg_type = Vector{Float64}
#             required = true
#         "--partition"
#             arg_type = Vector{Vector{Int64}}
#             required = true
#         "--fecundity"
#             arg_type = Int64
#             required = true
#         "--delta"
#             arg_type = Float64
#             required = true
#         "--alpha1"
#             arg_type = Float64
#             required = true
#         "--beta1"
#             arg_type = Float64
#             required = true
#         "--alpha2"
#             arg_type = Float64
#             required = true
#         "--beta2"
#             arg_type = Float64
#             required = true
#         "--p1"
#             arg_type = Float64
#             required = true
#     end
#     return settings
# end

function setup_parser()
    settings = ArgParseSettings()

    @add_arg_table! settings begin
        "minimize"
            action = :command
        "brute-force-minimize"
            action = :command
    end

    return settings
end

function parse()
    main_settings = setup_parser()

    main_settings["minimize"] = setup_minimize_parser()
    main_settings["brute-force-minimize"] = setup_brute_force_minimize_parser()

    args = parse_args(main_settings)
    return args
end

function main()

    all_args = parse()

    command = all_args["%COMMAND%"]

    if command == "minimize"
        args = all_args["minimize"]
        return minimize_extinction_probability(args["fecundity"], args["delta"], args["alpha1"], args["beta1"],args["alpha2"], args["beta2"], args["p1"], args["save-dir"], args["population-size"], args["partition-mutation-rate"], args["use-educated-guess"])

    elseif command == "brute-force-minimize"
        args = all_args["brute-force-minimize"]
        return minimize_extinction_probability(args["fecundity"], args["delta"], args["alpha1"], args["beta1"], args["alpha2"], args["beta2"], args["p1"], args["save-dir"], args["population-size"], args["use-educated-guess"])
    end
end

# This conditional ensures that the script runs main only if it is not being included as a module
if abspath(PROGRAM_FILE) == @__FILE__
    # main()
    main()
end








