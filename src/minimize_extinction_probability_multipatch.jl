# packages
using ArgParse
include("extinction_probability_multipatch_functions.jl")
using .ExtinctionMultipatch

# using Pkg
# Pkg.precompile()

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
            action = :store_true
        "--analyze-centers"
            action = :store_true
        "--tol"
            arg_type = Float64
            required = false
            default = 1e-8
        "--q0"
            arg_type = Float64
            required = false
        "--num-generations"
            arg_type = Int64
            required = false
        "--num-runs"
            arg_type = Int64
            required = false
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
            action = :store_true
        "--analyze-centers"
            action = :store_true
        "--tol"
            arg_type = Float64
            required = false
            default = 1e-8
        "--q0"
            arg_type = Float64
            required = false
        "--num-generations"
            arg_type = Int64
            required = false
        "--num-runs"
            arg_type = Int64
            required = false
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
    args = all_args[command]
    if command == "minimize"
        if "q0" ∈ keys(args)
            return minimize_extinction_probability(args["fecundity"], args["delta"], args["alpha1"], args["beta1"],args["alpha2"], args["beta2"], args["p1"], args["save-dir"], args["population-size"], args["partition-mutation-rate"], args["use-educated-guess"], args["analyze-centers"], args["q0"], args["num-generations"], args["num-runs"])
        else
            return minimize_extinction_probability(args["fecundity"], args["delta"], args["alpha1"], args["beta1"],args["alpha2"], args["beta2"], args["p1"], args["save-dir"], args["population-size"], args["partition-mutation-rate"], args["use-educated-guess"], args["analyze-centers"])
        end
    elseif command == "brute-force-minimize"
        if "q0" ∈ keys(args)
            return minimize_extinction_probability(args["fecundity"], args["delta"], args["alpha1"], args["beta1"], args["alpha2"], args["beta2"], args["p1"], args["save-dir"], args["population-size"], args["use-educated-guess"], args["analyze-centers"], args["q0"], args["num-generations"], args["num-runs"])
        else
            return minimize_extinction_probability(args["fecundity"], args["delta"], args["alpha1"], args["beta1"], args["alpha2"], args["beta2"], args["p1"], args["save-dir"], args["population-size"], args["use-educated-guess"], args["analyze-centers"], args["tol"])
        end
    end
end

#     main()
# end

main()







