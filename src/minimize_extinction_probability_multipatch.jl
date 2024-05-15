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
        "--use_educated_guess"
            arg_type = Bool
            required = false
            default = false
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
        "--save_dir"
            arg_type = String
            required = true
        "--population_size"
            arg_type = Int64
            required = false
            default = 100
        "--use_educated_guess"
            arg_type = Bool
            required = false
            default = false
    end
    return settings
end

function setup_parser()
    settings = ArgParseSettings()

    @add_arg_table! settings begin
        "minimize"
            action = :command
        "brute_force_minimize"
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

# function main_old()

#     settings = ArgParseSettings()  # Create a settings object
#     subcommands = ArgParse.add_subcommands(s)

#     # function 1
#     minimize_cmd = ArgParse.add_subcommand(subcommands, "minimize", help="Minimize extinction probability")

#     @add_arg_table minimize_cmd begin
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
#         "--save_dir"
#             arg_type = String
#             required = true
#         "--population_size"
#             arg_type = Int64
#             required = false
#             default = 100
#         "--partition_mutation_rate"
#             arg_type = Float64
#             required = false
#             default = 0.2
#         "--use_educated_guess"
#             arg_type = Bool
#             required = false
#             default = false
#     end

#     # function 2
#     brute_force_minimize_cmd = ArgParse.add_subcommand(subcommands, "brute-force-minimize", help="Minimize extinction probability for each partition, and pick minimum")

#     @add_arg_table brute_force_minimize_cmd begin
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
#         "--save_dir"
#             arg_type = String
#             required = true
#         "--population_size"
#             arg_type = Int64
#             required = false
#             default = 100
#         "--use_educated_guess"
#             arg_type = Bool
#             required = false
#             default = false
#     end

#     args, selected_cmd = ArgParse.parse_args(settings)

#     if selected_cmd == "minimize"

#         return minimize_extinction_probability(args["fecundity"], args["delta"], args["alpha1"], args["beta1"],args["alpha2"], args["beta2"], args["p1"], args["save_dir"], args["population_size"], args["partition_mutation_rate"], args["use_educated_guess"])

#     elseif selected_cmd == "brute_force_minimize"
#         return minimize_extinction_probability(args["fecundity"], args["delta"], args["alpha1"], args["beta1"], args["alpha2"], args["beta2"], args["p1"], args["save_dir"], args["population_size"], args["use_educated_guess"])
#     end
# end

function main()

    all_args = parse()

    command = all_args["%COMMAND%"]

    if command == "minimize"
        args = all_args["minimize"]
        return minimize_extinction_probability(args["fecundity"], args["delta"], args["alpha1"], args["beta1"],args["alpha2"], args["beta2"], args["p1"], args["save_dir"], args["population_size"], args["partition_mutation_rate"], args["use_educated_guess"])

    elseif command == "brute-force-minimize"
        args = all_args["brute-force-minimize"]
        return minimize_extinction_probability(args["fecundity"], args["delta"], args["alpha1"], args["beta1"], args["alpha2"], args["beta2"], args["p1"], args["save_dir"], args["population_size"], args["use_educated_guess"])
    end
end

# This conditional ensures that the script runs main only if it is not being included as a module
if abspath(PROGRAM_FILE) == @__FILE__
    # main()
    main()
end








