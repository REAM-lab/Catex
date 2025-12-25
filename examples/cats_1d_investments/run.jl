# Use Julia package
using Revise
using Infiltrator
using CATEX, MosekTools, Gurobi


# Set the main directory for the toy example
main_dir = @__DIR__

sys, mod = run_stocapex(;   main_dir = main_dir, 
                            solver = Gurobi.Optimizer, 
                            solver_settings = Dict("Crossover" => 0),
                            print_model = false,
                            model_settings = Dict(  "gen_costs" => "linear", 
                                                    "consider_shedding" => false)) 

#sys = init_system(main_dir)

println("Finished")

