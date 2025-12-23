"""
Objective function value: 752076.861
"""

# Use Julia package
using CATEX, MosekTools
using Profile
# Set the main directory for the toy example
main_dir = @__DIR__

sys, mod = run_stocapex(;   main_dir = main_dir, 
                            solver = Mosek.Optimizer, 
                            solver_settings = Dict(),
                            print_model = false,
                            model_settings = Dict(  "gen_costs" => "quadratic", 
                                                    "consider_shedding" => false)) 

println("Finished")

