# Use Julia package
using Revise
using Infiltrator
using CATEX, MosekTools, JuMP, DataFrames, CSV, NamedArrays


# Set the main directory for the toy example
main_dir ="/Users/paul/Documents/CATEX/examples/toy_example1"

#sys = init_system(main_dir= main_dir)
#pol = init_policies(main_dir= main_dir)
#model = solve_stochastic_capex_model(sys, pol, main_dir = main_dir)
sys, pol, mod = run_stocapex(; main_dir = main_dir, solver = Mosek.Optimizer, print_model = false) 

println("Finished running toy_example1")

