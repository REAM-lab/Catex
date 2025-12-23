"""
Decription of the case study:

It is the case5 from MATPOWER. The input data from MATPOWER has been re-formated and saved as csv files, 
which are the ones placed in the "inputs" folder. 

Objective function value: 17479.895
"""

# Use Julia package
using CATEX, MosekTools

# Set the main directory for the toy example
main_dir = @__DIR__

# Run stochastic capacity expansion 
sys, mod = run_stocapex(;   main_dir = main_dir, 
                            solver = Mosek.Optimizer, 
                            print_model = true,
                            gen_costs = "quadratic") 

#sys = init_system(main_dir = main_dir)

println("Finished")

