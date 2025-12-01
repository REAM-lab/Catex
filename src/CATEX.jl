"""
CATEX.jl - California Transmission System Expansion Model
A Julia package for modeling capacity expansion and operational 
optimization. Originally developed for the state of California,
it can be adapted for other regions.
"""
module CATEX

# Import Julia packages
using CSV, DataFrames, JuMP, MosekTools, NamedArrays

# Define internal modules
include("Utils.jl")
include("Scenarios.jl")
include("Transmission.jl")
include("Generators.jl")
include("EnergyStorages.jl")
include("Timepoints.jl")
include("Policies.jl")


# Use internal modules
using .Utils, .Scenarios, .Transmission, .Generators, .EnergyStorages, .Timepoints, .Policies

# Export the functions we want users to be able to access easily
export init_system, init_policies, solve_stochastic_capex_model, run_stocapex
export System, Scenario, Bus, Load, Generator, CapacityFactor, Line, EnergyStorage, Timepoint, Policy

"""
System represents the entire power system for the stochastic capacity expansion problem.
# Fields:
- sc: NamedArray of instances of Scenario structure
- buses: NamedArray of instances of Bus structure
- loads: multidimensional NamedArray of load data
- gen: NamedArray of instances of Generator structure
- cf: multidimensional NamedArray of capacity factors data
- line: NamedArray of instances of Line structure
- es: NamedArray of instances of EnergyStorage structure
- tp: NamedArray of instances of Timepoint structure
"""
struct System
    S:: NamedArray{Scenario}
    N:: NamedArray{Bus}
    load:: NamedArray{Union{Missing, Float64}}
    G:: NamedArray{Generator}
    cf:: NamedArray{Union{Missing, Float64}}
    L:: NamedArray{Line}
    E:: NamedArray{EnergyStorage}
    T:: NamedArray{Timepoint}
end

"""
This function defines how to display the System struct in the REPL or when printed in Julia console.
"""
function Base.show(io::IO, ::MIME"text/plain", sys::System)
    println(io, "CATEX System:")
    println(io, "├ N (buses) = ", names(sys.N, 1))
    println(io, "├ L (lines) = ", names(sys.L, 1))
    println(io, "├ G (generators) = ", names(sys.G, 1))
    println(io, "├ E (energy storages) = ", names(sys.E, 1))
    println(io, "├ S (scenarios) = ", names(sys.S, 1))
    println(io, "└ T (timepoints) = ", names(sys.T, 1))
end

"""
Initialize the System struct by loading data from CSV files in the inputs directory.
"""
function init_system(;main_dir = pwd())

    println("-------------------------") 
    println(" CATEX  - version 0.1.0") 
    println("-------------------------") 

    print("> Initializing system data...")
    # Define the inputs directory
    inputs_dir = joinpath(main_dir, "inputs")
    
    # Fill in the fields of the System struct with CSV data
    scs = Scenarios.load_data(inputs_dir)
    buses, load, lines = Transmission.load_data(inputs_dir)
    gens, cf = Generators.load_data(inputs_dir)
    ess = EnergyStorages.load_data(inputs_dir)
    tps = Timepoints.load_data(inputs_dir)

    # Create instance of System struct
    sys = System(scs, buses, load, gens, cf, lines, ess, tps)

    println("ok.")
    return sys
end

function init_policies(;main_dir = pwd())   

    print("> Initializing policies data...")
    # Define the inputs directory
    inputs_dir = joinpath(main_dir, "inputs")

    # Create instance of Policy struct
    pol = Policies.load_data(inputs_dir)

    println("ok.")
    return pol
end

"""
Solves a stochastic capacity expansion problem.
""" 
function solve_stochastic_capex_model(sys, pol    ;main_dir = pwd(), 
                                    solver = Mosek.Optimizer,
                                    print_model = false)


    println("> Building JuMP model:")

    # Create JuMP model
    mod = Model(optimizer_with_attributes(solver))

    print("> Generator vars and constraints ... ")
    tep = @elapsed Generators.stochastic_capex_model!(mod, sys, pol)
    println(" ok [$(round(tep, digits = 3)) seconds].")

    print("> Transmission vars and constraints ... ")
    tep = @elapsed Transmission.stochastic_capex_model!(mod, sys, pol)
    println(" ok [$(round(tep, digits = 3)) seconds].")

    print("> Policy vars and constraints ... ")
    tep = @elapsed Policies.stochastic_capex_model!(mod, sys, pol)
    println(" ok [$(round(tep, digits = 3)) seconds].")

    @objective(mod, Min, mod[:eTotalCosts])

    # Print model to a text file if print_model==true. 
    # By default, it is print_model is false.
    # Useful for debugging purposes.
    if print_model
        filename = "model.txt"
        println(" > $filename printed")
        open(joinpath(main_dir, "outputs", filename), "w") do f
            println(f, m)
        end
    end

    println("> JuMP model completed. Starting optimization: ")
                    
    optimize!(mod)

    print("\n")

    mod_status = termination_status(mod)
    mod_obj = round(value(mod[:eTotalCosts]); digits=3) 
    println("> Optimization status: $mod_status")
    println("> Objective function value: $(round(mod_obj, digits=3))")
    return mod

end

"""
Exports results of the stochastic capacity expansion model to CSV files.
"""
function print_stochastic_capex_results(mod:: Model; main_dir = pwd()) 

    # Define the outputs directory
    outputs_dir = joinpath(main_dir, "outputs")

    println("> Printing files in $outputs_dir")
    
    Generators.toCSV_stochastic_capex(mod, outputs_dir)

end

function run_stocapex(; main_dir = pwd(), 
                             solver = Mosek.Optimizer,
                             print_model = false)
    
    sys = init_system(main_dir = main_dir)
    pol = init_policies(main_dir = main_dir)
    mod = solve_stochastic_capex_model(sys, pol; main_dir = main_dir, solver = solver)
    print_stochastic_capex_results(mod; main_dir = main_dir)

    return sys, pol, mod
end

end # module CATEX
