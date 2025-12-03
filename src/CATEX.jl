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
- S: Vector containing instances of Scenario structure
- N: Vector containing instances of Bus structure
- loads: multidimensional array of load data
- G: Vector containing instances of Generator structure
- cf: multidimensional array of capacity factors data
- L: Vector containing instances of Line structure
- E: Vector containing instances of EnergyStorage structure
- T: Vector containing instances of Timepoint structure
"""
struct System
    S:: Vector{Scenario}
    T:: Vector{Timepoint}
    N:: Vector{Bus}
    load:: NamedArray{Union{Missing, Float64}}
    G:: Vector{Generator}
    cf:: NamedArray{Union{Missing, Float64}}
    L:: Vector{Line}
    E:: Vector{EnergyStorage}
end

"""
This function defines how to display the System struct in the REPL or when printed in Julia console.
"""
function Base.show(io::IO, ::MIME"text/plain", sys::System)
    println(io, "CATEX System:")
    println(io, "├ N (buses) = ", getfield.(sys.N, :name))
    println(io, "├ L (lines) = ", getfield.(sys.L, :name))
    println(io, "├ G (generators) = ", getfield.(sys.G, :name))
    println(io, "├ E (energy storages) = ", getfield.(sys.E, :name))
    println(io, "├ S (scenarios) = ", getfield.(sys.S, :name))
    println(io, "└ T (timepoints) = ", getfield.(sys.T, :name))
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
    S = to_structs(Scenario, joinpath(inputs_dir, "scenarios.csv"))
    T = to_structs(Timepoint, joinpath(inputs_dir, "timepoints.csv"))
    N = to_structs(Bus, joinpath(inputs_dir, "buses.csv"))
    L = to_structs(Line, joinpath(inputs_dir, "lines.csv"))
    G = to_structs(Generator, joinpath(inputs_dir, "generators.csv"))
    E = to_structs(EnergyStorage, joinpath(inputs_dir, "energy_storages.csv"))

    cf = process_cf(inputs_dir)
    load = process_load(inputs_dir)

    # Create instance of System struct
    sys = System(S, T, N, load, G, cf, L, E)

    println("ok.")
    return sys
end

function init_policies(;main_dir = pwd())   

    print("> Initializing policies data...")

    # Create instance of Policy struct
    pol = load_policies(joinpath(main_dir, "inputs"))

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

    # Under development ... 
    #print("> Energy storage vars and constraints ... ")
    #tep = @elapsed EnergyStorage.stochastic_capex_model!(mod, sys, pol)
    #println(" ok [$(round(tep, digits = 3)) seconds].")

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
