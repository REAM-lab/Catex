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
include("utils.jl")
include("scenarios.jl")
include("transmission.jl")
include("generators.jl")
include("storage.jl")
include("timescales.jl")
include("policies.jl")


# Use internal modules
using .Utils, .Scenarios, .Transmission, .Generators, .Storage, .Timescales, .Policies

# Export the functions we want users to be able to access easily
export init_system, solve_stochastic_capex_model, run_stocapex
export System, Scenario, Bus, Load, Generator, CapacityFactor, Line, StorageUnit, Timepoint, Timeseries, Policy

"""
System represents the entire power system for the stochastic capacity expansion problem.
# Fields:
- S: Vector containing instances of Scenario structure
- N: Vector containing instances of Bus structure
- loads: multidimensional array of load data
- G: Vector containing instances of Generator structure
- cf: multidimensional array of capacity factors data
- L: Vector containing instances of Line structure
- E: Vector containing instances of StorageUnit structure
- T: Vector containing instances of Timepoint structure
"""
struct System
    S:: Vector{Scenario}
    T:: Vector{Timepoint}
    TS:: Vector{Timeseries}
    N:: Vector{Bus}
    load:: Vector{Load}
    G:: Vector{Generator}
    cf:: NamedArray{Union{Missing, Float64}}
    L:: Vector{Line}
    E:: Vector{StorageUnit}
    policies:: Policy
end

"""
This function defines how to display the System struct in the REPL or when printed in Julia console.
"""
function Base.show(io::IO, ::MIME"text/plain", sys::System)
    println(io, "CATEX system:")
    println(io, "   ├ N (buses) = [", getfield.(sys.N, :bus)[1], ", ... ,", getfield.(sys.N, :bus)[end], "]")
    println(io, "   ├ L (lines) = [", getfield.(sys.L, :line)[1], ", ... ,", getfield.(sys.L, :line)[end], "]")
    println(io, "   ├ G (generators) = [", getfield.(sys.G, :generator)[1], ", ... ,", getfield.(sys.G, :generator)[end], "]")
    println(io, "   ├ E (storage) = [", getfield.(sys.E, :storage)[1], ", ... ,", getfield.(sys.E, :storage)[end], "]")
    println(io, "   ├ S (scenarios) = ", getfield.(sys.S, :scenario))
    println(io, "   ├ T (timepoints) = [", getfield.(sys.T, :timepoint)[1], ", ... ,", getfield.(sys.T, :timepoint)[end], "]")
    println(io, "   ├ TS (timeseries) = ", getfield.(sys.TS, :timeseries))
    println(io, "   ├ Policies: ", fieldnames(Policy))
    println(io, "   ├ Loads.")
    println(io, "   └ Capacity factors.")
end

"""
Initialize the System struct by loading data from CSV files in the inputs directory.
"""
function init_system(main_dir::String)

    println("-------------------------") 
    println(" Catex  - version 0.1.0") 
    println("-------------------------") 

    # Define the inputs directory
    full_start_time = time() 

    inputs_dir = joinpath(main_dir, "inputs")

    println("> Loading system data from $inputs_dir :")
    
    S = Scenarios.load_data(inputs_dir)

    T, TS = Timescales.load_data(inputs_dir)

    N, L, load = Transmission.load_data(inputs_dir, S, T)

    G, cf = Generators.load_data(inputs_dir)

    E = Storage.load_data(inputs_dir)

    policies = Policies.load_data(joinpath(main_dir, "inputs"))
     
    # Create instance of System struct
    sys = System(S, T, TS, N, load, G, cf, L, E, policies)

    total_time = round(time() - full_start_time; digits=2)
    println("> Total: $total_time seconds.\n")

    return sys
end

"""
Solves a stochastic capacity expansion problem.
""" 
function solve_stochastic_capex_model(sys, model_settings, main_dir, solver, solver_settings, print_model)


    println("> Building JuMP model:")
    full_start_time = time() 
    # Create JuMP model
    mod = Model(optimizer_with_attributes(solver, solver_settings...))

    # Initialize Costs for a period
    @expression(mod, eCostPerPeriod, 0)

    # Initialize Costs for a timepoint
    @expression(mod, eCostPerTp[t ∈ sys.T], 0)

    println(" > Generator vars and constraints ... ")
    tep = @elapsed Generators.stochastic_capex_model!(sys, mod,  model_settings)
    println(" └ Completed [$(round(tep, digits = 3)) seconds].")

    println(" > Storage vars and constraints ... ")
    tep = @elapsed Storage.stochastic_capex_model!(sys, mod, model_settings)
    println(" └ Completed [$(round(tep, digits = 3)) seconds].")

    println(" > Transmission vars and constraints ... ")
    tep = @elapsed Transmission.stochastic_capex_model!(sys, mod, model_settings)
    println(" └ Completed [$(round(tep, digits = 3)) seconds].")

    println(" > Policy vars and constraints ... ")
    tep = @elapsed Policies.stochastic_capex_model!(sys, mod)
    println(" └ Completed [$(round(tep, digits = 3)) seconds].")

    println(" > Objective function ... ")
    tep = @elapsed @expression(mod, eTotalCost, sum(mod[:eCostPerTp][t]*t.weight for t in sys.T) 
                                                    + mod[:eCostPerPeriod])
    println(" └ Completed [$(round(tep, digits = 3)) seconds].")

    @expression(mod, rescaling_factor_obj, 1e-6)  # Placeholder for potential rescaling
    @objective(mod, Min, eTotalCost*rescaling_factor_obj)

    # Print model to a text file if print_model==true. 
    # By default, it is print_model is false.
    # Useful for debugging purposes.
    if print_model
        filename = "model.txt"
        println(" > $filename printed")
        open(joinpath(main_dir, "outputs", filename), "w") do f
            println(f, mod)
        end
    end
    total_time = round(time() - full_start_time; digits=2)
    println("> Construction of model completed [$total_time seconds].")

    start_time = time()
    println("> JuMP model completed. Starting optimization: ")                
    optimize!(mod)

    print("\n")

    mod_status = termination_status(mod)
    mod_obj = round(objective_value(mod)/rescaling_factor_obj; digits=3) 
    println("> Optimization status: $mod_status")
    println("> Objective function value: $mod_obj USD")
    println("> Optimization elapsed time: ", round(time() - start_time, digits=2), " seconds.\n")
    return mod

end

"""
Exports results of the stochastic capacity expansion model to CSV files.
"""
function print_stochastic_capex_results(sys, mod:: Model, main_dir) 

    # Define the outputs directory
    outputs_dir = joinpath(main_dir, "outputs")
    if !isdir(outputs_dir)
        mkdir(outputs_dir)
    end

    full_start_time = time()
    println("> Printing files in $outputs_dir")
    
    println(" > Generators results ... ")
    tep = @elapsed Generators.toCSV_stochastic_capex(sys, mod, outputs_dir)
    println(" └ Completed [$(round(tep, digits = 3)) seconds].")

    println(" > Storage results ... ")
    tep = @elapsed Storage.toCSV_stochastic_capex(sys, mod, outputs_dir)
    println(" └ Completed [$(round(tep, digits = 3)) seconds].")
    println(" > Transmission results ... ")
    tep = @elapsed Transmission.toCSV_stochastic_capex(sys, mod, outputs_dir)
    println(" └ Completed [$(round(tep, digits = 3)) seconds].")

    # Print cost expressions
    start_time = time()
    filename = "costs_summary.csv"
    costs =  DataFrame(component  = ["CostPerTimepoint", "CostPerPeriod", "TotalCost"], 
                            cost  = [   value(sum(t.weight * mod[:eCostPerTp][t] for t in sys.T)), 
                                        value(mod[:eCostPerPeriod]), 
                                        value(mod[:eTotalCost])]) 
    CSV.write(joinpath(outputs_dir, filename), costs)
    println(" > $filename printed [", round(time() - start_time, digits=2), " seconds.]")

    println(" Total printing time: ", round(time() - full_start_time, digits=2), " seconds.\n")

end

function run_stocapex(; main_dir = pwd(), 
                             solver = Mosek.Optimizer,
                             solver_settings = Dict(),
                             print_model = false, 
                             model_settings = nothing)
    
    if isnothing(model_settings)
        model_settings = Dict(
                "generator_type_costs" => "linear",
                "load_shedding" => false,
                "single_storage_injection" => false,
                "line_capacity_expansion" => true,
                "line_capacity" => true,
                "bus_max_flow" => false,
                "angle_difference_limits" => true,
                "policies" => []
        )
    end
    start_time = time()
    sys = init_system(main_dir)
    mod = solve_stochastic_capex_model(sys, model_settings, main_dir, solver, solver_settings, print_model)
    print_stochastic_capex_results(sys, mod, main_dir)
    println(">>> Run completed in ", round(time() - start_time, digits=2), " seconds.")

    return sys, mod
end

end # module CATEX
