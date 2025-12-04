"""
Generators module defines structure and functions for handling generator data in a power system.
"""
module Generators

# Use Julia standard libraries and third-party packages
using NamedArrays, JuMP, DataFrames, CSV

# Use internal modules
using ..Utils

# Export variables and functions
export Generator, process_cf, stochastic_capex_model!, toCSV_stochastic_capex

"""
Generator represents a generation project or existing generator in the power system.
# Fields:
- gen_id: ID of the generation project
- gen_tech: technology type of the generator, for example, "solar", "wind", "gas_cc". It could be any string.
- bus_id: ID of the bus where the generator is connected to.
- c2: quadratic coefficient of the generation cost function (USD/MW²)
- c1: linear coefficient of the generation cost function (USD/MW)
- c0: fixed coefficient of the generation cost function (USD)
- invest_cost: investment cost per MW of capacity (USD/MW)
- exist_cap: pre-existing capacity of the generator (MW)
- cap_limit: maximum build capacity of the generator (MW)
- var_om_cost: variable operation and maintenance cost (USD/MW)
"""
struct Generator
    id:: Int64
    name:: String
    tech:: String
    bus_name:: String
    c2:: Float64
    c1:: Float64
    c0:: Float64
    invest_cost:: Float64
    exist_cap:: Float64
    cap_limit:: Float64
    var_om_cost:: Float64
end

"""
CapacityFactor represents the capacity factor of a generator
at a specific scenario and timepoint. 
# Fields:
- gen_id: ID of the generation project
- tp_id: ID of the timepoint
- sc_id: ID of the scenario
- capacity_factor: capacity factor (between 0 and 1)
"""
struct CapacityFactor
    gen_name:: String
    sc_name:: String
    tp_name:: String
    capacity_factor:: Float64
end

"""
Load generator data from a CSV file and return it as a NamedArray of Generator structures.
"""
function process_cf(inputs_dir:: String) :: NamedArray{Union{Missing, Float64}}

    # Load capacity factor data
    cf = to_structs(CapacityFactor, joinpath(inputs_dir, "capacity_factors.csv"); add_id_col = false)

    # Transform capacity factor data into a multidimensional NamedArray
    cf = to_multidim_array(cf, [:gen_name, :sc_name, :tp_name], :capacity_factor)

    return cf
end

function stochastic_capex_model!(mod:: Model, sys, pol)

    S = @views sys.S
    T = @views sys.T
    G = @views sys.G
    cf = @views sys.cf
    N = @views sys.N

    """
    - GV is a vector of instances of generators with capacity factor profiles.
        Power generation and capacity of these generators are considered random variables 
        in the second-stage of the stochastic problem.
    """
    GV = filter(g -> g.name in names(cf, 1), G)

    """
    - GN is a vector of instances of generators without capacity factor profiles.
      The power generation and capacity are considered as part of 
      first-stage of the stochastic problem.
    """
    GN = setdiff( G, GV )

    G_AT_BUS = [filter(g -> g.bus_name == n.name, G) for n in N]
    GV_AT_BUS= intersect.(G_AT_BUS, fill(GV, length(N)))  
    GN_AT_BUS = intersect.(G_AT_BUS, fill(GN, length(N)))    
    
    # Define generation variables
    @variables(mod, begin
            vGEN[GN, T] ≥ 0       
            vCAP[GN] ≥ 0 
            vGENV[GV, S, T] ≥ 0
            vCAPV[GV, S] ≥ 0    
            end)

    # Minimum capacity of generators
    @constraint(mod, cMinCapGenVar[g ∈ GV, s ∈ S], 
                    vCAPV[g, s] ≥ g.exist_cap)

    @constraint(mod, cMinCapGenNonVar[g ∈ GN], 
                    vCAP[g] ≥ g.exist_cap)

    # Maximum build capacity 
    @constraint(mod, cMaxCapNonVar[g ∈ GN], 
                    vCAP[g] ≤ g.cap_limit)

    @constraint(mod, cMaxCapVar[g ∈ GV, s ∈ S], 
                    vCAPV[g, s] ≤ g.cap_limit)
    
    # Maximum power generation
    @constraint(mod, cMaxGenNonVar[g ∈ GN, t ∈ T], 
                    vGEN[g, t] ≤ vCAP[g])

    @constraint(mod, cMaxGenVar[g ∈ GV, s ∈ S, t ∈ T], 
                    vGENV[g, s, t] ≤ cf[g.name, s.name, t.name]*vCAPV[g, s])

    # Power generation by bus
    @expression(mod, eGenAtBus[n ∈ N, s ∈ S, t ∈ T], 
                    sum(vGEN[g, t] for g ∈ GN_AT_BUS[n.id]) 
                    + sum(vGENV[g, s, t] for g ∈ GV_AT_BUS[n.id]) )

    
    # The weighted operational costs of running each generator
    @expression(mod, eGenVariableCosts,
                    (sum(t.duration * g.var_om_cost * vGEN[g, t] 
                        + t.duration * g.c1 * vGEN[g, t] for g ∈ GN, t ∈ T) +
                    + 1/length(S)*(sum(s.prob * t.duration * g.c1 * vGENV[g, s, t] 
                                     + s.prob * t.duration * g.var_om_cost * vGENV[g, s, t] for g ∈ GV, s ∈ S, t ∈ T) ) ) )

    # Fixed costs 
	@expression(mod, eGenFixedCosts,
                    sum(g.invest_cost * vCAP[g] for g ∈ GN) 
                    + 1/length(S)*sum( (s.prob * g.invest_cost * vCAPV[g, s]) for g ∈ GV, s ∈ S ))

    # Total costs
    @expression(mod, eGenTotalCosts,
                        eGenVariableCosts + eGenFixedCosts)
end


function toCSV_stochastic_capex(sys, pol, mod:: Model, outputs_dir:: String)
    
    # Print vGEN variable solution
    to_df(mod[:vGEN], [:gen_name, :tp_name, :DispatchGen_MW]; 
                        struct_fields=[:name, :name], csv_dir = joinpath(outputs_dir,"gen_dispatch.csv"))
    
    # Print vCAP variable solution
    to_df(mod[:vCAP], [:gen_name, :GenCapacity]; 
                        struct_fields=[:name], csv_dir = joinpath(outputs_dir,"gen_cap.csv"))

    # Print vGENV variable solution
    to_df(mod[:vGENV], [:gen_name, :sc_name, :tp_name, :DispatchGen_MW]; 
                        struct_fields=[:name, :name, :name], csv_dir = joinpath(outputs_dir,"var_gen_dispatch.csv"))

    # Print vCAPV variable solution
    to_df(mod[:vCAPV], [:gen_name, :sc_name, :GenCapacity]; 
                        struct_fields=[:name, :name], csv_dir = joinpath(outputs_dir,"var_gen_cap.csv"))

    # Print cost expressions
    filename = "gen_costs_itemized.csv"
    costs =  DataFrame(component  = ["variable_costs", "fixed_costs", "total_costs"], 
                            cost  = [value(mod[:eGenVariableCosts]), value(mod[:eGenFixedCosts]), value(mod[:eGenTotalCosts])]) 
    CSV.write(joinpath(outputs_dir, filename), costs)
    println(" > $filename printed.")

end

end # module Generators