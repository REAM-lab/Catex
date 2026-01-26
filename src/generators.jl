"""
Generators module defines structure and functions for handling generator data in a power system.
"""
module Generators

# Use Julia standard libraries and third-party packages
using NamedArrays, JuMP, DataFrames, CSV

# Use internal modules
using ..Utils

# Export variables and functions
export Generator, load_data, stochastic_capex_model!, toCSV_stochastic_capex

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
mutable struct Generator
    id:: Int64
    name:: String
    technology:: String
    bus:: String
    site:: String
    cap_existing_power_MW:: Float64
    cap_max_power_MW:: Float64
    cost_fixed_power_USDperkW:: Float64
    cost_variable_USDperMWh:: Float64
    c0_USD:: Float64
    c1_USDperMWh:: Float64
    c2_USDperMWh2:: Float64
    expand_capacity:: Bool
end

function Generator(id, generator, technology, bus, site, cap_existing_power_MW, cap_max_power_MW, cost_fixed_power_USDperkW, cost_variable_USDperMWh, c0_USD, c1_USDperMWh, c2_USDperMWh2; expand_capacity = true)
    return Generator(id, generator, technology, bus, site, cap_existing_power_MW, cap_max_power_MW, cost_fixed_power_USDperkW, cost_variable_USDperMWh, c0_USD, c1_USDperMWh, c2_USDperMWh2, expand_capacity)
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
    site:: String
    scenario:: String
    timepoint:: String
    capacity_factor:: Float64
end

function load_data(inputs_dir::String)

    # Load generator data from CSV file
    start_time = time()
    filename = "generators.csv"
    println(" > $filename ...")
    G = to_structs(Generator, joinpath(inputs_dir, filename))
    println("   └ Completed, loaded ", length(G), " generators. Elapsed time ", round(time() - start_time, digits=2), " seconds.")

    # Load capacity factor data
    start_time = time()
    filename = "capacity_factors.csv"
    println(" > $filename ...")
    cf = to_structs(CapacityFactor, joinpath(inputs_dir, filename); add_id_col = false)
   
    # Transform capacity factor data into a multidimensional NamedArray
    cf = to_multidim_array(cf, [:site, :scenario, :timepoint], :capacity_factor; asNamedArray=true)
    println("   └ Completed, loaded ", length(cf), " capacity factor entries. Elapsed time ", round(time() - start_time, digits=2), " seconds.")

    # Extra calculations or checks can be added here
    for g in G
        if g.cap_existing_power_MW >= g.cap_max_power_MW
            g.expand_capacity = false
        end
    end

    return G, cf
end

function stochastic_capex_model!(sys, mod:: Model,  model_settings:: Dict)

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
    GV = filter(g -> g.site != "no_capacity_factor", G)

    """
    - GN is a vector of instances of generators without capacity factor profiles.
      The power generation and capacity are considered as part of 
      first-stage of the stochastic problem.
    """
    GN = setdiff( G, GV )

    G_AT_BUS = [filter(g -> g.bus == n.name, G) for n in N]
    GV_AT_BUS= intersect.(G_AT_BUS, fill(GV, length(N)))  
    GN_AT_BUS = intersect.(G_AT_BUS, fill(GN, length(N)))    
    
    # Define generation variables
    @variables(mod, begin
            vGEN[GN, T] ≥ 0       
            vCAP[GN] ≥ 0 
            vGENV[GV, S, T] ≥ 0
            vCAPV[GV, S] ≥ 0    
            end)
    
    for g in GN
        if !g.expand_capacity
            fix(vCAP[g], 0; force=true)
        end
    end

    for g in GV
        if !g.expand_capacity
            for s in S
                fix(vCAPV[g, s], 0; force=true)
            end
        end
    end

    # Maximum build capacity 
    println(" - Maximum built capacity constraints") 
    @constraint(mod, cMaxCapNonVar[g ∈ GN; g.expand_capacity == true], 
                    vCAP[g] ≤ g.cap_max_power_MW - g.cap_existing_power_MW)

    @constraint(mod, cMaxCapVar[g ∈ GV, s ∈ S; g.expand_capacity == true], 
                    vCAPV[g, s] ≤ g.cap_max_power_MW - g.cap_existing_power_MW)
    
    # Maximum power generation
    println(" - Maximum dispatch constraints")    
    @constraint(mod, cMaxGenNonVar[g ∈ GN, t ∈ T], 
                    vGEN[g, t] ≤ vCAP[g] + g.cap_existing_power_MW)

    println(" - Capacity factor constraints")
    @constraint(mod, cMaxGenVar[g ∈ GV, s ∈ S, t ∈ T], 
                    vGENV[g, s, t] ≤ cf[g.site, s.name, t.name] * (vCAPV[g, s] + g.cap_existing_power_MW))

    # Power generation by bus
    println(" - Dispatch at bus expressions")
    @expression(mod, eGenAtBus[n ∈ N, s ∈ S, t ∈ T], 
                    sum(vGEN[g, t] for g ∈ GN_AT_BUS[n.id]) 
                    + sum(vGENV[g, s, t] for g ∈ GV_AT_BUS[n.id]))

    
    # The weighted operational costs of running each generator
    if model_settings["generator_type_costs"] == "quadratic"
    @expression(mod, eGenCostPerTp[t ∈ T],
                        sum(g.c2_USDperMWh2 * vGEN[g, t]* vGEN[g, t] + g.c1_USDperMWh * vGEN[g, t] + g.c0_USD for g ∈ GN) + 
                        1/length(S) * sum(s.probability * (g.c2_USDperMWh2 * vGENV[g, s, t]* vGENV[g, s, t] + g.c1_USDperMWh * vGENV[g, s, t] + g.c0_USD ) for g ∈ GV, s ∈ S))

    elseif model_settings["generator_type_costs"] == "linear"
    @expression(mod, eGenCostPerTp[t ∈ T],
                        sum(g.cost_variable_USDperMWh * vGEN[g, t] for g ∈ GN) + 
                        1/length(S) * sum(s.probability * g.cost_variable_USDperMWh * vGENV[g, s, t] for g ∈ GV, s ∈ S))
    end
    
        
    eCostPerTp =  @views mod[:eCostPerTp]
    unregister(mod, :eCostPerTp)
    @expression(mod, eCostPerTp[t ∈ T], eCostPerTp[t] + eGenCostPerTp[t])

    # Fixed costs 
    println(" - Generation cost per period expression")
    @expression(mod, eGenCostPerPeriod,
                    sum(g.cost_fixed_power_USDperkW * vCAP[g] * 1000 for g ∈ GN) 
                    + 1/length(S) * sum( (s.probability * g.cost_fixed_power_USDperkW * vCAPV[g, s] * 1000) for g ∈ GV, s ∈ S ))

    eCostPerPeriod =  @views mod[:eCostPerPeriod]
    unregister(mod, :eCostPerPeriod)
    @expression(mod, eCostPerPeriod, eCostPerPeriod + eGenCostPerPeriod)

    @expression(mod, eGenTotalCost, sum(eGenCostPerTp[t] * t.weight for t in T) + eGenCostPerPeriod)

end


function toCSV_stochastic_capex(sys, mod:: Model, outputs_dir:: String)
    
    # Print vGEN variable solution
    to_df(mod[:vGEN], [:generator, :timepoint, :DispatchGen_MW]; csv_dir = joinpath(outputs_dir,"generator_dispatch.csv"), struct_fields=[:name, :name])
    
    # Print vCAP variable solution
    to_df(mod[:vCAP], [:generator, :GenCapacity]; csv_dir = joinpath(outputs_dir,"generator_capacity.csv"), struct_fields=[:name])

    # Print vGENV variable solution
    to_df(mod[:vGENV], [:generator, :scenario, :timepoint, :DispatchGen_MW]; csv_dir = joinpath(outputs_dir,"variable_generator_dispatch.csv"), struct_fields=[:name, :name, :name])

    # Print vCAPV variable solution
    to_df(mod[:vCAPV], [:generator, :scenario, :GenCapacity]; csv_dir = joinpath(outputs_dir,"variable_generator_capacity.csv"), struct_fields=[:name, :name, :name])


    # Print cost expressions
    filename = "generator_costs_summary.csv"
    costs =  DataFrame(component  = ["CostPerTimepoint_USD", "CostPerPeriod_USD", "TotalCost_USD"], 
                            cost  = [   value(sum(t.weight * mod[:eGenCostPerTp][t] for t in sys.T)), 
                                        value(mod[:eGenCostPerPeriod]), 
                                        value(mod[:eGenTotalCost])]) 
    CSV.write(joinpath(outputs_dir, filename), costs)
    println("   - $filename printed.")

end

end # module Generators