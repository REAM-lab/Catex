module CATSExpand

# Import Julia packages
using CSV, DataFrames, DataStructures, Tables, JuMP, MosekTools, NamedArrays

# Define internal modules
include("Utils.jl")
include("Scenarios.jl")
include("Buses.jl")
include("Generators.jl")
include("Lines.jl")
include("EnergyStorages.jl")
include("Timepoints.jl")
include("Policies.jl")

# Use internal modules
using .Utils, .Scenarios, .Buses, .Generators, .Lines, .EnergyStorages, .Timepoints, .Policies

# Export the functions we want users to be able to access easily
export initialize, System, Scenario, Bus, Load, Generator, CapacityFactor, Line, EnergyStorage, Timepoint, Policy

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
    sc:: NamedArray{Scenario}
    bus:: NamedArray{Bus}
    load:: NamedArray{Union{Missing, Float64}}
    gen:: NamedArray{Generator}
    cf:: NamedArray{Union{Missing, Float64}}
    line:: NamedArray{Line}
    es:: NamedArray{EnergyStorage}
    tp:: NamedArray{Timepoint}
end

"""
This function defines how to display the System struct in the REPL or when printed in Julia console.
"""
function Base.show(io::IO, ::MIME"text/plain", s::System)
    println(io, "System:")
    println(io, " Scenarios (sc) = ", names(s.sc, 1))
    #print(io, " y = ", p.y)
end


"""
Initialize the System struct by loading data from CSV files in the inputs directory.
"""
function init_system(;main_dir = pwd())

    # Define the inputs directory
    inputs_dir = joinpath(main_dir, "inputs")
    
    # Fill in the fields of the System struct with CSV data
    sc = Scenarios.load_data(inputs_dir)
    bus, load = Buses.load_data(inputs_dir)
    gen, cf = Generators.load_data(inputs_dir)
    line = Lines.load_data(inputs_dir)
    es = EnergyStorages.load_data(inputs_dir)
    tp = Timepoints.load_data(inputs_dir)

    # Create instance of System struct
    sys = System(sc, bus, load, gen, cf, line, es, tp)

    return sys
end

function init_policies(;main_dir = pwd())   

    # Define the inputs directory
    inputs_dir = joinpath(main_dir, "inputs")

    # Create instance of Policy struct
    pol = Policies.load_data(inputs_dir)

    return pol
end

"""
Solves a stochastic capacity expansion problem.
""" 
function stoch_capex(sys, pol    ;main_dir = pwd(), 
                                  solver = Mosek.Optimizer,
                                    print_model = false)
    
    # Initialize the system by loading data
    sys = init_system(main_dir = main_dir)

    # Get scenarios
    sc = sys.sc
    S = names(sc,1) # IDs

    # Get buses
    bus = sys.bus
    N = names(bus,1) # IDs

    # Get slack bus
    slack_bus = buses[ findfirst([n.slack == true for n in buses]) ]

    # Get load
    load = sys.load

    # Get generators
    gen = sys.gen
    G = names(gen, 1) # IDs

    # Get capacity_factors
    cf = sys.cf

    # Get lines
    line = sys.line

    Y = build_admittance_matrix(N, line)
    B = imag(Y) # take susceptance matrix
    maxFlow = get_maxFlow(N, line)

    # Get storage units
    es = sys.es
    E = names(es, 1) # IDs

    # Get timepoints
    tp = sys.tp
    T = names(tp, 1) # IDs


   """
    - GV is a list of generator IDs that has capacity factor profile.
    - gensv is a list of instances of generators with capacity factor profiles.
        Power generation and capacity of these generators are considered random variables 
        in the second-stage of the stochastic problem.
    """
    GV = names(cf, 1) # get IDs of generators with capacity factor profiles
    gensv = gens[GV] # get list of instances of generators with capacity factor profiles


    """
    - GN is a list generator IDs that does not have capacity factor profile.
    - gensn is a list of instances of generators without capacity factor profiles.
      The power generation and capacity are considered as part of 
      first-stage of the stochastic problem.
    """

    GN = setdiff(G, GV) # get IDs of generators using set difference.
    gensn = gens[GN] # get list of instances of generators without capacity factor profiles


    G_AT_BUS = NamedArray( [getfield.(filter(g -> g.bus_id == n, gens).array, :gen_id) for n in N], 
                               (N), :bus_id )

    GV_AT_BUS = intersect.(G_AT_BUS, fill(GV, length(N)))
    GN_AT_BUS = intersect.(G_AT_BUS, fill(GN, length(N)))

    println("> Building JuMP model:")

    # Create JuMP model
    mod = Model(optimizer_with_attributes(solver))
    
    # Decision variables   
    @variables(m, begin
        GEN[GN, TPS] ≥ 0       
        CAP[GN] ≥ 0 
        GENV[GV, S, TPS] ≥ 0
        CAPV[GV, S] ≥ 0    
        THETA[N, S, TPS] 
    end)

    print(" > Generation constraints ...")
    
    # Minimum capacity of generators
    @constraint(mod, cFixCapGenVar[g ∈ GV, s ∈ S], 
                    CAPV[g, s] ≥ gen[g].exist_cap)

    @constraint(mod, cFixCapGenNonVar[g ∈ GN], 
                    CAP[g] ≥ gen[g].exist_cap)

    # Maximum build capacity 
    @constraint(mod, cMaxCapNonVar[g ∈ GN], 
                    CAP[g] ≤ gen[g].cap_limit)

    #Main.@infiltrate
    @constraint(mod, cMaxCapVar[g ∈ GV, s ∈ S], 
                    CAPV[g, s] ≤ gen[g].cap_limit)
    
    # Maximum power generation
    @constraint(mod, cMaxGenNonVar[g ∈ GN, t ∈ T], 
                    GEN[g, t] ≤ CAP[g])

    @constraint(mod, cMaxGenVar[g ∈ GV, s ∈ S, t ∈ T], 
                    GENV[g, s, t] ≤ cf[g, s, t]*CAPV[g, s])

    print(" ok.\n")

    print(" > Bus constraints ...")
    # Fix bus angle of slack bus
    fix.(THETA[slack_bus, S, T], 0)

    # Maximum power transfered by bus
    @constraint(mod, cMaxFlowAtBus[n ∈ setdiff(N, slack_bus), s ∈ S, t ∈ T],
                    -θlim ≤ THETA[n, s, t] ≤ θlim)

    # Power generation by bus
    @expression(mod, eGenAtBus[n ∈ N, s ∈ S, t ∈ T], 
                    sum(GEN[g, t] for g ∈ GN_AT_BUS[n]) + sum(GENV[g, s, t] for g ∈ GV_AT_BUS[n]) )

    # DC Power flow transfered from a bus
    @expression(mod, eFlowAtBus[n ∈ N, s ∈ S, t ∈ T], 
                    sum(B[n,m] * (THETA[n, s, t] - THETA[m, s, t]) for m in N))

    # Maximum power transfered by bus
    @constraint(mod, cMaxFlowAtBus[n ∈ N, s ∈ S, t ∈ T],
                    -maxFlow[n] ≤ eFlowAtBus[n, s, t] ≤ maxFlow[n])

    # Power balance
    @constraint(mod, cGenBalance[n ∈ N, s ∈ S, t ∈ T], 
                    eGenAtBus[n, s, t] ≥ load[n, s, t] + eFlowAtBus[n, s, t])    

    print(" ok.\n")

    print(" > Objective function ...")

    # The weighted operational costs of running each generator
    @expression(mod, eVariableCosts[t ∈ T],
                    sum(tp[t] * gen_variable_om[g] * GEN[g, t] for g ∈ GN, t ∈ TPS_IN_TS[ts]) + 
                    1/num_scen*(sum(prob[s]*tp_duration_hrs[t]* gen_variable_om[g] * GENV[g, s, t] for g ∈ GV, s ∈ S, t ∈ TPS_IN_TS[ts]) ))
    #Main.@infiltrate
    # Fixed costs 
	@expression(mod, eFixedCosts,
                    sum(gen_inv_cost[g] * CAP[g] for g ∈ GN) + 1/num_scen*sum( (prob[s] * gen_inv_cost[g] * CAPV[g, s]) for g ∈ GV, s ∈ S ))

    # Total costs
    @expression(mod, eTotalCosts,
                    sum(eTimeSeriesVariableCosts[ts] for ts ∈ TS) + eFixedCosts)
    
    @objective(mod, Min, eTotalCosts)
    print(" ok.\n")

    println("> JuMP model completed. Starting optimization: ")
                    
    optimize!(mod)

    print("\n")

    m_status = termination_status(m)
    m_obj = round(value(eTotalCosts); digits=3) 
    println("> Optimization status: $m_status")
    println("> Objective function: $m_obj")

    return mod

end


#=
function stoch_capex( ; main_dir = pwd(), 
                             solver = Mosek.Optimizer,
                             print_model = false)
    

    inputs_dir = joinpath(main_dir, "inputs")

    """
    Timepoint represents a timepoint in the optimization horizon for the 
    stochastic capacity expansion problem.

    # Fields:
    - tp_id: ID of the timepoint
    - timestamp: timestamp of the timepoint
    - duration: duration of the timepoint in hours
    """
    struct Timepoint
        tp_id:: Int64
        timestamp:: Float64
        duration:: Float64
    end 

    # Get a list of Scenario structures
    tps = to_Structs(Timepoint, inputs_dir, "timepoints.csv")

    # Get a list of the timepoint IDs
    TPS = getfield.(tps, :tp_id)

    # Transform tps into NamedArray, so we can access timepoints by their IDs
    tps = NamedArray(tps, (TPS))


    """
    Bus represents a bus or node in the power system.

    # Fields:
    - bus_id: ID of the bus
    - kv: voltage level of the bus in kilovolts
    - type: type of the bus (e.g., Substation). It could be any string.
    - lat: latitude of the bus location
    - lon: longitude of the bus location
    - slack: boolean (true or false) indicating if the bus is a slack bus. 
             At least there must be one slack bus in the system.
    """
    struct Bus
        bus_id:: String
        kv:: Float64
        type:: String
        lat:: Float64
        lon:: Float64
        slack:: Bool
    end

    # Get a list of Bus structures
    buses = to_Structs(Bus, inputs_dir, "buses.csv")

    # Get a list of the bus IDs
    N = getfield.(buses, :bus_id)

    # Transform buses into NamedArray, so we can access buses by their IDs
    buses = NamedArray(buses, (N))

    # Get slack bus
    slack_bus = buses[ findfirst([n.slack == true for n in buses]) ]


    """
    Load represents the load demand at a specific bus, scenario, and timepoint.
    # Fields:
    - bus_id: ID of the bus
    - sc_id: ID of the scenario
    - t_id: ID of the timepoint
    - load: load demand in megawatts (MW)
    """
    struct Load
        bus_id:: String
        sc_id:: String
        t_id:: String
        load:: Float64
    end
   
    l = to_Structs(Load, inputs_dir, "loads.csv")
    load = to_multidim_NamedArray(l, [:bus_id, :sc_id, :t_id], :load)


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
        gen_id:: String
        gen_tech:: String
        bus_id:: String
        c2:: Float64
        c1:: Float64
        c0:: Float64
        invest_cost:: Float64
        exist_cap:: Float64
        cap_limit:: Float64
        var_om_cost:: Float64
    end

    # Get a list of instances of generators structures
    gens = to_Structs(Generator, inputs_dir, "generators.csv")

    # Get a list of the generator IDs
    G = getfield.(gens, :gen_id)

    # Transform gens into NamedArray, so we can access generators by their IDs
    gens = NamedArray(gens, (G))
    

    """
    CapacityFactor represents the capacity factor of a generator. 
    at a specific scenario and timepoint. 

    # Fields:
    - gen_id: ID of the generation project
    - tp_id: ID of the timepoint
    - sc_id: ID of the scenario
    - capacity_factor: capacity factor (between 0 and 1)
    """

    struct CapacityFactor
        gen_id:: String
        tp_id:: String
        sc_id:: String
        capacity_factor:: Float64
    end

    c = to_Structs(CapacityFactor, inputs_dir, "capacity_factors.csv")
    cf = to_multidim_NamedArray(c, [:gen_id, :tp_id, :sc_id], :capacity_factor)



    """
    - GV is a list of generator IDs that has capacity factor profile.
    - gensv is a list of instances of generators with capacity factor profiles.
        Power generation and capacity of these generators are considered random variables 
        in the second-stage of the stochastic problem.
    """
    GV = names(cf, 1) # get IDs of generators with capacity factor profiles
    gensv = gens[GV] # get list of instances of generators with capacity factor profiles


    """
    - GN is a list generator IDs that does not have capacity factor profile.
    - gensn is a list of instances of generators without capacity factor profiles.
      The power generation and capacity are considered as part of 
      first-stage of the stochastic problem.
    """

    GN = setdiff(G, GV) # get IDs of generators using set difference.
    gensn = gens[GN] # get list of instances of generators without capacity factor profiles


    G_AT_BUS = NamedArray( [getfield.(filter(g -> g.bus_id == n, gens).array, :gen_id) for n in N], 
                               (N), :bus_id )

    GV_AT_BUS = intersect.(G_AT_BUS, fill(GV, length(N)))
    GN_AT_BUS = intersect.(G_AT_BUS, fill(GN, length(N)))

    
    """
    Line is a π-model transmission line connecting two buses in the power system.
    # Fields:
    - line_id: ID of the line
    - from_bus: ID of the bus where the line starts
    - to_bus: ID of the bus where the line ends
    - rate: thermal rating of the line (MW)
    - r: resistance of the line (p.u.)
    - x: reactance of the line (p.u.)
    - g: conductance of the shunt at one extreme of the line (p.u.)
    - b: susceptance of the shunt at one extreme of the line (p.u.)
    """
    struct Line
        line_id:: String
        from_bus:: String
        to_bus:: String
        rate:: Float64
        r:: Float64
        x:: Float64
        g:: Float64
        b:: Float64
    end

    lines = to_Structs(Line, inputs_dir, "lines.csv")

    Y = build_admittance_matrix(N, lines)

    B = imag(Y) # take susceptance matrix

    maxFlow = get_maxFlow(lines, N)

    """
    Energy storage represents an energy storage system in the power system.
    # Fields:
    - es_id: ID of the storage system
    - es_tech: technology type of the storage system, for example, "battery", "pumped_hydro". It could be any string.
    - bus_id: ID of the bus where the storage system is connected to.
    - invest_cost: investment cost per MW of power capacity (USD/MW)
    - exist_power_cap: pre-existing power capacity of the storage system (MW)
    - exist_energy_cap: pre-existing energy capacity of the storage system (MWh)
    - var_om_cost: variable operation and maintenance cost (USD/MW)
    - efficiency: round-trip efficiency of the storage system (between 0 and 1)
    - duration: duration of the storage system at full power (hours)
    """
    struct EnergyStorage 
        es_id:: String
        es_tech:: String
        bus_id:: String
        invest_cost:: Float64
        exist_power_cap:: Float64
        exist_energy_cap:: Float64
        var_om_cost:: Float64
        efficiency:: Float64
        duration:: Float64
    end

    # Get a list of instances of EnergyStorage structures
    ess = to_Structs(EnergyStorage, inputs_dir, "storages.csv")

    # Get a list of the storage IDs
    E = getfield.(gens, :gen_id)

    # Transform gens into NamedArray, so we can access storages by their IDs
    ess = NamedArray(ess, (E))
    

    # Additional settings
    max_diffangle = CSV.read(joinpath(inputs_dir, "max_diffangle.csv"),DataFrame;
                            header=["deg"],
                            types=[Float64])

    θlim = max_diffangle[1, :deg] * π/180

    
    println("> Building JuMP model:")

    # Create JuMP model
    m = Model(optimizer_with_attributes(solver))
    
    # Decision variables   
    @variables(m, begin
        GEN[GN, TPS] ≥ 0       
        CAP[GN] ≥ 0 
        GENV[GV, S, TPS] ≥ 0
        CAPV[GV, S] ≥ 0    
        THETA[N, S, TPS] 
    end)

    print(" > Generation constraints ...")
    
    # Minimum capacity of generators
    @constraint(m, cFixCapGenVar[g ∈ GV, s ∈ S], 
                    CAPV[g, s] ≥ gens[g].exist_cap)

    @constraint(m, cFixCapGenNonVar[g ∈ GN], 
                    CAP[g] ≥ gens[g].exist_cap)

    # Maximum build capacity 
    @constraint(m, cMaxCapNonVar[g ∈ GN], 
                    CAP[g] ≤ gens[g].cap_limit)

    #Main.@infiltrate
    @constraint(m, cMaxCapVar[g ∈ GV, s ∈ S], 
                    CAPV[g, s] ≤ gens[g].cap_limit)
    
    # Maximum power generation
    @constraint(m, cMaxGenNonVar[g ∈ GN, t ∈ TPS], 
                    GEN[g, t] ≤ CAP[g])

    @constraint(m, cMaxGenVar[g ∈ GV, s ∈ S, t ∈ TPS], 
                    GENV[g, s, t] ≤ cf[g, s, t]*CAPV[g, s])

    print(" ok.\n")

    print(" > Bus constraints ...")
    # Fix bus angle of slack bus
    fix.(THETA[slack_bus, S, T], 0)

    # Maximum power transfered by bus
    @constraint(m, cMaxFlowAtBus[n ∈ setdiff(N, slack_bus), s ∈ S, t ∈ TPS],
                    -θlim ≤ THETA[n, s, t] ≤ θlim)

    # Power generation by bus
    @expression(m, eGenAtBus[n ∈ N, s ∈ S, t ∈ TPS], 
                    sum(GEN[g, t] for g ∈ GN_AT_BUS[n]) + sum(GENV[g, s, t] for g ∈ GV_AT_BUS[n]) )

    # DC Power flow transfered from a bus
    @expression(m, eFlowAtBus[n ∈ N, s ∈ S, t ∈ TPS], 
                    sum(B[n,m] * (THETA[n, s, t] - THETA[m, s, t]) for m in N))

    # Maximum power transfered by bus
    @constraint(m, cMaxFlowAtBus[n ∈ N, s ∈ S, t ∈ TPS],
                    -maxFlow[n] ≤ eFlowAtBus[n, s, t] ≤ maxFlow[n])

    # Power balance
    @constraint(m, cGenBalance[n ∈ N, s ∈ S, t ∈ TPS], 
                    eGenAtBus[n, s, t] ≥ load[n, s, t] + eFlowAtBus[n, s, t])    

    print(" ok.\n")

    print(" > Objective function ...")

    # The weighted operational costs of running each generator
    @expression(m, eTimeSeriesVariableCosts[ts ∈ TS],
                    sum(tp_duration_hrs[t] * gen_variable_om[g] * GEN[g, t] for g ∈ GN, t ∈ TPS_IN_TS[ts]) + 
                    1/num_scen*(sum(prob[s]*tp_duration_hrs[t]* gen_variable_om[g] * GENV[g, s, t] for g ∈ GV, s ∈ S, t ∈ TPS_IN_TS[ts]) ))
    #Main.@infiltrate
    # Fixed costs 
	@expression(m, eFixedCosts,
                    sum(gen_inv_cost[g] * CAP[g] for g ∈ GN) + 1/num_scen*sum( (prob[s] * gen_inv_cost[g] * CAPV[g, s]) for g ∈ GV, s ∈ S ))

    # Total costs
    @expression(m, eTotalCosts,
                    sum(eTimeSeriesVariableCosts[ts] for ts ∈ TS) + eFixedCosts)
    
    @objective(m, Min, eTotalCosts)
    print(" ok.\n")

    println("> JuMP model completed. Starting optimization: ")
                    
    optimize!(m)

    print("\n")

    m_status = termination_status(m)
    m_obj = round(value(eTotalCosts); digits=3) 
    println("> Optimization status: $m_status")
    println("> Objective function: $m_obj")

    outputs_dir = joinpath(main_dir, "outputs")
    
    println("> Printing files in $outputs_dir")
    to_Df(GEN, [:generation_project, :timepoint, :DispatchGen_MW], outputs_dir , "dispatch.csv")
    to_Df(CAP, [:generation_project, :GenCapacity], outputs_dir , "gen_cap.csv")
    to_Df(GENV, [:generation_project, :scenario, :timepoint, :DispatchGen_MW], outputs_dir , "v_dispatch.csv")
    to_Df(CAPV, [:generation_project, :scenario, :GenCapacity], outputs_dir , "v_gen_cap.csv")
    
    filename = "costs_itemized.csv"
    varcost_data = to_Df(eTimeSeriesVariableCosts, [:timeseries, :cost], outputs_dir , "variable_costs.csv"; print_csv=false)
    insertcols!(varcost_data, 1, :component .=> "variable_costs")
    more_costs =  DataFrame(component  = ["fixed_costs", "total_costs"], 
                            timeseries = ["all", "all"], 
                            cost = [value(eFixedCosts), value(eTotalCosts)]) 
    costs_itemized = vcat(varcost_data, more_costs)
    CSV.write(joinpath(outputs_dir, filename), costs_itemized)
    println(" > $filename printed.")
      
    if print_model
        filename = "model.txt"
        println(" > $filename printed")
        open(joinpath(outputs_dir, filename), "w") do f
            println(f, m)
        end
    end

    #df = DataFrame(value_type.(model[name]).data, :auto)
    #CSV.write(joinpath(outputs_dir, "generation.csv"), df)

    print("Completed")
    
end 
=#

end # module CATSExpand
