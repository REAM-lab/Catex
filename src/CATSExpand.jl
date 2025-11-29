module CATSExpand

# Import Julia packages
using CSV, DataFrames, DataStructures, Tables, JuMP, MosekTools

include("utils.jl")

using .utils

# Export the functions you want users to be able to access easily
export stochastic_capex


"""
Solves a stochastic capacity expansion problem.
"""  
function stochastic_capex( ; main_dir = pwd(), 
                             solver = Mosek.Optimizer,
                             print_model = false)
    

    inputs_dir = joinpath(main_dir, "inputs")

    # Scenarios
    struct Scenario
        sc_id:: Int64
        name:: String
        prob:: Float64
    end

    #scen_data = CSV.read(joinpath(inputs_dir, "scenarios.csv"),DataFrame, types=[String, Float64]);

    #S = scen_data[:, :scenario]
    #num_scen = length(S)
    #prob = to_Dict(scen_data, :scenario, :probability)   

    S = to_Structs(Scenario, inputs_dir, "scenarios.csv")

    # Timepoints
    struct Timepoint
        t_id:: Int64
        timestamp:: Float64
        duration:: Float64
    end 

    T = to_Structs(Timepoint, inputs_dir, "timepoints.csv")


    #tps_data =  CSV.read(joinpath(inputs_dir, "timepoints.csv"), DataFrame, 
    #                        types=[Int64, Int64, String]);

    #ts_data = CSV.read(joinpath(inputs_dir, "timeseries.csv"), DataFrame, 
    #                        types=[String, Float64]);

    #TPS = tps_data[:, :timepoint]
    #TS = ts_data[:, :timeseries]

    #tps_ts_data = leftjoin(tps_data, ts_data, on = :timeseries)
    #tp_duration_hrs =  to_Dict(tps_ts_data, :timepoint, :tp_duration_hrs)   

    #TPS_IN_TS = to_stacked_Dict(tps_ts_data, "timeseries", "timepoint")

    # Buses and loads
    #buses_data =  CSV.read(joinpath(inputs_dir, "buses.csv"),DataFrame, 
    #                        types=[String, Bool]);
    # Buses 
    struct Bus
        bus_id:: Int64
        name:: String
        kv:: Float64
        type:: String
        lat:: Float64
        lon:: Float64
        slack:: Bool
    end

    N = to_Structs(Bus, inputs_dir, "buses.csv")


    loads_data = CSV.read(joinpath(inputs_dir, "loads.csv"),DataFrame, 
                            types=[String, String, Int64, Float64]);

    #N = buses_data[:, :bus]

    #load = to_tupled_Dict(loads_data, [:bus, :scenario, :timepoint], :bus_demand_mw)
  
    slack_bus = filter(row -> row.slack , buses_data) # get slack bus that has the column slack equal true
    slack_bus = slack_bus[1, :bus] # get the first bus 

    # Generation projects
    gens_data = CSV.read(joinpath(inputs_dir, "generation_projects_info.csv"), DataFrame, 
                                types=[String, String, String, Float64, Float64, Float64])

    # Capacity factors
    cf_data = CSV.read(joinpath(inputs_dir, "variable_capacity_factors.csv"),DataFrame,
                                types=[String, Int64, String, Float64])

    cf = to_tupled_Dict(cf_data, 
                                [:generation_project, :scenario, :timepoint], 
                                :gen_max_capacity_factor)
    
    G = gens_data[:, :generation_project]
    GV = unique(cf_data[:, :generation_project])
    GN = setdiff(G, GV)

    GENS_AT_BUS = to_stacked_Dict(gens_data, "bus", "generation_project")
    GV_AT_BUS = Dict(n => intersect(GV,gens) for (n, gens) in GENS_AT_BUS)
    GN_AT_BUS = Dict(n => intersect(GN,gens) for (n, gens) in GENS_AT_BUS)

    gen_variable_om= to_Dict(gens_data, :generation_project, :gen_variable_om)
    gen_capacity_limit= to_Dict(gens_data, :generation_project, :gen_capacity_limit)
    gen_min_build_capacity =  to_Dict(gens_data, :generation_project, :gen_min_build_capacity)

    # Build generation costs
    build_costs_data = CSV.read(joinpath(inputs_dir, "gen_build_costs.csv"),DataFrame,
                                types=[String, Float64])

    gens_build_costs = leftjoin(gens_data, build_costs_data, on = :gen_tech)

    gen_inv_cost = to_Dict(gens_build_costs, :generation_project, :investment_cost)
 
    
    # Pre-existing generation projects
    pregens_data = CSV.read(joinpath(inputs_dir, "gen_build_predetermined.csv"),DataFrame,
                            types=[String, Float64])

    PREGENS = pregens_data[:, :generation_project]

    PREGN = intersect(PREGENS, GN)
    PREGV = intersect(PREGENS, GV)
    
    precap = to_Dict(pregens_data, :generation_project, :gen_predetermined_cap)

    # Transmission lines
    trans_line_data = CSV.read(joinpath(inputs_dir, "transmission_lines.csv"),DataFrame;
                        types=[String, String, String, Float64, Float64, Float64, Float64])

    Y, maxFlow = build_admittance_matrix(N, trans_line_data)

    B = imag(Y) # take susceptance matrix

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

end # module CATSExpand
