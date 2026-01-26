"""
Transmission Module for handling bus data in a power system.
"""
module Transmission

# Use Julia standard libraries and third-party packages
using NamedArrays, JuMP, DataFrames, CSV

# Use internal modules
using ..Utils

# Export variables and functions
export Bus, Line, Load, load_data, stochastic_capex_model!

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
mutable struct Bus
    id:: Int64
    name:: String
    bus_type:: String
    max_flow_MW:: Union{Float64, Nothing}
end

function Bus(id, name, bus_type; max_flow_MW=nothing)
    return Bus(id, name, bus_type, max_flow_MW)
end

"""
Load represents the load demand at a specific bus, scenario, and timepoint.
# Fields:
- bus_id: ID of the bus
- sc_id: ID of the scenario
- t_id: ID of the timepoint
- load_MW: load demand in megawatts (MW)
"""
mutable struct Load
    bus:: String
    scenario:: String
    timepoint:: String
    load_MW:: Float64
    bus_id:: Int64
    scenario_id:: Int64
    timepoint_id:: Int64
end

function Load(bus, scenario, timepoint, load_MW; bus_id=0, scenario_id=0, timepoint_id=0)
    return Load(bus, scenario, timepoint, load_MW, bus_id, scenario_id, timepoint_id)
end

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
mutable struct Line
    id:: Int64
    name:: String
    from_bus:: String
    to_bus:: String
    cap_existing_power_MW:: Float64
    cost_fixed_power_USDperkW:: Float64
    r_pu:: Float64
    x_pu:: Float64
    g_pu:: Float64
    b_pu:: Float64
    angle_max_deg:: Float64
    angle_min_deg:: Float64
    expand_capacity:: Bool
    bus_id_from:: Int64
    bus_id_to:: Int64
    max_flow_MW:: Union{Float64, Nothing}
end

# Default values for Line
function Line(id, name, from_bus, to_bus, cap_existing_power_MW, r_pu, x_pu, g_pu, b_pu, angle_max_deg, angle_min_deg, expand_capacity; bus_id_from=0, bus_id_to=0, cost_fixed_power_USDperkW=0.0, max_flow_MW=nothing)
       return Line(id, name, from_bus, to_bus, cap_existing_power_MW, cost_fixed_power_USDperkW, r_pu, x_pu, g_pu, b_pu, angle_max_deg, angle_min_deg, expand_capacity, bus_id_from, bus_id_to, max_flow_MW)
end

"""
Load bus data from a CSV file and return it as a NamedArray of Bus structures.
"""
function load_data(inputs_dir:: String, S, T)

    start_time = time()
    filename = "buses.csv"
    println(" > $filename ...")
    N = to_structs(Bus, joinpath(inputs_dir, filename))
    println("   └ Completed, loaded ", length(N), " buses. Elapsed time ", round(time() - start_time, digits=2), " seconds.")

    start_time = time()
    filename = "lines.csv"
    println(" > $filename ...")
    L = to_structs(Line, joinpath(inputs_dir, filename))
    println("   └ Completed, loaded ", length(L), " lines. Elapsed time ", round(time() - start_time, digits=2), " seconds.")
    for l in L
        # Find ids of from_bus and to_bus from line instance
        l.bus_id_from = findfirst(n -> n.name == l.from_bus, N)
        l.bus_id_to = findfirst(n -> n.name == l.to_bus, N)
    end

    for n in N
        # Exit if the bus has not constraint on max power flow
        if n.max_flow_MW !== nothing
            continue
        end
        
        # Otherwise, we deduce max power flow based on the lines
        # the current bus is connected to.
        n.max_flow_MW = 0.0
        connected_lines = [line for line in L if (n.name in [line.from_bus, line.to_bus])]

        for line in connected_lines
            # There is no constraint on max power flow on the line,
            # thus the bus should also inherit no constraint (and we can exit).
            if (line.cap_existing_power_MW === nothing)
                n.max_flow_MW = nothing
                continue
            # Otherwise
            else
                n.max_flow_MW += line.cap_existing_power_MW
            end
        end
    end

    # Load load data
    start_time = time()
    filename = "loads.csv"
    println(" > $filename ...")
    load = to_structs(Load, joinpath(inputs_dir, filename); add_id_col = false)
    println("   └ Completed, loaded ", length(load), " load entries. Elapsed time ", round(time() - start_time, digits=2), " seconds.")

    for l in load
        # Find ids of bus, scenario, and timepoint from load instance
        if !isnothing(findfirst(n -> n.name == l.bus, N))
            l.bus_id = findfirst(n -> n.name == l.bus, N) 
        end

        if !isnothing(findfirst(s -> s.name == l.scenario, S))
            l.scenario_id = findfirst(s -> s.name == l.scenario, S)
        end
        
        
        if !isnothing(findfirst(t -> t.name == l.timepoint, T))
            l.timepoint_id = findfirst(t -> t.name == l.timepoint, T)
        end
    end

    
    # Transform load data into a multidimensional NamedArray
    #load = to_multidim_array(load, [:bus_id, :scenario_id, :timepoint_id], :load_MW)
    
    return N, L, load
end


"""
`build_admittance_matrix(N:: Vector{String}, lines:: Vector{Any}; include_shunts=false) 
                         :: Matrix{ComplexF64}`

This function builds the admittance matrix of any power system.

## Args:
    - buses: a vector containing the buses of the system. For example, buses=["san_diego", "lima"]
    - lines: A vector or list of instances of the structure Line. The struct Line must
             have the following attributes from_bus, to_bus, r, x, g, b. Note that g and b are the 
             conductance and susceptance, respectively, in one extreme of the line.

## Optional Args:
    - include_shunts: if yes, the conductance (g) and susceptance (b) are considered in the calculation of
                      the admittance matrix.

## Returns:
    - Y: 

    NamedArray(Y, (bus_names, bus_names), (:bus_name, :bus_name))
    bus_names = getfield.(N, :bus_name)
    a NamedArray that contains the admittance matrix. Y is commonly defined as a pure array, 
         but here we use a NamedArray, so the user can access entries of Y by two options:
         using strings like "san_diego", "lima", or numerical indices 1, 2 .. 
        for example: these combinations to access Y data work:
            Y["san_diego", "lima"]   = 0+0im
            Y["lima", "san_diego"]  = 0+0im
            Y["lima", "lima"] = 0+0im
            Y["lima", "lima"] =  0+0im 

TODO: add hint type to the lines argument. We may need to import the Line Struct.

"""
function build_admittance_matrix(N:: Vector{Bus}, L:: Vector{Line}; include_shunts=false):: Matrix{ComplexF64}

    # Define admittance matrix (actually it is NamedArray)
    # Note: we opt to use a NamedArray so N does not have to be a vector of numbers
    #       then, the user has more flexibility to access the admittance matrix, for example, Y["sandiego", "lima"]
    Y = zeros(Complex, length(N), length(N))
    
    for line in L
        # Calculate branch admittance and shunt admittance
        z_branch = complex(line.r_pu, line.x_pu)
        y_branch = 1.0 / z_branch
        y_shunt = complex(line.g_pu, line.b_pu)
        
        bus_id_from = line.bus_id_from
        bus_id_to = line.bus_id_to

        # Off-diagonal elements. Y_ij = -y_ij
        Y[bus_id_from, bus_id_to] -= y_branch
        Y[bus_id_to, bus_id_from] -= y_branch

        # Diagonal elements. Note: Y_ii = y_1i + y2i + ... + yii + ...
        y_at_bus = y_branch + (include_shunts ? y_shunt : 0)
        Y[bus_id_from, bus_id_from] += y_at_bus 
        Y[bus_id_to, bus_id_to] += y_at_bus

    end

    return Y
end


"""
    - maxFlow: a dictionary that contains maximum power transfer ber bus. For example:
        Dict{String, Float64} with 2 entries:
            "san_diego" => 500
            "lima"    => 1000

"""
function get_maxFlow(N:: Vector{Bus}, L:: Vector{Line}):: Vector{Float64}

    maxFlow =  zeros(Float64, length(N))

    for line in L
        # Find ids of from_bus and to_bus from line instance
        from_bus = findfirst(n -> n.bus == line.bus_from, N)
        to_bus = findfirst(n -> n.bus == line.bus_to, N)

        maxFlow[from_bus] += line.rating_MVA
        maxFlow[to_bus] += line.rating_MVA

    end

    return maxFlow
end

function stochastic_capex_model!(sys, mod:: Model, model_settings:: Dict)

    # Extract system data
    N = @views sys.N
    L = @views sys.L
    S = @views sys.S
    T = @views sys.T
    load = @views sys.load

    # Build admittance matrix
    Y = build_admittance_matrix(N, L)
    B = imag(Y) # take susceptance matrix

    # Get slack bus
    slack_bus_id = findfirst(n -> n.bus_type == "slack", N)
    slack_bus = N[slack_bus_id]

    # Define bus angle variables
    @variable(mod, vTHETA[N, S, T])

    # Fix bus angle of slack bus
    fix.(vTHETA[slack_bus, S, T], 0)

    println(" - Load data lookup")
    load_lookup = {(ld.bus, ld.scenario, ld.timepoint) => ld.load_MW for ld in load}
    load_buses = [ld.bus for ld in load if abs(ld.load_MW) != 0]
    N_load = [n for n in N if n.name in load_buses]

    if model_settings["load_shedding"] == true
        println(" - Load shedding variables")
        @variables(mod, begin
                vSHED[N_load, S, T] ≥ 0       
        end)
    end    

    # Extracting expressions from other submodules
    eGenAtBus = mod[:eGenAtBus]
    eNetDischargeAtBus = mod[:eNetDischargeAtBus]
    
    # DC Power flow transfered from each bus, assumes sbase = 100 MVA for all lines
    N_at_bus = Dict(n.id => [N[k.id] for k in N if ( (B[n.id, k.id] != 0) && (n.id != k.id))] for n in N) # dictionary to access bus neighors

    @expression(mod, eFlowAtBus[n ∈ N, s ∈ S, t ∈ T], 
                    100 * sum(B[n.id, m.id] * (vTHETA[n, s, t] - vTHETA[m, s, t]) for m in N_at_bus[n.id]))

    if model_settings["line_capacity_expansion"] == true

        println(" - Set of expandable and non-expandable lines")
        L_expandable = [l for l in L if ((l.expand_capacity == true) && (l.cap_existing_power_MW !== nothing))]
        L_nonexpandable = [l for l in L if ((l.expand_capacity == false) && (l.cap_existing_power_MW !== nothing))]

        println(" - Decision variables of line capacity expansion")
        @variable(mod, vCAPL[L_expandable] ≥ 0)

        println(" - Maximum and minimum flow constraints per expandable line")
        @constraint(mod, cMaxFlowPerLine[l ∈ L_expandable, s ∈ S, t ∈ T],
                100 * l.x_pu/(l.x_pu^2 + l.r_pu^2) *  (vTHETA[N[l.bus_id_from], s, t] - vTHETA[N[l.bus_id_to], s, t]) ≤ +l.cap_existing_power_MW + vCAPL[l] )

        @constraint(mod, cMinFlowPerLine[l ∈ L_expandable, s ∈ S, t ∈ T],
            100 * l.x_pu/(l.x_pu^2 + l.r_pu^2) *  (vTHETA[N[l.bus_id_from], s, t] - vTHETA[N[l.bus_id_to], s, t]) >= -(l.cap_existing_power_MW + vCAPL[l]) )

        println(" - Maximum flow constraints per non-expandable line")
        @constraint(mod, cFlowPerNonExpLine[l ∈ L_nonexpandable, s ∈ S, t ∈ T],
            -l.cap_existing_power_MW ≤ 100 * l.x_pu/(l.x_pu^2 + l.r_pu^2) *  (vTHETA[N[l.bus_id_from], s, t] - vTHETA[N[l.bus_id_to], s, t]) ≤ +l.cap_existing_power_MW )

        println(" - Line cost per period expression")
        @expression(mod, eLineCostPerPeriod,
                     sum(l.cost_fixed_power_USDperkW * vCAPL[l] * 1000 for l in L_expandable))
        
        eCostPerPeriod =  @views mod[:eCostPerPeriod]
        unregister(mod, :eCostPerPeriod)
        @expression(mod, eCostPerPeriod, eCostPerPeriod + eLineCostPerPeriod)

    end

    if (model_settings["line_capacity"] == true) && (model_settings["line_capacity_expansion"] == false)
        println(" - Maximum and minimum flow constraints per line")
        L_cap_constrained = [l for l in L if l.cap_existing_power_MW !== nothing]
        @constraint(mod, cFlowPerNonExpLine[l ∈ L_cap_constrained, s ∈ S, t ∈ T],
                           - l.cap_existing_power_MW ≤ 100 * l.x_pu/(l.x_pu^2 + l.r_pu^2) *  (vTHETA[N[l.bus_id_from], s, t] - vTHETA[N[l.bus_id_to], s, t]) ≤ + l.cap_existing_power_MW )

    end

    if model_settings["bus_max_flow"] == true

        buses_with_max_flow = [n for n in N if (n.max_flow_MW  !== nothing) ]
        @constraint(mod, cFlowPerBus[n ∈ buses_with_max_flow, s ∈ S, t ∈ T],
                           - n.max_flow_MW ≤ eFlowAtBus[n, s, t] ≤ n.max_flow_MW)
    end

    if model_settings["angle_difference_limits"] == true

        lines_with_angle_limits = [l for l in L if ((l.angle_min_deg > -360) && (l.angle_max_deg < 360))]
        @constraint(mod, cDiffAngle[l ∈ lines_with_angle_limits, s ∈ S, t ∈ T],
                           l.angle_min_deg * pi/180 ≤ (vTHETA[N[l.bus_id_from], s, t] - vTHETA[N[l.bus_id_to], s, t])  ≤ l.angle_max_deg * pi/180)
    end

    # Power balance at each bus
    @constraint(mod, cEnergyBalance[n ∈ N, s ∈ S, t ∈ T], 
                    (eGenAtBus[n, s, t]
                    + eNetDischargeAtBus[n, s, t] 
                    + ( ((model_settings["load_shedding"]) && (n.name in load_buses)) ? vSHED[n, s, t] : 0 )) * t.weight == 
                     (get(load_lookup, (n.name, s.name, t.name), 0.0) + eFlowAtBus[n, s, t]) * t.weight )

    if model_settings["load_shedding"] == true
        # Total load shedding cost expression
        println(" - Load shedding cost expressions")
        @expression(mod, eShedCostPerTp[t ∈ T],
                            1/length(S) * sum(s.probability * 5000 * vSHED[n, s, t] for n ∈ N_load, s ∈ S)) 
        @expression(mod, eShedTotalCost,
                                sum(eShedCostPerTp[t] * t.weight for t ∈ T))

        eCostPerTp =  @views mod[:eCostPerTp]
        unregister(mod, :eCostPerTp)
        @expression(mod, eCostPerTp[t ∈ T], eCostPerTp[t] + eShedCostPerTp[t])
    end

end

function toCSV_stochastic_capex(sys, mod:: Model, outputs_dir:: String)

    # Set of line instances
    L = @views sys.L 
    lines_df = DataFrame(L) #as dataframe

    # Export line capacities
    if haskey(mod, :vCAPL)
        to_df(mod[:vCAPL], [:line, :capacity_MW]; struct_fields=[:name], csv_dir = joinpath(outputs_dir,"line_built_capacity.csv"))
        
        costs = DataFrame( component = ["CostPerPeriod_USD"],
                                cost = [ value(mod[:eLineCostPerPeriod]) ])
        CSV.write(joinpath(outputs_dir, "line_costs_summary.csv"), costs)
    end

    # Print shedding if applicable
    if haskey(mod, :vSHED)
        to_df(mod[:vSHED], [:bus, :scenario, :timepoint, :load_shedding_MW]; csv_dir = joinpath(outputs_dir,"load_shedding.csv"), struct_fields=[:name, :name, :name])
    
        costs = DataFrame(component  = ["TotalCost_USD"],
                              cost  =  [value(mod[:eShedTotalCost])])
        CSV.write(joinpath(outputs_dir, "load_shedding_costs_summary.csv"), costs)
    end

    # Dataframe of bus angle. Columns: bus, scenario, timepoint, rad
    angle_df = to_df(mod[:vTHETA], [:bus, :scenario, :timepoint, :rad]; struct_fields=[:name, :name, :name])

    # Join
    df = rightjoin(lines_df, angle_df, on=[:from_bus => :bus])
    dropmissing!(df) # drop rows with missing values
    rename!(df, [:rad => :bus_from_angle]) # rename col

    # Join
    df = rightjoin(df, angle_df, on=[:to_bus => :bus, :scenario, :timepoint])
    dropmissing!(df) # drop rows with missing values
    rename!(df, [:rad => :bus_to_angle]) # rename col

    df.y_pu = 1 ./ (df.r_pu + 1im * df.x_pu) # compute branch admittance
    df.pflow_MW = -100*(df.bus_from_angle - df.bus_to_angle) .* imag.(df.y_pu) # compute DC power flowing in the branch, assumes Sbase = 100 MVA
    df.fict_loss_MW = df.r_pu .* abs.(df.pflow_MW) .* abs.(df.pflow_MW) # compute losses (not considered in the optimization)
    df.bus_from_angle = df.bus_from_angle * 180/pi # transform to deg
    df.bus_to_angle = df.bus_to_angle * 180/pi # transform to deg
    select!(df, Not(:y_pu, :r_pu, :x_pu, :g_pu, :b_pu))
    
    #filename = "line_flows.csv"
    #CSV.write(joinpath(outputs_dir, filename), df)
    #println("   - $filename printed.")
end
   
end # module Transmission