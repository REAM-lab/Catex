"""
Transmission Module for handling bus data in a power system.
"""
module Transmission

# Use Julia standard libraries and third-party packages
using NamedArrays, JuMP, DataFrames, CSV

# Use internal modules
using ..Utils

# Export variables and functions
export Bus, Line, process_load, stochastic_capex_model!

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
    id:: Int64
    name:: String
    kv:: Float64
    type:: String
    lat:: Float64
    lon:: Float64
    slack:: Bool
end

"""
This function defines how to display the Bus struct in the REPL or when printed in Julia console.
"""
function Base.show(io::IO, bus::Bus)
    #print(io, "Bus(", bus.bus_id, ")")
end

"""
Load represents the load demand at a specific bus, scenario, and timepoint.
# Fields:
- bus_id: ID of the bus
- sc_id: ID of the scenario
- t_id: ID of the timepoint
- load: load demand in megawatts (MW)
"""
struct Load
    bus_name:: String
    sc_name:: String
    tp_name:: String
    load:: Float64
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
struct Line
    id:: Int64
    name:: String
    from_bus:: String
    to_bus:: String
    rate:: Float64
    r:: Float64
    x:: Float64
    g:: Float64
    b:: Float64
end

"""
Load bus data from a CSV file and return it as a NamedArray of Bus structures.
"""
function process_load(inputs_dir:: String):: NamedArray{Union{Missing, Float64}}

    # Load load data
    load = to_structs(Load, joinpath(inputs_dir,"loads.csv"); add_id_col = false)
    
    # Transform load data into a multidimensional NamedArray
    load = to_multidim_array(load, [:bus_name, :sc_name, :tp_name], :load)

    return load
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
        z_branch = complex(line.r, line.x)
        y_branch = 1.0 / z_branch
        y_shunt = complex(line.g, line.b)
        
        # Find ids of from_bus and to_bus from line instance
        from_bus = findfirst(n -> n.name == line.from_bus, N)
        to_bus = findfirst(n -> n.name == line.to_bus, N)

        # Off-diagonal elements. Y_ij = -y_ij
        Y[from_bus, to_bus] -= y_branch
        Y[to_bus, from_bus] -= y_branch

        # Diagonal elements. Note: Y_ii = y_1i + y2i + ... + yii + ...
        y_at_bus = y_branch + (include_shunts ? y_shunt : 0)
        Y[from_bus, from_bus] += y_at_bus 
        Y[to_bus, to_bus] += y_at_bus

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
        from_bus = findfirst(n -> n.name == line.from_bus, N)
        to_bus = findfirst(n -> n.name == line.to_bus, N)
        rate = line.rate
        
        maxFlow[from_bus] += rate
        maxFlow[to_bus] += rate

    end

    return maxFlow
end

function stochastic_capex_model!(mod:: Model, sys, pol)

    # Extract system data
    N = @views sys.N
    L = @views sys.L
    S = @views sys.S
    T = @views sys.T
    load = @views sys.load

    # Build admittance matrix and maxFlow
    Y = build_admittance_matrix(N, L)
    B = imag(Y) # take susceptance matrix
    maxFlow = get_maxFlow(N, L)

    # Get slack bus
    slack_bus_id = findfirst(n -> n.slack == true, N)
    slack_bus = N[slack_bus_id]

    # Define bus angle variables
    @variable(mod, vTHETA[N, S, T]) 

    # Fix bus angle of slack bus
    fix.(vTHETA[slack_bus, S, T], 0)

    # Extracting expressions from other submodules
    eGenAtBus = mod[:eGenAtBus]
    eNetDischargeAtBus = mod[:eNetDischargeAtBus]
    
    # DC Power flow transfered from each bus, assumes sbase = 100 MVA for all lines
    @expression(mod, eFlowAtBus[n ∈ N, s ∈ S, t ∈ T], 
                    100 * sum(B[n.id, m.id] * (vTHETA[n, s, t] - vTHETA[m, s, t]) for m in N))

    # Maximum power transfered at each bus
    @constraint(mod, cMaxFlowAtBus[n ∈ N, s ∈ S, t ∈ T],
                    -maxFlow[n.id] ≤ eFlowAtBus[n, s, t] ≤ maxFlow[n.id])

    # Power balance at each bus
    @constraint(mod, cGenBalance[n ∈ N, s ∈ S, t ∈ T], 
                    eGenAtBus[n, s, t] + eNetDischargeAtBus[n, s, t] ≥ 
                    load[n.name, s.name, t.name] + eFlowAtBus[n, s, t])    

end

function toCSV_stochastic_capex(sys, pol, mod:: Model, outputs_dir:: String)

    # Set of line instances
    L = @views sys.L 
    lines_df = DataFrame(L) #as dataframe

    # Dataframe of bus angle. Columns: bus_name, sc_name, tp_name, rad
    angle_df = to_df(mod[:vTHETA], [:bus_name, :sc_name, :tp_name, :rad]; struct_fields=[:name, :name, :name])

    # Join
    df = rightjoin(lines_df, angle_df, on=[:from_bus => :bus_name])
    dropmissing!(df) # drop rows with missing values
    rename!(df, [:rad => :from_bus_angle]) # rename col

    # Join
    df = rightjoin(df, angle_df, on=[:to_bus => :bus_name, :sc_name, :tp_name])
    dropmissing!(df) # drop rows with missing values
    rename!(df, [:rad => :to_bus_angle]) # rename col

    df.y = 1 ./ (df.r + 1im * df.x) # compute branch admittance
    df.pflow_MW = 100*(df.from_bus_angle - df.to_bus_angle) .* imag.(df.y) # compute DC power flowing in the branch, assumes Sbase = 100 MVA
    df.fict_loss_MW = df.r .* abs.(df.pflow_MW) .* abs.(df.pflow_MW) # compute losses (not considered in the optimization)
    select!(df, Not(:y, :r, :x, :g, :b))
    
    filename = "transmission_flows.csv"
    println(" > $filename printed.")
    CSV.write(joinpath(outputs_dir, filename), df)
end
   
end # module Transmission