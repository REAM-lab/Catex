"""
Buses Module for handling bus data in a power system.
"""
module Buses

# Use Julia standard libraries and third-party packages
using NamedArrays, JuMP

# Use internal modules
using ..Utils

# Export variables and functions
export Bus, load_data, stoch_model!

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
    tp_id:: String
    load:: Float64
end

"""
Load bus data from a CSV file and return it as a NamedArray of Bus structures.
"""
function load_data(inputs_dir:: String):: Tuple{NamedArray{Bus}, NamedArray{Union{Missing, Float64}}}

    # Get a list of Bus structures
    buses = to_Structs(Bus, inputs_dir, "buses.csv")

    # Get a list of the bus IDs
    N = getfield.(buses, :bus_id)

    # Transform buses into NamedArray, so we can access buses by their IDs
    buses = NamedArray(buses, (N))

    # Load load data
    l = to_Structs(Load, inputs_dir, "loads.csv")
    
    # Transform load data into a multidimensional NamedArray
    load = to_multidim_NamedArray(l, [:bus_id, :sc_id, :tp_id], :load)

    return buses, load
end


function stochastic_capex_model!(sys, mod:: Model)

    N = sys.N
    S = sys.S
    T = sys.T

    # Get slack bus
    slack_bus = N[ findfirst([n.slack == true for n in N]) ]

    # Define bus angle variables
    @variable(mod, THETA[N, S, T]) 

    # Fix bus angle of slack bus
    fix.(THETA[slack_bus, S, T], 0)
    
end
   


end # module Buses