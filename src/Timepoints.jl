module Timepoints    
    
# Use Julia standard libraries and third-party packages
using NamedArrays

# Use internal modules
using ..utils

# Export variables and functions
export Timepoint, load_data

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

"""
Load timepoint data from a CSV file and return it as a NamedArray of Timepoint structures.
"""
function load_data(inputs_dir:: String)
    # Get a list of Scenario structures
    tps = to_Structs(Timepoint, inputs_dir, "timepoints.csv")

    # Get a list of the timepoint IDs
    TPS = getfield.(tps, :tp_id)

    # Transform tps into NamedArray, so we can access timepoints by their IDs
    tps = NamedArray(tps, (TPS))

    return tps
end

end # module Timepoints