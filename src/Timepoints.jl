"""
Module for handling timepoints in the stochastic capacity expansion problem.
"""
module Timepoints    
    
# Use Julia standard libraries and third-party packages
using NamedArrays

# Use internal modules
using ..Utils

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
    tp_id:: String
    timestamp:: Int64
    duration:: Float64
end 

"""
Load timepoint data from a CSV file and return it as a NamedArray of Timepoint structures.
"""
function load_data(inputs_dir:: String):: NamedArray{Timepoint}

    # Get a list of Scenario structures
    T = to_Structs(Timepoint, inputs_dir, "timepoints.csv")

    # Get a list of the timepoint IDs
    Tids = getfield.(tps, :tp_id)

    # Transform T into NamedArray, so we can access timepoints by their IDs
    T = NamedArray(T, Tids, :tp_id)

    return T
end

function generate_prev_tps(T:: NamedArray{Timepoint}):: NamedArray{Timepoint}
    Tids = names(T, 1)

    prev = NamedArray{Timepoint}(undef, Tids, :tp_id)

    for (pos, t) in enumerate(T)
        if pos == 1
            prev[t.tp_id] = T[end]
        else
            prev[t.tp_id] = T[pos - 1]
        end
    end

    condition ? value_if_true : value_if_false 

    prev = NamedArray([pos == 1 ? T[end] : T[pos - 1] for (pos ,t) in enumerate(T)], Tids, :tp_id)
    return prev
end
end # module Timepoints