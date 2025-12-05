"""
Module for handling timepoints in the stochastic capacity expansion problem.
"""
module Timescales 
    
# Use Julia standard libraries and third-party packages
using NamedArrays

# Use internal modules
using ..Utils

# Export variables and functions
export Timepoint, Timeseries, load_data

"""
Timepoint represents a timepoint in the optimization horizon for the 
stochastic capacity expansion problem.
# Fields:
- tp_id: ID of the timepoint
- timestamp: timestamp of the timepoint
- duration: duration of the timepoint in hours
"""
mutable struct Timepoint
    id:: Int64
    name:: String
    timestamp:: Int64
    ts_name:: String
    ts_id:: Int64 
    weight:: Float64
    duration:: Float64
    prev_id:: Int64
end 

# Default values for Timepoint
function Timepoint(id, name, timestamp, ts_name; ts_id=0, weight=0, duration=0, prev_tp_id=0)
       return Timepoint(id, name, timestamp, ts_name, ts_id, weight, duration, prev_tp_id)
end

mutable struct Timeseries
    id:: Int64
    name:: String
    duration_of_tps:: Float64
    num_tps:: Int64
    ts_scale_to_period:: Float64
    tps_ids:: Vector{Int64}
end 

# Default values for Timeseries
function Timeseries(id, name, duration_of_tps, num_tps, ts_scale_to_period; tps_ids = Vector{Int64}())
       return Timeseries(id, name, duration_of_tps, num_tps, ts_scale_to_period, tps_ids)
end

function load_data(inputs_dir:: String)

    T = to_structs(Timepoint, joinpath(inputs_dir, "timepoints.csv"))
    TS = to_structs(Timeseries, joinpath(inputs_dir, "timeseries.csv"))

    for t in T
        # Timeseries id for each timepoint
        t.ts_id = findfirst(x -> x.name == t.ts_name, TS)

        # Duration (hrs) for each timepoint
        t.duration = (TS[t.ts_id]).duration_of_tps

        # Weight for each timepoint
        t.weight = t.duration * (TS[t.ts_id]).ts_scale_to_period

    end

    # For each timeseries, get the list (or vector) of ids of 
    # the timepoints that belong to each timeseries
    for ts in TS
        ts.tps_ids = findall(t -> t.ts_id == ts.id, T)
    end

    for t in T
        # Get ids of the timepoint that are within the timeseries to which the timepoint t belongs to
        tps_in_ts = TS[t.ts_id].tps_ids

        # If the id of the timepoint is the beginning timepoint of the timeseries,
        # then assign the end, otherwise just gives previous id
        if t.id == tps_in_ts[begin]
            t.prev_id = tps_in_ts[end]
        else
            t.prev_id = t.id - 1
        end
    end

    return T, TS
end
end # module Timepoints