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
    timeseries:: String
    timeseries_id:: Int64 
    weight:: Float64
    duration_hr:: Float64
    prev_timepoint_id:: Int64
end 

# Default values for Timepoint
function Timepoint(id, name, timeseries; timeseries_id=0, weight=0, duration_hr=0, prev_timepoint_id=0)
       return Timepoint(id, name, timeseries, timeseries_id, weight, duration_hr, prev_timepoint_id)
end

mutable struct Timeseries
    id:: Int64
    name:: String
    timepoint_duration_hr:: Float64
    number_of_timepoints:: Int64
    timeseries_scale_to_period:: Float64
    timepoints_ids:: Vector{Int64}
end 

# Default values for Timeseries
function Timeseries(id, name, timepoint_duration_hr, number_of_timepoints, timeseries_scale_to_period; timepoints_ids = Vector{Int64}())
       return Timeseries(id, name, timepoint_duration_hr, number_of_timepoints, timeseries_scale_to_period, timepoints_ids)
end

function load_data(inputs_dir:: String)

    filename = "timepoints.csv"
    print(" > $filename ...")
    T = to_structs(Timepoint, joinpath(inputs_dir, filename))
    println(" ok, loaded ", length(T), " timepoints.")

    filename = "timeseries.csv"
    print(" > $filename ...")
    TS = to_structs(Timeseries, joinpath(inputs_dir, filename))
    println(" ok, loaded ", length(TS), " timeseries.")

    print(" > Timescale calculations ...")
    for t in T
        # Timeseries id for each timepoint
        t.timeseries_id = findfirst(x -> x.name == t.timeseries, TS)

        # Duration (hrs) for each timepoint
        t.duration_hr = (TS[t.timeseries_id]).timepoint_duration_hr 

        # Weight for each timepoint
        t.weight = t.duration_hr * (TS[t.timeseries_id]).timeseries_scale_to_period

    end

    # For each timeseries, get the list (or vector) of ids of 
    # the timepoints that belong to each timeseries
    for ts in TS
        ts.timepoints_ids = findall(t -> t.timeseries_id == ts.id, T)
    end

    for t in T
        # Get ids of the timepoint that are within the timeseries to which the timepoint t belongs to
        tps_in_ts = TS[t.timeseries_id].timepoints_ids

        # If the id of the timepoint is the beginning timepoint of the timeseries,
        # then assign the end, otherwise just gives previous id
        if t.id == tps_in_ts[begin]
            t.prev_timepoint_id = tps_in_ts[end]
        else
            t.prev_timepoint_id = t.id - 1
        end
    end
    println(" ok.")

    return T, TS
end
end # module Timescales