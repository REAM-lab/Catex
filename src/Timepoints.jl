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
    id:: Int64
    name:: String
    timestamp:: Int64
    duration:: Float64
end 


end # module Timepoints