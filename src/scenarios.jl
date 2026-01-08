"""
Scenarios module for managing scenario data in the stochastic capacity expansion problem.
"""
module Scenarios   

# Use Julia standard libraries and third-party packages


# Use internal modules
using ..Utils

# Export variables and functions
export Scenario, load_data

"""
Scenario represents a scenario in the stochastic capacity expansion problem.
# Fields:
- sc_id: ID of the scenario
- prob: probability of the scenario
"""
struct Scenario
    id:: Int64
    name:: String
    probability:: Float64
end

function load_data(inputs_dir:: String)
        
    # Load scenarios using CSV files
    filename = "scenarios.csv"
    print(" > $filename ...")
    S = to_structs(Scenario, joinpath(inputs_dir, filename))
    println(" ok, loaded ", length(S), " scenarios.")
    
    return S
end

end # module Scenarios