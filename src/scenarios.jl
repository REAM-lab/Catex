"""
Scenarios module for managing scenario data in the stochastic capacity expansion problem.
"""
module Scenarios   

# Use Julia standard libraries and third-party packages
using NamedArrays

# Use internal modules
using ..utils

# Export variables and functions
export Scenario, load_data


"""
Scenario represents a scenario in the stochastic capacity expansion problem.
# Fields:
- sc_id: ID of the scenario
- prob: probability of the scenario
"""
struct Scenario
    sc_id:: String
    prob:: Float64
end

"""
Load scenario data from a CSV file and return it as a NamedArray of Scenario structures.
"""
function load_data(inputs_dir:: String)

    # Get a list of Scenario structures
    scen = to_Structs(Scenario, inputs_dir, "scenarios.csv")

    # Get a list of the scenario IDs
    S = getfield.(scen, :sc_id)

    # Transform scen into NamedArray, so we can access scenarios by their IDs
    scen = NamedArray(scen, (S))

    return scen

end

end # module Scenarios