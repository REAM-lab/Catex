module Policies 

# Use Julia standard libraries and third-party packages
using CSV, DataFrames

"""
Policy struct to hold policy parameters, additional restrictions, etc for the power system.
# Fields:
- budget: total budget available for investments
- bus_angle_diff: maximum allowable bus angle difference (in radians)
- max_CO2_emissions: maximum allowable CO2 emissions (in tons) (currently commented out)
"""
struct Policy
    # budget:: Float64
    max_diffangle:: Float64
    # max_CO2_emissions:: Float64
end

function load_data(inputs_dir:: String):: Policy

    # Read policies from CSV files
    # It is suggested to keep policies in different files as they can have different formats
    # or indices. For example, budget is a single value, while max CO2 emissions could be 
    # defined for certain time periods.

    max_diffangle = CSV.read(joinpath(inputs_dir, "max_diffangle.csv"), DataFrame;
                            header=["deg"], types=[Float64])

    max_diffangle = max_diffangle[1, :deg] * π/180


    return Policy(bus_angle_diff)
    # return Policy(budget, bus_angle_diff, max_CO2_emissions)
end
end # module Policies