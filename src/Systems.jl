module Systems

# Use Julia standard libraries and third-party packages
using NamedArrays

# Use internal modules
using ..Scenarios: Scenario
using ..Buses: Bus, Load
using ..Generators: Generator, CapacityFactor
using ..Lines: Line
using ..EnergyStorages: EnergyStorage
using ..Timepoints: Timepoint

# Export variables and functions
export System

"""
System represents the entire power system for the stochastic capacity expansion problem.
# Fields:
- sc: NamedArray of instances of Scenario structure
- buses: NamedArray of instances of Bus structure
- loads: multidimensional NamedArray of load data
- gen: NamedArray of instances of Generator structure
- cf: multidimensional NamedArray of capacity factors data
- line: NamedArray of instances of Line structure
- es: NamedArray of instances of EnergyStorage structure
- tp: NamedArray of instances of Timepoint structure
"""
struct System
    S:: NamedArray{Scenario}
    N:: NamedArray{Bus}
    load:: NamedArray{Union{Missing, Float64}}
    G:: NamedArray{Generator}
    cf:: NamedArray{Union{Missing, Float64}}
    L:: NamedArray{Line}
    E:: NamedArray{EnergyStorage}
    T:: NamedArray{Timepoint}
end

"""
This function defines how to display the System struct in the REPL or when printed in Julia console.
"""
function Base.show(io::IO, ::MIME"text/plain", s::System)
    println(io, "System:")
    println(io, " Scenarios (sc) = ", names(s.sc, 1))
    #print(io, " y = ", p.y)
end

end # module Systems