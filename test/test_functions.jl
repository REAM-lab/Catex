using CATEX, CSV
using CATEX.Utils
using DataFrames

# Pre-compile package again
Base.compilecache(Base.PkgId(CATEX))


# Set the main directory for the toy example
main_dir ="/Users/paul/Documents/CATEX/examples/toy_example1"
inputs_dir = joinpath(main_dir, "inputs")
#S = to_structs(Scenario, joinpath(inputs_dir, "scenarios.csv"))
#TS = to_structs(Timeseries, joinpath(inputs_dir, "timeseries.csv"))
#file_dir = joinpath(inputs_dir, "timepoints.csv")
#=
    df = CSV.read(file_dir, DataFrame)
    field_type_pairs = Dict(fieldnames(Timepoint) .=> fieldtypes(Timepoint))
    field_type_pairs = Dict((k, v) for (k,v) in field_type_pairs if k in propertynames(df))
    insertcols!(df, 1, :id => 1:nrow(df))
    convert_df_columns!(df, field_type_pairs)
    cols = Tuple(df[!, col] for col in names(df))
V = Timepoint.(cols...)   
=#
T = to_structs(Timepoint, joinpath(inputs_dir, "timepoints.csv"))
#df = DataFrame(name=["Alice", "Bob"], age=["25", "30"], score=[85, 92])
#convert_df_columns!(df, Dict(:name => String, :age => Int, :score => Float64))

#Dict(fieldnames(Scenario) .=> fieldtypes(Scenario))