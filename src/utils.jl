module Utils

# Use Julia standard libraries and third-party packages
using DataFrames, JuMP, CSV, NamedArrays

# Export variables and functions
export to_structs, to_multidim_array, to_df, convert_df_columns!


function to_structs(structure::DataType, file_dir:: String; add_id_col = true):: Vector{structure}
    
    df = CSV.read(file_dir, DataFrame)
    if add_id_col
        insertcols!(df, 1, :id => 1:nrow(df))
    end

    field_type_pairs = Dict(fieldnames(structure) .=> fieldtypes(structure))
    field_type_pairs = Dict((k, v) for (k,v) in field_type_pairs if k in propertynames(df))
    convert_df_columns!(df, field_type_pairs)
    
    cols = Tuple(df[!, col] for col in names(df))

    V = structure.(cols...)    

    return V
end

function to_immutable_structs(structure::DataType, file_dir:: String; add_id_col = true):: Vector{structure}
    struct_names = fieldnames(structure)
    struct_types = fieldtypes(structure)

    first_twolines = CSV.read(file_dir, DataFrame; limit=1)
    csv_header = Tuple(propertynames(first_twolines))
    
    if add_id_col
        csv_header = (:id, csv_header...) 
    end
    
    @assert csv_header == struct_names """Incorrect column names of $file_dir.
                                          Column names must be $struct_names"""

    df = CSV.read(file_dir, DataFrame; types=Dict(zip(struct_names, struct_types)), validate=false)

    if add_id_col
        insertcols!(df, 1, :id => 1:nrow(df))
    end

    cols = Tuple(df[!, col] for col in names(df))
    V = structure.(cols...)    

    return V
end

function to_multidim_array(structures:: Vector{T}, dims:: Vector{Symbol}, value:: Symbol):: NamedArray{Union{Missing, Float64}} where {T} 
    
    # Get unique values in each dimensions, for example if dims=[:gen], then vals_in_dim=[["sd", "lima"]]
    vals_in_dim = [unique(getfield.(structures, d)) for (_, d) in enumerate(dims)]
    
    arr = NamedArray(Array{Union{Missing, Float64}}(missing, length.(vals_in_dim)...), vals_in_dim, dims)

    for s in structures
        arr[getfield.(Ref(s), dims)...] = getfield(s, value)
    end

    return arr
end

function convert_df_columns!(df::DataFrame, col_type_pairs::Dict{Symbol, DataType})
    for (col_name, target_type) in col_type_pairs
        if target_type == Int && eltype(df[!, col_name]) <: AbstractString
            df[!, col_name] = parse.(target_type, df[!, col_name])
        elseif target_type == Float64 && eltype(df[!, col_name]) <: AbstractString
            df[!, col_name] = parse.(target_type, df[!, col_name])
        elseif target_type == String && eltype(df[!, col_name]) == Int64
            df[!, col_name] = string.(df[!, col_name])
        else
            df[!, col_name] = convert.(target_type, df[!, col_name])
        end
    end
    return df
end


"""
`to_Df(var_name:: JuMP.Containers.DenseAxisArray, header:: Vector, outputs_dir:: String, filename:: String; print_csv = true)`

This function returns a csv file with the numerical solution of a JuMP variable once the optimization has been finished.

## Args:
    - var_name: variable name defined in the JuMP model, e.g., GEN, CAP.
    - header:   a vector with headers for the dataframe. It should be consistent with the dimensions.
                For example, if it is GEN[G, TPS], then header is [:generation_project, :timepoint, :DispatchGen_MW]
    - dir_file: directory where the csv will be saved. It must include the name, for example: /Users/paul/Documents/CATSExpand/examples/toy_example1/dispatch.csv
## Returns:
    - A csv file in the specified directory.
"""
function to_df(var_name:: JuMP.Containers.DenseAxisArray, 
                df_colnames:: Vector{Symbol}
                ; struct_fields=df_colnames[1:end-1], csv_dir=nothing) :: DataFrame

    # Transform DenseAxisArray into DataFrame.
    # The argument 'value' is a reserved word in JuMP to get the numerical value of a variable.
    # The argument 'var_name' is the name of the variable defined in the JuMP model. For example, GEN, CAP.
    # The argument 'header' is a vector with the headers we want for the dataframe. It should be consistent with the dimensions.
    df = DataFrame(Containers.rowtable(value, var_name; header = df_colnames))

    # The following loop is to replace the columns of the df with the actual IDs.
    # If we omit this loop, the columns will contain the entire struct, e.g., Timepoint(tp_id="tp1", period_id="p1", weight=1.0)
    # Actually, it is useful to have the entire struct in each cell of the Dataframe, because we can access its bus_id, and then make another dataframe.
    cols = @view df_colnames[1:end-1]
    for (c, f) in collect(zip(cols, struct_fields))
        df[!, c] = getfield.(df[!, c], f)
    end

    # Print csv file in the directory csv_dir
    if !isnothing(csv_dir)
        CSV.write(csv_dir, df)
        println(" > $(basename(csv_dir)) printed.")
    end

    return df # return the dataframe 
end


"""
Rehashing and Resizing: When a Dict grows beyond its current allocated capacity, 
it needs to be resized, which often involves rehashing all existing key-value pairs 
and moving them to a larger memory location. This process can be computationally expensive, 
especially with a large number of elements
"""
function to_stacked_Dict(data:: DataFrame, key:: String, value:: String)
    col_key = data[!, Symbol(key)]
    col_value = data[!, Symbol(value)]

    tuples = collect(zip(col_key, col_value))

    stacked_dict = Dict()

    for (key, value) in tuples
        if haskey(stacked_dict, key)
            push!(stacked_dict[key], value)
        else
            stacked_dict[key] = [value]
        end
    end
    
    return stacked_dict
end

function to_Dict(data:: DataFrame, key:: Symbol, value:: Symbol)
    return Dict(Pair.(data[:, key], data[:, value]))
end

function to_tupled_Dict(data:: DataFrame, keys:: Vector{Symbol}, value:: Symbol)
    col_keys = Tuple.(eachrow(data[:, keys]))
    col_value = data[:, value]
    return Dict(Pair.(col_keys, col_value))
end




end # ends utils module