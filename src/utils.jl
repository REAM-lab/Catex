module utils

using DataFrames, JuMP, CSV

export to_stacked_Dict, to_Dict, to_tupled_Dict, to_Df, build_admittance_matrix, to_Structs

function to_Structs(structure::DataType, inputs_dir:: String, filename:: String):: Vector{structure}
    file_dir = joinpath(inputs_dir, filename)
    struct_names = fieldnames(structure)
    struct_types = fieldtypes(structure)

    first_csvlines = CSV.File(joinpath(inputs_dir, filename); limit=1)
    csv_header = Tuple(propertynames(first_csvlines))

    @assert csv_header == struct_names """Column names of $filename does not match the fields of the structure $structure."""

    df = CSV.read(file_dir, DataFrame; types=Dict(zip(struct_names, struct_types)))
    cols = Tuple(df[!, col] for col in names(df))
    V = structure.(cols...)    
    
    return V
end

function to_multidim_NamedArray(structures:: Vector{T}, dims:: Vector{Symbol}, value:: Symbol):: NamedArray{Union{Missing, Float64}} where {T} 
  
    vals_in_dim = [unique(getfield.(structures, d)) for (i, d) in enumerate(dims)]
    
    arr = Array{Union{Missing, Float64}}(missing, length.(vals_in_dim)...)
    arr = NamedArray(arr, vals_in_dim, dims)

    for s in structures
        arr[getfield.(Ref(s), dims)...] = getfield(s, value)
    end

    return arr
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
function to_Df(var_name:: JuMP.Containers.DenseAxisArray, header:: Vector, outputs_dir:: String, filename:: String; print_csv = true)
    dir_file = joinpath(outputs_dir, filename)
    df = DataFrame(Containers.rowtable(value, var_name; header = header))
    if print_csv
        CSV.write(dir_file, df)
        println(" > $filename printed.")
    end
    return df
end


"""
`build_admittance_matrix(N:: Vector{String}, lines:: Vector{Any}; include_shunts=false) 
                         :: NamedArray{ComplexF64}`

This function builds the admittance matrix of any power system.

## Args:
    - buses: a vector containing the buses of the system. For example, buses=["san_diego", "lima"]
    - lines: A vector or list of instances of the structure Line. The struct Line must
             have the following attributes from_bus, to_bus, r, x, g, b. Note that g and b are the 
             conductance and susceptance, respectively, in one extreme of the line.

## Optional Args:
    - include_shunts: if yes, the conductance (g) and susceptance (b) are considered in the calculation of
                      the admittance matrix.

## Returns:
    - Y: a NamedArray that contains the admittance matrix. Y is commonly defined as a pure array, 
         but here we use a NamedArray, so the user can access entries of Y by two options:
         using strings like "san_diego", "lima", or numerical indices 1, 2 .. 
        for example: these combinations to access Y data work:
            Y["san_diego", "lima"]   = 0+0im
            Y["lima", "san_diego"]  = 0+0im
            Y["lima", "lima"] = 0+0im
            Y["lima", "lima"] =  0+0im 

TODO: add hint type to the lines argument. We may need to import the Line Struct.

"""
function build_admittance_matrix(N:: Vector{String}, lines; include_shunts=false) :: NamedArray{ComplexF64}

    # Define admittance matrix (actually it is NamedArray)
    # Note: we opt to use a NamedArray so N does not have to be a vector of numbers
    #       then, the user has more flexibility to access the admittance matrix, for example, Y["sandiego", "lima"]
    num_buses = length(N)

    Y =  NamedArray( zeros(Complex, num_buses, num_buses), (N, N), (:bus_id, :bus_id))
    
    for line in lines
        # Calculate branch admittance
        z = complex(line.r, line.x)
        y = 1.0 / z
        
        # Extract from_bus and to_bus from line instance
        from_bus = line.from_bus
        to_bus = line.to_bus

        # Off-diagonal elements. Y_ij = -y_ij
        Y[from_bus, to_bus] -= y
        Y[to_bus, from_bus] -= y

        # Diagonal elements. Note: Y_ii = y_1i + y2i + ... + yii + ...
        Y[from_bus, from_bus] += y
        Y[to_bus, to_bus] += y
    end

    if include_shunts
        for line in lines
            # Calculate shunt admittance 
            y_shunt = complex(line.g, line.b)

            # Extract bus 
            from_bus = line.from_bus
            to_bus = line.to_bus
            
            # Add shunt admittance to the current admittance matrix
            Y[from_bus, from_bus] += y_shunt
            Y[to_bus, to_bus] += y_shunt
        end
    end

    return Y
end

"""
    - maxFlow: a dictionary that contains maximum power transfer ber bus. For example:
        Dict{String, Float64} with 2 entries:
            "san_diego" => 500
            "lima"    => 1000

"""
function get_maxFlow(N:: Vector{String}, lines):: NamedArray{Float64}

    num_buses = length(N)
    maxFlow =  NamedArray( zeros(Float64, num_buses), (N), :bus_id )

    for line in lines
        # Extract from_bus and to_bus from line instance
        from_bus = line.from_bus
        to_bus = line.to_bus
        rate = line.rate
        
        maxFlow[from_bus] += rate
        maxFlow[to_bus] += rate

    end
    return maxFlow
end


end # ends utils module