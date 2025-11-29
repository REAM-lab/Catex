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
`build_admittance_matrix(buses:: Vector{String}, branch_data:: DataFrame; shunt_data=nothing)`

This function builds the admittance matrix of any power system.

## Args:
    - buses: a vector containing the buses of the system. For example, buses=["san_diego", "lima"]
    - branch_data (dataframe): It must include these columns:
            from_bus | to_bus | rate | r | l
        0   

    - shunt_data (dataframe): It must include these columns:
            bus | g | b
        0   
## Returns:
    - Y: a dictionary that contains admittance matrix information. Y is commonly defined as a matrix, but here we use
         a dictionary so the user can enter buses as a set of strings like "san_diego", "lima", instead of 1, 2 ..
        for example: 
        Dict{Tuple{String, String}, Complex{Float64}} with 4 entries:
            ("san_diego", "lima")   => 0+0im
            ("lima", "san_diego")  => 0+0im
            ("lima", "lima") => 0+0im
            ("san_diego", "san_diego")  => 0+0im 

    - maxFlow: a dictionary that contains maximum power transfer ber bus. For example:
        Dict{String, Float64} with 2 entries:
            "san_diego" => 500
            "lima"    => 1000

"""
function build_admittance_matrix(buses:: Vector{String}, branch_data:: DataFrame; shunt_data=nothing) 
                                :: Tuple{Dict{Tuple{String, String}, ComplexF64}, Dict{String, Float64}}

    # Define admittance matrix (actually it is dictionary here, but it acts as a matrix)
    # Note: we opt to use a dictionary so buses does not have to be a vector of numbers
    #       then, the user has more flexibility to enter a set of buses like "san_diego", "lima", instead of 1, 2 ..
    Y =  Dict((i, j) => 0.0 + 0.0im for i in buses for j in buses)
    maxFlow = Dict( buses .=> 0.0 )

    # Transform dataframe to a vector of tuples (for efficient iterations)
    branch_data = Tuple.(eachrow(branch_data)) # assume order (from_bus, to_bus, rate, r [p.u], x [p.u.])

    for (from_bus, to_bus, rate, r, x) in branch_data
        # Calculate branch admittance
        z = complex(r, x)
        y = 1.0 / z
        
        # Off-diagonal elements. Y_ij = -y_ij
        Y[(from_bus, to_bus)] -= y
        Y[(to_bus, from_bus)] -= y

        # Diagonal elements. Note: Y_ii = y_1i + y2i + ... + yii + ...
        Y[(from_bus, from_bus)] += y
        Y[(to_bus, to_bus)] += y

        maxFlow[from_bus] += rate
        maxFlow[to_bus] += rate
    end

    if shunt_data !== nothing
        shunt_data = Tuple.(eachrow(shunt_data)) # assume order (bus, g [p.u], b [p.u.])

        for (bus, g, b) in shunt_data

            # Calculate shunt admittance 
            y_shunt = complex(g, b)
            
            # Add shunt admittance to the current admittance matrix
            Y[(bus, bus)] += y_shunt

        end
    end
    return Y, maxFlow
end

end
