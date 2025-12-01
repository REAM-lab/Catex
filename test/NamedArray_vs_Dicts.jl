# This script benchmarks the performance of using Dicts versus NamedArrays
# These benchmarks serve to inspire the use of NamedArrays instead of commonly used Dicts
# The scripts was modified from: https://discourse.julialang.org/t/performance-of-namedarrays-vs-dictionaries-of-tuples/1854

using NamedArrays, BenchmarkTools, Random

I = [randstring(4) for i=1:100]
J = [randstring(4) for j=1:100]

d1 = Dict((i,j) => rand() for i in I, j in J)
d2 = Dict((i,j) => rand() for i in I, j in J)

println("\nDict comprehension\n", @benchmark d3 = Dict((i,j) => d1[i,j]*d2[i,j] for i in I, j in J))
println("\nDict assignment\n", @benchmark begin
    d4 = Dict{Tuple{String,String}, Float64}()
    for i in I, j in J
        d4[i,j] = d1[i,j]*d2[i,j]
    end
end)

n1 = NamedArray([rand() for i in I, j in J], (I,J))
n2 = NamedArray([rand() for i in I, j in J], (I,J))

println("\nNamedArray comprehension\n", @benchmark n3 = NamedArray([n1[i,j]*n2[i,j] for i in I, j in J], (I,J)))
println("\nNamedArray assignment\n", @benchmark begin
    n4 = NamedArray(zeros(100,100), (I,J))
    for i in I, j in J
        n4[i,j] = n1[i,j]*n2[i,j]
    end
end)

b1 = [rand() for i=1:100, j=1:100]
b2 = [rand() for i=1:100, j=1:100]

println("\nBasic Array comprehension (no name lookup)\n", @benchmark b3 = [b1[i,j]*b2[i,j] for i=1:100, j=1:100])

I = [randstring(4) for i=1:100]
J = [randstring(4) for j=1:100]
K = [randstring(4) for j=1:100]

println("\nNamedArray lookup\n", @benchmark begin
    n4 = NamedArray(zeros(100,100, 100), (I, J, K))
    for i in I, j in J, k in K
        n4[i,j,k]
    end
end)

d4 = Dict{Tuple{String,String, String}, Float64}()
d1 = Dict((i,j) => rand() for i in I, j in J)
d2 = Dict((i,j) => rand() for i in I, j in J)
for i in I, j in J, k in K
    d4[i,j,k] = d1[i,j]*d2[i,j]*2
end

println("\nDict lookup\n", @benchmark begin
    for i in I, j in J, k in K
        d4[i,j,k]
    end
end)