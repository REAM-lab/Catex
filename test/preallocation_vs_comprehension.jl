# In this benchmark test, we test preallocation vs comprehension for NamedArrays
# As of now Tue, Dec 2, 2025, Version 1.12.2 (2025-11-20), 1 wins in speed and allocation among 1, 2, 3
# Still 1 cannot beat 4 and 5, which are pure matrices. Among 4 and 5, 4 wins
using BenchmarkTools, NamedArrays

ids = collect(1:10^6)

# 1
@benchmark S = NamedArray([i+1 for i in ids], (ids), :i)

# 2
@benchmark begin 
    S = NamedArray(Vector{Int64}(undef, length(ids)), (ids), :i)
    for i in ids
        S[i] = i + 1
    end
end

# 3
@benchmark begin 
    S = Vector{Int64}(undef, length(ids))
    for i in ids
        S[i] = i + 1
    end
    S = NamedArray(S, (ids), :i)
end

# Here we 
# 4
@benchmark T = [i+1 for i in ids]

# 5 
@benchmark begin 
    T = Vector{Int64}(undef, length(ids))
    for i in ids
        T[i] = i + 1
    end
end