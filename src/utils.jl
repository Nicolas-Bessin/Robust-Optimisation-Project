include("instance.jl")

"""
Translate the resulting values of y to human-readable partition
"""
function rebuild_partition(yval, data :: Data)
    N = data.N
    K = data.K
    # Re-build the partitions in a printable way
    partitions :: Vector{Vector{Int}} = repeat([[]], K)
    # println(length(partitions))
    for i in 1:N
        k = findfirst(y -> y == 1, yval[i, :])
        # println("Vertex $i is in partition $k")
        push!(partitions[k], i)
    end

    partitions = filter(p -> !isempty(p), partitions)

    return partitions
end