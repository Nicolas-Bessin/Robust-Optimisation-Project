include("instance.jl")

"""
Translate the resulting values of y to human-readable partition
"""
function rebuild_partition(yval, instance :: Data)
    N = instance.N
    K = instance.K
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

"""
Compute the static cost of a cluster
"""
function cluster_static_cost(
    cluster :: Vector{Int},
    instance :: Data
)
    cost = 0.0

    for i in cluster
        for j in cluster
            if i >j
                cost += instance.edge_lengths[i, j]
            end
        end
    end

    return cost
end