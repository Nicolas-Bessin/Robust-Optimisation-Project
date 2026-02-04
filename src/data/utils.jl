include("instance.jl")

using JuMP, Gurobi

"""
Translate the resulting values of y to human-readable partition
"""
function rebuild_partition(yval, instance :: Data) :: Vector{Vector{Int}}
    N = instance.N
    K = instance.K
    # Re-build the partitions in a printable way
    partitions :: Vector{Vector{Int}} = repeat([[]], K)
    # println(length(partitions))
    # println(yval)
    for i in 1:N
        k = findfirst(val -> abs(val - 1) < 1e-6, yval[i, :])
        if isnothing(k)
            println("Error : vertex $i is not in a partition : we got k = $k")
        end
        # println("Vertex $i is in partition $k")
        push!(partitions[k], i)
    end

    partitions = filter(p -> !isempty(p), partitions)

    return partitions
end

"""
Compute the static cost of a clusterer
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


"""
Checks the feasibility of a solution
"""
function check_feasability(
    instance :: Data,
    solution :: Vector{Vector{Int}};
    robust :: Bool = true
)

    # Use the fact that we know what is the worst assignment of delta_2 for a given cluster
    for cluster in solution
        ncluster = length(cluster)
        perm = sortperm(1:ncluster, by = x -> instance.weights[cluster[x]], rev = true)

        cluster_cost = sum(instance.weights[cluster])
        available_uncert = instance.W
        if !robust
            available_uncert = 0
        end
        idx = 1
        while idx <= ncluster && available_uncert > 0
            node = cluster[perm[idx]]
            available_delta = min(available_uncert, instance.delta_2_max[node])
            cluster_cost += instance.weights[node] * available_delta
            available_uncert -= available_delta

            idx += 1
        end

        if cluster_cost > instance.B
            println("INFEASIBILITY : Cluster $cluster ")
            println("Worst case cost of $cluster_cost ; Budget is $(instance.B)")
            return false
        end
    end
        
    return true
end


"""
Compute the robust cost of this clustering assignment
"""
function compute_robust_cost(
    instance :: Data,
    solution :: Vector{Vector{Int}}
)
    N = instance.N

    # Build the list of interior edges
    interior_edges :: Vector{Tuple{Int, Int}} = []

    for clust in solution
        # We must have generic builder since no guarantee is expected on the sortedness of clusters
        for v1 in clust
            for v2 in clust
                if v2 > v1
                    push!(interior_edges, (v1, v2))
                end
            end
        end
    end

    model = Model(Gurobi.Optimizer)
    set_silent(model)
    
    @variable(model, 0 <= delta1[i = 1:N, j = i+1:N] <= instance.delta_1_max)

    @constraint(model, sum(delta1) <= instance.L)

    @objective(model, Max,
        sum(delta1[v1, v2] * (instance.l_hat[v1] + instance.l_hat[v2]) for (v1, v2) in interior_edges)
    )

    optimize!(model)
    @assert is_solved_and_feasible(model) "Robust computation was not feasible, something went very wrong"

    static_cost = sum(
        cluster_static_cost(clust, data) for clust in solution
    )

    return static_cost + objective_value(model)
end