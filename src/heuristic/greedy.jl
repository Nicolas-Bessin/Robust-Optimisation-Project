include("../data/instance.jl")
include("../data/utils.jl")

using Clustering

function dist(
    x1 :: Tuple{Float64, Float64},
    x2 :: Tuple{Float64, Float64}
    )
    return sqrt(sum((x1 .- x2).^2))
end

"""
 - Perform a K-means clustering to initialize the clusters
 - Distribute the nodes by decreasing weight to their nearest cluster with available budget
 - We know how to ensure the robust budget constraint : distribute the uncertainty in priority to the highest static weight nodes
"""
function greedy_init(
    instance :: Data;
    verbose :: Bool = false
) :: Vector{Vector{Int}}
    N = instance.N
    K = instance.K

    # coords_flat = reshape(reinterpret(Float64, instance.coordinates), (2, :))

    # clusters = Clustering.kmeans(coords_flat, K)
    # centers = [Tuple(col) for col in eachcol(clusters.centers)]

    # Initialize the cluster
    partition = [Int[] for k in 1:K]
    used_budget_by_cluster = zeros(K)
    remaining_uncertainty_by_cluster = instance.W * ones(K)

    # Sort the nodes by weights
    perm = sortperm(1:N, by = i -> instance.weights[i], rev = true)

    for node in perm

        node_weight = instance.weights[node]
        # centers_dist = [dist(instance.coordinates[node], centers[i]) for i in 1:K]
        # clusters_by_distance = sortperm(centers_dist)
        parts_by_uncertainty = sortperm(remaining_uncertainty_by_cluster)

        # Find the first part that has no uncertainty left and has budget left
        # Else, find the first part that has budget left
        assigned = false
        index = 1
        while index <= K && !assigned
            # Try the cluster
            part = parts_by_uncertainty[index]
            used_budg = used_budget_by_cluster[part]
            rem_uncert = remaining_uncertainty_by_cluster[part]
            available_delta = min(instance.delta_2_max[node], rem_uncert)

            if node_weight * (1 + available_delta) + used_budg < instance.B
                # Assign the node to this cluster
                push!(partition[part], node)
                used_budget_by_cluster[part] += node_weight * (1 + available_delta)
                remaining_uncertainty_by_cluster[part] -= available_delta
                assigned = true
            else
                index += 1
            end
        end

        if index > K
            error("Assignment failed, either algorithm is false or no assignment is possible")
        end

    end

    if verbose
        # Computing the cost
        cost = 0.0
        for clust in partition
            cost += cluster_static_cost(clust, instance)
        end

        println("Partition is :")
        println(partition)
        println("Budget is $(instance.B), worst case scenario budget usage by partition is :")
        println(used_budget_by_cluster)
        println("Total cost is : $cost")
    end
    
    return partition
end