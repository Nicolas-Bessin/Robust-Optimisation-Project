include("../data/instance.jl")

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
    instance :: Data
)
    N = instance.N
    K = instance.K

    coords_flat = reshape(reinterpret(Float64, instance.coordinates), (2, :))

    clusters = Clustering.kmeans(coords_flat, K)
    centers = [Tuple(col) for col in eachcol(clusters.centers)]

    # Initialize the cluster
    partition = [[] for k in 1:K]
    used_budget_by_cluster = zeros(K)
    remaining_uncertainy_by_cluster = instance.W * ones(K)

    # Sort the nodes by weights
    perm = sortperm(1:N, by = i -> instance.weights[i], rev = true)
    println(perm)

    for node in perm

        node_weight = instance.weights[node]
        centers_dist = [dist(instance.coordinates[node], centers[i]) for i in 1:K]
        clusters_by_distance = sortperm(centers_dist)

        # Find the closest cluster we can assign this node to
        assigned = false
        index = 1
        while index <= K && !assigned
            # Try the cluster
            cluster = clusters_by_distance[index]
            used_budg = used_budget_by_cluster[cluster]
            rem_uncert = remaining_uncertainy_by_cluster[cluster]
            available_delta = min(instance.delta_2_max[node], rem_uncert)

            if node_weight * (1 + available_delta) + used_budg < instance.B
                # Assign the node to this cluster
                push!(partition[cluster], node)
                used_budget_by_cluster[cluster] += node_weight * (1 + available_delta)
                remaining_uncertainy_by_cluster[cluster] -= available_delta
                assigned = true
            else
                index += 1
            end
        end

        if index > K
            error("Assignment failed, either algorithm is false or no assignment is possible")
        end

    end

    println("Partition is :")
    println(partition)
    println("Budget is $(instance.B), worst case scenario budget usage by partition is :")
    println(used_budget_by_cluster)
    
end


data = parse_file("data/22_ulysses_6.tsp");
greedy_init(data)