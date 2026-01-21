include("../data/instance.jl")
include("../data/utils.jl")

using JuMP, Gurobi


"""
Creates the master problem for the non compact formulation
"""
function master_problem(
    instance :: Data,
    initial_sol :: Vector{Vector{Int}}
)
    N = instance.N
    K = instance.K

    static_costs = [
        cluster_static_cost(clust, instance) for clust in initial_sol
    ]
    n_parts =  length(initial_sol)

    println("Initializing the master problem with $n_parts clusters")

    model = Model(Gurobi.Optimizer)
    # set_silent(model)

    @variable(model, 0 <= z[1:n_parts] <= 1)
    @variable(model, alpha >= 0)
    @variable(model, beta[i = 1:N, j = i+1:N] >= 0)

    # Number of clusters we can use
    @constraint(model, used_clusters, sum(z) <= K)
    # Covering all vertices
    @constraint(model, cover[i = 1:N], 
        sum(z[p] for p in 1:n_parts if i ∈ initial_sol[p]) >= 1
        )
    # Robust distances
    @constraint(model, rdist[i = 1:N, j = i+1:N],
        alpha + beta[i, j] - (instance.l_hat[i] + instance.l_hat[j]) * 
            sum(z[p] for p in 1:n_parts if (i ∈ initial_sol[p] && j ∈ initial_sol[p]) )
            >= 0
    )

    @objective(model, Min, 
        instance.L * alpha + instance.delta_1_max * sum(beta) + sum(z[p] * static_costs[p] for p in 1:n_parts)
    )

    return model
end

"""
Modifies the given master problem to add a nex cluster to the list of available clusters
"""
function add_cluster_to_master(
    instance :: Data,
    master :: JuMP.GenericModel{Float64},
    cluster :: Vector{Int}
)
    N = instance.N

    # Create the new variable
    z = master[:z]
    var_index = length(z) + 1
    push!(z, @variable(master, lower_bound = 0, upper_bound = 1, base_name = "z"))

    # Set the objective & constraint coefficients of this variable
    used_clusters = master[:used_clusters]
    cover = master[:cover]
    rdist = master[:rdist]

    # Constraints
    set_normalized_coefficient(used_clusters, z[var_index], 1)
    for i in cluster
        set_normalized_coefficient(cover[i], z[var_index], 1)
        for j in cluster
            if j > i
                set_normalized_coefficient(rdist[i, j], z[var_index],
                    -instance.l_hat[i] - instance.l_hat[j]
                )
            end
        end
    end

    # Objective
    set_objective_coefficient(master, z[var_index], cluster_static_cost(cluster, instance))

    return
end