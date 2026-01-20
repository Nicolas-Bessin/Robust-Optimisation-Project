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

    available_partitions = deepcopy(initial_sol)
    available_partitions_static_costs = [
        cluster_static_cost(clust, data) for clust in available_partitions
    ]
    n_parts =  length(available_partitions)

    println("Initializing the master problem with $n_parts clusters")

    model = Model(Gurobi.Optimizer)
    # set_silent(model)

    @variable(model, 0 <= z[1:n_parts] <= 1)
    @variable(model, alpha >= 0)
    @variable(model, beta[i = 1:N, j = i+1:N] >= 0)

    # Number of clusters we can use
    @constraint(model, sum(z) <= K)
    # Covering all vertices
    @constraint(model, [i = 1:N], 
        sum(z[p] for p in 1:n_parts if i ∈ available_partitions[p]) >= 1
        )
    # Robust distances
    @constraint(model, [i = 1:N, j = i+1:N],
        alpha + beta[i, j] >= (instance.l_hat[i] + instance.l_hat[j]) * 
            sum(z[p] for p in 1:n_parts if (i ∈ available_partitions[p] && j ∈ available_partitions[p]) )
    )

    @objective(model, Min, 
        instance.L * alpha + instance.delta_1_max * sum(beta) + sum(z[p] * available_partitions_static_costs[p] for p in 1:n_parts)
    )

    return model
end


# include("../heuristic/greedy.jl")
# data = parse_file("data/22_ulysses_6.tsp");

# master_problem(
#     data,
#     greedy_init(data)
# )