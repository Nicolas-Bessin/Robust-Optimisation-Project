include("../instance/instance.jl")

using JuMP, Gurobi

function pricing_model_quadratic(
    instance :: Data
)
    N = instance.N
    
    model = Model(Gurobi.Optimizer)
    set_silent(model)

    @variable(model, x[i = 1:N], Bin)
    @variable(model, ξ[i = 1:N] >= 0)
    @variable(model, t >= 0)

    @constraint(model, [i  = 1:N],
        instance.weights[i] * x[i] <= t + ξ[i]
    )
    @constraint(model, 
        sum(instance.weights[i] * x[i] for i ∈ 1:N) <= instance.B - instance.W * t - sum(instance.delta_2_max[i] * ξ[i] for i ∈ 1:N)
    )

    return model
end

function find_best_candidate_quadratic(
    instance :: Data,
    pricing_model :: JuMP.GenericModel{Float64},
    lambda,
    mu;
    eps = 1e-6
)
    N = instance.N

    x = pricing_model[:x]

    @objective(pricing_model, Min,
        - sum(x[i] * lambda[i] for i ∈ 1:N) 
            + sum(
                x[i] * x[j] * (instance.edge_lengths[i, j] + mu[i, j] * (instance.l_hat[i] + instance.l_hat[j])) for i in 1:N, j in i+1:N
                )
    )

    optimize!(pricing_model)

    @assert is_solved_and_feasible(pricing_model)

    # Build the cluster
    x_val = value.(x)
    cluster = [i for i ∈ 1:N if x_val[i] > 1 - eps]

    reduced_cost = objective_value(pricing_model)

    return reduced_cost, cluster
end