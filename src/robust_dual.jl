include("instance.jl")
include("results_manager.jl")
include("utils.jl")

using JuMP, Gurobi

function robustdual_problem(data :: Data, timelimit :: Int = 600)
    METHOD = "robust_dual"

    model = Model(Gurobi.Optimizer)
    set_time_limit_sec(model, timelimit)
    set_silent(model)

    N = data.N
    K = data.K

    @variable(model, x[i = 1:N, j = i+1:N], Bin)
    @variable(model, y[i = 1:N, k = 1:K], Bin)
    @variable(model, z >= 0)
    @variable(model, eta[i = 1:N, j = i+1:N] >= 0)
    @variable(model, t[k = 1:K] >= 0)
    @variable(model, xi[i = 1:N, k = 1:K] >= 0)

    @constraint(model, [i = 1:N, j = i+1:N, k = 1:K],
        x[i,j] >= y[i,k] + y[j, k] - 1
    )
    @constraint(model, [i = 1:N],
        sum(y[i,k] for k in 1:K) == 1
    )
    @constraint(model, [i = 1:N, j = i+1:N],
        z + eta[i,j] >= x[i,j]*(data.l_hat[i] + data.l_hat[j])
    )
    @constraint(model, [k = 1:K],
        data.W*t[k] + sum(y[i,k] * data.weights[i] + xi[i,k]*data.delta_2_max[i] for i in 1:N) <= data.B
    )
    @constraint(model, [i = 1:N, k = 1:K],
        t[k] + xi[i,k] >= y[i,k] * data.weights[i]
    )

    @objective(model, Min,
        data.L * z + sum(x[i,j] * data.edge_lengths[i,j] + 3*eta[i,j] for i in 1:N, j in i+1:N)
    )

    optimize!(model)

    status = string(termination_status(model))
    solving_time = solve_time(model)

    if !is_solved_and_feasible(model)
        sol = SolutionInfo(
            data.instance_name,
            METHOD,
            Inf,
            solving_time,
            status,
            Inf,
            [[]]
        )
        write_solution_info_to_raw_file(sol)
        return
    end

    gap = relative_gap(model)
    cost = objective_value(model)

    x_val = value.(x)
    y_val = value.(y)
    partitions = rebuild_partition(y_val, data)

    println("Partition is $partitions")
    println("With a cost of $cost")

    sol = SolutionInfo(
        data.instance_name,
        METHOD,
        gap,
        solving_time,
        status,
        cost,
        partitions
    )

    write_solution_info_to_raw_file(sol)

end

data = parse_file("data/22_ulysses_3.tsp");

@time robustdual_problem(data);