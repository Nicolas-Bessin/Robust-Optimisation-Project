include("../data/instance.jl")
include("../data/results_manager.jl")
include("../data/utils.jl")

using JuMP, Gurobi

function static_problem(instance :: Data, timelimit :: Int = 600)
    METHOD = "static"

    model = Model(Gurobi.Optimizer)
    set_time_limit_sec(model, timelimit)
    set_silent(model)

    N = instance.N
    K = instance.K

    @variable(model, x[i = 1:N, j = 1:N], Bin)
    @variable(model, y[i = 1:N, k = 1:K], Bin)

    @constraint(model, [i = 1:N, j = 1:N, k = 1:K],
        x[i,j] >= y[i,k] + y[j, k] - 1
    )
    @constraint(model, [i = 1:N],
        sum(y[i,k] for k in 1:K) == 1
    )
    @constraint(model, [k in 1:K],
        sum(y[i,k] * instance.weights[i] for i in 1:N) <= instance.B
    )

    @objective(model, Min,
        sum(x[i,j] * instance.edge_lengths[i,j] for i in 1:N, j in i+1:N)
    )

    optimize!(model)

    status = string(termination_status(model))
    solving_time = solve_time(model)

    if !is_solved_and_feasible(model)
        sol = SolutionInfo(
            instance.instance_name,
            METHOD,
            Inf,
            solving_time,
            status,
            Inf,
            [[]]
        )
        write_solution_info_to_raw_file(sol)
        return status
    end

    gap = relative_gap(model)
    cost = objective_value(model)

    x_val = value.(x)
    y_val = value.(y)
    partitions = rebuild_partition(y_val, instance)

    println("Partition is $partitions")

    sol = SolutionInfo(
        instance.instance_name,
        METHOD,
        gap,
        solving_time,
        status,
        cost,
        partitions
    )

    write_solution_info_to_raw_file(sol)

    return status
end