include("../data/instance.jl")
include("../data/results_manager.jl")
include("../data/utils.jl")

include("../heuristic/greedy.jl")

using JuMP, Gurobi

function robustdual_problem(
    instance :: Data,
    timelimit :: Int = 600;
    initial_sol = nothing
    )

    METHOD = "robust_dual"

    model = Model(Gurobi.Optimizer)
    set_time_limit_sec(model, timelimit)
    set_silent(model)

    N = instance.N
    K = instance.K

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
        z + eta[i,j] >= x[i,j]*(instance.l_hat[i] + instance.l_hat[j])
    )
    @constraint(model, [k = 1:K],
        instance.W*t[k] + sum(y[i,k] * instance.weights[i] + xi[i,k]*instance.delta_2_max[i] for i in 1:N) <= instance.B
    )
    @constraint(model, [i = 1:N, k = 1:K],
        t[k] + xi[i,k] >= y[i,k] * instance.weights[i]
    )

    # Breaking symetries

    # First method : imposing the biggest partition to be first 
    # @constraint(model, [k = 1:(K-1)],
    #     sum(y[i, k] for i ∈ 1:N) >= sum(y[i, k+1] for i ∈ 1:N)
    # )
    # Res (22_ulysses_3.tsp): Absolutely horrible : 4s to solve without, 24s to solve with
   
    # Second method : Imposing the first vertex to be in partition 1
    # @constraint(model, y[N, 1] == 1)
    # Res (22_ulysses_3.tsp) : Seems to make it slightly faster, from ~4s to ~3.6s
    # However : if we impose the last vertex to be in the first partition : ~4s -> ~6s
    

    @objective(model, Min,
        instance.L * z + sum(x[i,j] * instance.edge_lengths[i,j] + 3*eta[i,j] for i in 1:N, j in i+1:N)
    )

    # Adding the initial solution
    if !isnothing(initial_sol)
        println("Setting the initial solution of the problem")

        @assert typeof(initial_sol) == Vector{Vector{Int}}
        @assert length(initial_sol) <= instance.K "Solution has to many partitions for the current instance"

        for (k, cluster) in enumerate(initial_sol)
            # Starting x values
            for i in cluster
                for j in cluster 
                    if j > i 
                        set_start_value(x[i, j], 1)
                    end
                end
            end
            # Starting y values
            for i in cluster
                set_start_value(y[i, k], 1)
            end
        end
    end

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
    println("With a cost of $cost")

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


# instance = parse_file("instance/22_ulysses_6.tsp");
# @time robustdual_problem(instance)
# @time robustdual_problem(instance, initial_sol = greedy_init(instance))