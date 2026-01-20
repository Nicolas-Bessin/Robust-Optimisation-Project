include("../data/instance.jl")
include("../data/results_manager.jl")
include("../data/utils.jl")

include("cutting_planes_separation.jl")

using JuMP, Gurobi


""" 
Creates the initial master problem of the cutting planes formulation
"""
function cutting_planes_initial_model(data :: Data, gurobi_env)
    N = data.N
    K = data.K

    model = Model(() -> Gurobi.Optimizer(gurobi_env))
    set_attribute(model, "Method", 1)
    set_silent(model)

    @variable(model, x[i = 1:N, j = i+1:N], Bin)
    @variable(model, y[i = 1:N, k = 1:K], Bin)
    @variable(model, z >= 0)

    @constraint(model, [i = 1:N, j = i+1:N, k = 1:K],
        x[i,j] >= y[i,k] + y[j, k] - 1
    )
    @constraint(model, [i = 1:N],
        sum(y[i,k] for k in 1:K) == 1
    )
    @constraint(model, [k in 1:K],
        sum(y[i,k] * data.weights[i] for i in 1:N) <= data.B
    )

    @objective(model, Min,
        sum(x[i,j] * data.edge_lengths[i,j] for i in 1:N, j in i+1:N) + z
    )

    return model
end

function cutting_planes_method(data :: Data, timelimit :: Int = 600, eps :: Float64 = 1e-6)
    N = data.N
    K = data.K

    gurobi_env = Gurobi.Env()
    METHOD = "cutting_planes"

    t0 = time()
    
    model = cutting_planes_initial_model(data, gurobi_env)
    x = model[:x]
    y = model[:y]
    z = model[:z]

    # Prepare the separation problems 
    model_SPL = cutting_planes_length_separation_model(data, gurobi_env)
    model_SPW = cutting_planes_weights_separation_model(data, gurobi_env)

    optimality_reached = false
    iter = 0
    cutting_planes_length = 0
    cutting_planes_weight = 0

    total_time_master = 0.
    total_time_separation = 0.

    while !optimality_reached && time() - t0 < timelimit

        master_time = @elapsed begin
        # Re-optimize model
        optimize!(model)
        @assert is_solved_and_feasible(model)
        optimality_reached = true # we will set it back to false as needed
        x_val = value.(x)
        y_val = value.(y)
        end
        total_time_master += master_time

        separation_time = @elapsed begin
        # Solve both separation problems, update optimality and add the cuts
        # 1) Length separation
        remaining = floor(Int, timelimit - (time() - t0)) + 1
        if remaining < 1
            optimality_reached = false
            break
        end
        val, delta1 = separation_length(model_SPL, x_val, data, remaining)
        if val > value.(z) + eps 
            optimality_reached = false
            # Add the length cut
            # println("Adding a length cut")
            @constraint(model, 
                sum(x[i, j] * delta1[i, j] * (data.l_hat[i] + data.l_hat[j]) for i in 1:N, j in i+1:N) <= z
            )
            cutting_planes_length += 1
        end
        # 2) Weight separation
        # Parallelisation ?
        for k in 1:K
            remaining = floor(Int, timelimit - (time() - t0)) + 1
            if remaining < 1
                break
            end
            val, delta2 = separation_weights(model_SPW, y_val, k, data, remaining)
            if val > data.B
                optimality_reached = false
                # Add the weight cut for all values of k !
                # println("Adding a weight cut")
                @constraint(model, [k = 1:K],
                    sum( (1 + delta2[i]) * y[i, k] * data.weights[i] for i in 1:N) <= data.B
                )
                cutting_planes_weight += K
            end
        end

        end
        total_time_separation += separation_time 
        iter += 1
    end

    status = string(termination_status(model))
    # We can't simply measure resolution time of the master problem
    # solving_time = solve_time(model) + solve_time(model_SPL) 
    #    + sum(solve_time(model_SPW_k) for model_SPW_k in models_SPW) # Innacurate due to all the non-solving stuff
    solving_time = time() - t0 # Innacurate on first run due to compilation :(

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
        return status
    end

    # WARNING : does not mean anything if cutting plane method didn't finish, because primal solution is thus not admissible
    gap = relative_gap(model) 
    cost = objective_value(model)

    x_val = value.(x)
    y_val = value.(y)
    partitions = rebuild_partition(y_val, data)
    
    println("Partition is $partitions")
    println("Objective value is $(objective_value(model))")
    println("Added a total of $cutting_planes_length length related cutting planes")
    println("Added a total of $cutting_planes_weight weight related cutting planes")

    # Time info
    master_percentage = trunc(100 * total_time_master / solving_time, digits = 3)
    separation_percentage = trunc(100 * total_time_separation / solving_time, digits = 3)
    println("Time spent solving master problem : $(trunc(total_time_master, digits = 3)) ($master_percentage %)")
    println("Time spent solving separation problem : $(trunc(total_time_separation, digits = 3)) ($separation_percentage %)")

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

    return status
end