include("instance.jl")
include("results_manager.jl")
include("cutting_planes_separation.jl")
include("cutting_planes.jl")
include("utils.jl")

using JuMP, Gurobi



function cutting_planes_with_callbacks(data :: Data, timelimit :: Int = 600, eps :: Float64 = 1e-6)
    N = data.N
    K = data.K

    gurobi_env = Gurobi.Env()
    METHOD = "cutting_planes_callbacks"

    t0 = time()
    
    model = cutting_planes_initial_model(data, gurobi_env)
    set_time_limit_sec(model, timelimit)

    x = model[:x]
    y = model[:y]
    z = model[:z]

    # Prepare the separation problems 
    model_SPL = cutting_planes_length_separation_model(data, gurobi_env)
    model_SPW = cutting_planes_weights_separation_model(data, gurobi_env)

    cutting_planes_length = 0
    cutting_planes_weight = 0
    total_time_separation = 0.

    function lazy_cuts(cb_data)
        t0sep = time()
        status = callback_node_status(cb_data, model)
        if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
            return
        end
        x_val = callback_value.(cb_data, x)
        y_val = callback_value.(cb_data, y)
        z_val = callback_value(cb_data, z)

        # Length-related cutting planes 
        val, delta1 = separation_length(model_SPL, x_val, data)
        if val > z_val + eps 
            # Build the length cut
            # println("Adding a length cut")
            lcut = @build_constraint( 
                sum(x[i, j] * delta1[i, j] * (data.l_hat[i] + data.l_hat[j]) for i in 1:N, j in i+1:N) <= z
            )
            MOI.submit(model, MOI.LazyConstraint(cb_data), lcut)
            cutting_planes_length += 1
        end

        # Weight related cutting planes
        for k in 1:K
            val, delta2 = separation_weights(model_SPW, y_val, k, data)
            if val > data.B
                # println("Adding a weight cut")             
                for kprime in 1:K
                    wcut = @build_constraint(
                        sum( (1 + delta2[i]) * y[i, kprime] * data.weights[i] for i in 1:N) <= data.B
                    )
                    MOI.submit(model, MOI.LazyConstraint(cb_data), wcut)
                end
                cutting_planes_weight += K
            end
        end

        total_time_separation += time() - t0sep
        return
    end

    set_attribute(model, MOI.LazyConstraintCallback(), lazy_cuts)

    optimize!(model)

    status = string(termination_status(model))
    solving_time = time() - t0 # Innacurate on first run due to compilation :(
    println("Status is $status")
    # Time info
    separation_percentage = trunc(100 * total_time_separation / solving_time, digits = 3)
    println("Time spent in the callback function : $(trunc(total_time_separation, digits = 3)) ($separation_percentage %)")
    # Cutting planes info
    println("Added a total of $cutting_planes_length length related cutting planes")
    println("Added a total of $cutting_planes_weight weight related cutting planes")

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