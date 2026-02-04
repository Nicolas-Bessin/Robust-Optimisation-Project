include("master.jl")
include("pricing.jl")

include("../heuristic/greedy.jl")
include("../data/results_manager.jl")

function CG_solver_non_compact(
    instance :: Data;
    find_integer_sol :: Bool = true,
    timelimit :: Int = 600,
    itermax :: Int = 1000,
    eps = 1e-6,
    verbose = false
)
    total_master_time = 0.0
    total_pricing_time = 0.0

    initial_sol = greedy_init(instance, verbose = true)
    available_partitions = deepcopy(initial_sol)

    master_model = master_problem(instance, initial_sol)
    set_silent(master_model)
    pricing_model = pricing_model_quadratic(instance)
    set_silent(pricing_model)

    # Retrieving references to the variables & constraints
    z = master_model[:z]
    used_clusters = master_model[:used_clusters]
    cover = master_model[:cover]
    rdist = master_model[:rdist]

    lower_bound_obj = 0.
    last_rc = -Inf
    iter = 0
    stop = false
    while !stop && iter < itermax && total_master_time + total_pricing_time < timelimit
        if verbose
            println("Starting iteration $iter, solving master problem")
        end
        # master_model step
        T0master = time()
        optimize!(master_model)
        @assert is_solved_and_feasible(master_model) "master resolution failed at iter $iter"

        curr_obj = objective_value(master_model)

        total_master_time += time() - T0master

        if verbose
            println("Current master objective is $curr_obj, starting pricing resolution")
        end
        # Pricing step
        T0pricing = time()

        theta = dual(used_clusters)
        lambda = dual.(cover)
        mu = dual.(rdist)

        rc, clust = find_best_candidate_quadratic(instance, pricing_model, lambda, mu)

        if verbose
            println("Found a candidate with reduced cost of $(rc - theta)")
        end
        # Note that the pricing is implemented as a min problem (to go with the master_model which is also a min problem)
        # So the reduced cost differs a little from the problem described in the Overleaf
        if rc - theta < -eps
            add_cluster_to_master(instance, master_model, clust)
            push!(available_partitions, clust)
        else 
            stop = true
        end

        lower_bound_obj = curr_obj + instance.K * (rc - theta) 
        last_rc = (rc - theta)
        total_pricing_time += time() - T0pricing

        iter += 1
    end

    if !find_integer_sol
        return
    end

    # Solve with integer constraint using only these columns (No Branch & Price)
    println("------------")
    println("Relaxation solved after $iter iterations, last reduced cost  is $last_rc, solving with integer constraint")
    T0integer = time()

    set_binary.(z)

    optimize!(master_model)

    @assert is_solved_and_feasible(master_model)

    integer_obj = objective_value(master_model)

    gap_estimation = (integer_obj - lower_bound_obj) / lower_bound_obj
    status = gap_estimation < eps ? string(MOI.OPTIMAL) : "CG_REL_ENDED"
    z_val = value.(z)

    solution = [clust for (p, clust) in enumerate(available_partitions) if z_val[p] > 1 - eps]

    total_integer_time = time() - T0integer


    feas = check_feasability(instance, solution, robust = true)
    println("Feasibility : $feas")

    println("Integer solution has a value of $integer_obj, estimated gap (to relaxation) is $gap_estimation")
    println("------ Time Info ------")
    println("Master : $total_master_time, Pricing : $total_pricing_time, Integer : $total_integer_time")
    println("Total solving time : $(total_integer_time + total_pricing_time + total_master_time)")

    solution_info = SolutionInfo(
        instance.instance_name,
        "CG-NC",
        gap_estimation,
        total_master_time + total_pricing_time,
        status,
        integer_obj,
        solution
    )

    write_solution_info_to_raw_file(solution_info)

    return status
end

# data = parse_file("data/22_ulysses_6.tsp");

# @time sol_info = CG_solver_non_compact(data)

# check_feasability(data, sol_info.solution)
# compute_robust_cost(data, sol_info.solution)