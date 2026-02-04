include("cutting_planes/cutting_planes.jl")
include("cutting_planes/cutting_planes_callbacks.jl")
include("compact/robust_dual.jl")
include("compact/static.jl")
include("non-compact/column_generation.jl")

"""
Methods : 1 for static, 2 for cutting planes, 3 for robust via dualization, ...
"""
function run_all_instances(
    instances_dir :: String,
    timelimit :: Int,
    methods :: Vector{Int} = [1, 3, 4];
    stop_at_first_failure :: Bool = true
)
    instances_name = readdir(instances_dir)
    # The solution.txt file is not an instance
    deleteat!(instances_name, findfirst(x -> x == "solution.txt", instances_name))
    all_instances = ["$instances_dir/$name" for name in instances_name]
    run_list_of_instances(
        all_instances,
        timelimit,
        methods,
        stop_at_first_failure = stop_at_first_failure
    )
end

"""
Methods : 1 for static, 2 for cutting planes, 3 for robust via dualization, ...
"""
function run_list_of_instances(
    instances_list :: Vector{String},
    timelimit :: Int64,
    methods :: Vector{Int} = [1, 3, 4];
    stop_at_first_failure :: Bool = true
)
    if stop_at_first_failure
        @assert length(methods) == 1 "When stopping at a failure, only compute for 1 method at a time"
    end

    OPT = string(MOI.OPTIMAL)
    CG_ENDED = "CG_REL_ENDED"
    # Sort the instances by increasing size Nvertices (filenames are 'Nvertices_name_Kclust.tsp')
    filesizes = [parse(Int, split(basename(x), "_")[1]) for x in instances_list]
    perm = sortperm(filesizes)
    println(instances_list[perm])
    println("------------")

    for filepath in instances_list[perm]
        println("--------------")
        println("Running instance $filepath")
        instance = parse_file(filepath)

        if 1 in methods
            println("----- Static -----")
            status = static_problem(instance, timelimit)
        end
        
        if 2 in methods
            println("----- Cutting Planes -----")
            status = cutting_planes_method(instance, timelimit)

        end

        if 3 in methods 
            println("----- Dualization -----")
            status = robustdual_problem(instance, timelimit)
        end
        
        if 4 in methods 
            println("----- Cutting Planes - Lazy Callback -----")
            status = cutting_planes_with_callbacks(instance, timelimit)
        end

        if 5 in methods
            println("----- Non Compact by CG -----")
            status = CG_solver_non_compact(instance, timelimit = timelimit)
        end

        if status ∉ [OPT, CG_ENDED] && stop_at_first_failure
            print("Status is $status, expected $OPT")
            break
        end

    end
    
end

run_all_instances("data/", 600, [1], stop_at_first_failure = false)