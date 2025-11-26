include("cutting_planes.jl")
include("robust_dual.jl")
include("static.jl")


"""
Methods : 1 for static, 2 for cutting planes, 3 for robust via dualization, ...
"""
function run_all_instances(
    instance_folder :: String,
    timelimit :: Int64,
    methods :: Vector{Int} = [1, 2, 3]
)

    for instance in readdir(instance_folder)
        println("--------------")
        println("Running instance $instance")
        filepath = "$instance_folder/$instance"
        data = parse_file(filepath)

        if 1 in methods
            println("----- Static -----")
            static_problem(data, timelimit)
        end
        
        if 2 in methods
            println("----- Cutting Planes -----")
            cutting_planes_method(data, timelimit)
        end

        if 3 in methods 
            println("----- Dualization -----")
            robustdual_problem(data, timelimit)
        end

    end
    
end

run_all_instances("data/", 5, [1, 3])