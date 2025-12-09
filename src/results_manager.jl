using JSON, LibGit2, Dates, CSV, DataFrames

struct SolutionInfo
    instance :: String
    method :: String
    gap :: Float64
    time :: Float64
    status :: String
    cost :: Float64
    solution :: Vector{Vector{Int}}
end

function write_solution_info_to_raw_file(
    sol :: SolutionInfo,
    raw_results_folder :: String = "results/raw"
)
    instance_res_dir = "$raw_results_folder/$(sol.instance)"
    if !isdir(instance_res_dir)
        mkdir(instance_res_dir)
    end

    raw_results_file = "$instance_res_dir/$(sol.method).csv"
    if !isfile(raw_results_file)
        line = "instance,method,gap,time,status,cost,date,commit-id\n"
        write(raw_results_file, line)
    end

    cmd :: Cmd= `git log --pretty=tformat:'%h' -n1 .`
    commit_id = read(cmd, String)
    date = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")

    newline = "$(sol.instance),$(sol.method),$(sol.gap),$(sol.time),$(sol.status),$(sol.cost),$date,$commit_id"

    open(raw_results_file, "a") do f
        write(f, newline)
    end
    return
end

function compile_raw_results(
    raw_results_folder :: String,
    output_res_file :: String,
    criteria :: Vector{Symbol};
    descending :: Bool = false,
    instance_folder :: String = "data/"
)
    
    output = Dict()

    # Build the output by collecting the best result along the given criterion
    for instance in readdir(instance_folder)
        if instance == "solution.txt"
            continue
        end

        instance_raw_results = "$raw_results_folder/$instance"
        if !isdir(instance_raw_results)
            println("No results for instance $instance, skipping")
            continue
        end
        
        # Needed only in a future version where we are not creating the output from scratch
        if instance ∉ keys(output)
            output[instance] = Dict()
        end

        # Process every method
        for method_file in readdir(instance_raw_results)
            method = splitext(method_file)[1]
            filepath = "$instance_raw_results/$method_file"
            method_data = CSV.read(filepath, DataFrame, header = true)
            sort!(method_data, criteria, rev = descending)
            # We want to keep the best method according to these criteria
            best_data = method_data[1, :]
            method_output = Dict()
            method_output["time"] = best_data[:time]
            method_output["gap"] = best_data[:gap]
            method_output["status"] = best_data[:status]
            method_output["cost"] = best_data[:cost]
            method_output["date"] = best_data[:date]
            method_output["commit-id"] = best_data[Symbol("commit-id")]

            output[instance][method] = method_output
        end

    end

    open(output_res_file, "w") do f
        content = JSON.json(output, pretty = true, allownan = true)
        write(f, content)
    end

    return
end

RUN_SCRIPTING = false
if RUN_SCRIPTING
    compile_raw_results(
        "results/raw/",
        "results/static_dual_5sec.json",
        [Symbol("gap")],
        descending = false
    )
end