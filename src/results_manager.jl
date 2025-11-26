using JSON, LibGit2, Dates

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
        line = "instance, method, gap, time, status, cost, date, commit-id \n"
        write(raw_results_file, line)
    end

    cmd :: Cmd= `git log --pretty=tformat:'%h' -n1 .`
    commit_id = read(cmd, String)
    date = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")

    newline = "$(sol.instance), $(sol.method), $(sol.gap), $(sol.time), $(sol.status), $(sol.cost), $date, $commit_id"

    open(raw_results_file, "a") do f
        write(f, newline)
    end
    return
end