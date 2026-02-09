using JSON

json_file = "results/all_10_minutes_new.json"
tex_file  = "results/results_list.txt"

data = JSON.parsefile(json_file, allownan = true)

# Escape LaTeX special characters
function tex_escape(s::AbstractString)
    replace(s,
        "_" => "\\_",
        "%" => "\\%",
        "&" => "\\&",
        "#" => "\\#",
        "{" => "\\{",
        "}" => "\\}"
    )
end

open(tex_file, "w") do io
    println(io, "% Auto-generated")
    println(io, "\\begin{itemize}")

    data = JSON.parsefile(json_file, allownan = true)
    instances = collect(keys(data))
    sort!(instances, by = x -> (parse(Int, split(x, "_")[1]), x)) # Sort by instance number

    for instance in instances
        methods = data[instance]
        for (method, res) in methods
            cost     = res["cost"]
            time     = res["time"]
            gap      = res["gap"]
            solution = res["solution"]
            solution_str = replace(solution, ";" => ",")[2:end-1]

            println(io, "  \\item \\textbf{Instance:} $(tex_escape(instance))")
            println(io, "    \\begin{itemize}")
            println(io, "      \\item Method: $(tex_escape(method))")
            println(io, "      \\item Cost: $(cost)")
            println(io, "      \\item Time: $(time)")
            println(io, "      \\item Gap: $(gap)")
            println(io, "      \\item Solution: $(solution_str)")
            println(io, "    \\end{itemize}")
        end
    end

    println(io, "\\end{itemize}")
end

println("LaTeX list written to $tex_file")
