using JSON, Plots, Printf

function import_data(filepath :: String)

    return JSON.parsefile(filepath, allownan = true)
end

function solve_time_graph(data; only_opt_CG = false, tmax = 600)
    instances = keys(data)

    solve_times_per_method = Dict{String, Vector{Float64}}()
    for instance in instances
        results = data[instance]
        methods = keys(results)
        for method in methods
            if !haskey(solve_times_per_method, method)
                solve_times_per_method[method] = Float64[]
            end
            status = results[method]["status"]
            if status == "TIME_LIMIT"
                println("Warning: Time limit reached for instance $instance, method $method. Skipping.")
                continue
            end
            time = results[method]["time"]
            if isnan(time)
                println("Warning: NaN solve time for instance $instance, method $method. Skipping.")
                continue
            end
            if time > tmax
                println("Warning: Solve time $time exceeds maximum threshold $tmax for instance $instance, method $method. Skipping.")
                continue
            end

            # Treat the CG case
            if only_opt_CG && method == "CG-NC"
                cg_is_optimal = false
                if haskey(results, "robust_dual")
                    opt_value = results["robust_dual"]["cost"]
                    if abs(results[method]["cost"] - opt_value) < 1e-3
                        cg_is_optimal = true
                    end
                end
                # if results[method]["gap"] < 1e-3
                #     cg_is_optimal = true
                # end
                if !cg_is_optimal
                    println("Warning: CG-NC is not optimal for instance $instance. Skipping.")
                    continue
                end
            end
            push!(solve_times_per_method[method], results[method]["time"])
        end
    end

    # Sort them by solve times
    for method in keys(solve_times_per_method)
        sort!(solve_times_per_method[method])
    end
    
    # Add a dummy point at the end of each method's solve times to ensure they all have the same length for plotting
    max_length = maximum(length.(values(solve_times_per_method)))
    for method in keys(solve_times_per_method)
        push!(solve_times_per_method[method], tmax)
    end

    # Plot
    p = plot()
    max_number_solved = 0
    for method in keys(solve_times_per_method)
        number_solved = [1:length(solve_times_per_method[method]) - 1; length(solve_times_per_method[method]) - 1]
        max_number_solved = max(max_number_solved, length(solve_times_per_method[method]) - 1)
        plot!(solve_times_per_method[method], number_solved, label = method)
    end


    xlabel!("Solve Time (s)")
    ylabel!("Number of Instances Solved")
    plot!(legend = :topright, yticks = 0:2:max_number_solved, xlim = (0, tmax))

    suffix = only_opt_CG ? "_only_opt_CG" : ""

    savefig(p, "results/plots/solve_time_graph$suffix.png")

    plot!(xscale = :log10, legend = :topleft, yticks = 0:2:max_number_solved, xlim = (1e-2, tmax))
    savefig(p, "results/plots/solve_time_graph_log$suffix.png")
    return p
end

# -----------------------------
# Helper formatting functions
# -----------------------------

fmt_time(t) = isfinite(t) ? @sprintf("%.2f", t) : "--"
fmt_gap(g)  = isfinite(g) ? @sprintf("%.1f\\%%", 100g) : "--"
fmt_num(x)  = isfinite(x) ? @sprintf("%.2f", x) : "--"

function price_of_robustness(static_cost, robust_cost)
    if isfinite(static_cost) && isfinite(robust_cost) && static_cost > 0
        return (robust_cost - static_cost) / static_cost
    else
        return Inf
    end
end

# -----------------------------
# Main table generator
# -----------------------------

function json_to_latex_table(
    json_file::String,
    output_file::String;
    methods = ["cutting_planes", "cutting_planes_callbacks", "robust_dual", "CG-NC"]
)

    data = JSON.parsefile(json_file, allownan = true)
    instances = collect(keys(data))
    sort!(instances, by = x -> (parse(Int, split(x, "_")[1]), x)) # Sort by instance number

    open(output_file, "w") do io
        # ---- LaTeX header ----
        println(io, "\\begin{table}[ht]")
        println(io, "\\centering")
        println(io, "\\small")
        println(io, "\\begin{tabular}{l c" * " r r"^length(methods) * "}")
        println(io, "\\hline")
        println(io, "Instance & PR " *
                    join(["& \\multicolumn{2}{c}{$(replace(m, "_" => "\\_"))}" for m in methods], " ") *
                    " \\\\")
        println(io, " & " *
                    join(["& Time & Gap" for _ in methods], " ") *
                    " \\\\")
        println(io, "\\hline")

        # ---- Rows ----
        for instance in instances
            results = data[instance]
            if length(keys(results)) == 1
                println("Warning: Only static results for instance $instance. Skipping.")
                continue
            end
            if !haskey(results, "static")
                println("Warning: No static solution for instance $instance. Skipping.")
                continue
            end
            static_cost = get(results, "static", NaN)["cost"]

            # pick first robust method to compute PR
            robust_cost = Inf
            for m in methods
                if haskey(results, m)
                    robust_cost = min(robust_cost, results[m]["cost"])
                end
            end

            pr = price_of_robustness(static_cost, robust_cost)

            instance_formatted = replace(instance, "_" => "\\_")
            print(io, instance_formatted, " & ", fmt_gap(pr))

            for m in methods
                if haskey(results, m)
                    # Specific computation of the gap for CG-NC
                    gap = results[m]["gap"]
                    if m == "CG-NC" && haskey(results, "robust_dual")
                        opt_value = results["robust_dual"]["cost"]
                        gap = min(gap, price_of_robustness(opt_value, results[m]["cost"]))
                    end
                    r = results[m]
                    print(io,
                        " & ", fmt_time(r["time"]),
                        " & ", fmt_gap(abs(gap))
                    )
                else
                    print(io, " & -- & --")
                end
            end

            println(io, " \\\\")
        end

        # ---- Footer ----
        println(io, "\\hline")
        println(io, "\\end{tabular}")
        println(io, "\\caption{Comparison of solution methods}")
        println(io, "\\end{table}")
    end
end
