include("instance.jl")

using JuMP, Gurobi

function static_problem(data :: Data)

    model = Model(Gurobi.Optimizer)

    N = data.N
    K = data.K

    @variable(model, x[i = 1:N, j = 1:N], Bin)
    @variable(model, y[i = 1:N, k = 1:K], Bin)

    @constraint(model, [i = 1:N, j = 1:N, k = 1:K],
        x[i,j] >= y[i,k] + y[j, k] - 1
    )
    @constraint(model, [i = 1:N],
        sum(y[i,k] for k in 1:K) == 1
    )
    @constraint(model, [k in 1:K],
        sum(y[i,k] * data.weights[i] for i in 1:N) <= data.B
    )

    @objective(model, Min,
        sum(x[i,j] * data.edge_lengths[i,j] for i in 1:N, j in i+1:N)
    )

    optimize!(model)

    @assert is_solved_and_feasible(model)

    x_val = value.(x)
    y_val = value.(y)

    # Re-build the partitions in a printable way
    partitions :: Vector{Vector{Int}} = repeat([[]], K)
    # println(length(partitions))
    for i in 1:N
        k = findfirst(y -> y == 1, y_val[i, :])
        # println("Vertex $i is in partition $k")
        push!(partitions[k], i)
    end

    partitions = filter(p -> !isempty(p), partitions)
    println("Partition is $partitions")

end

data = parse_file("data/10_ulysses_3.tsp")

static_problem(data)