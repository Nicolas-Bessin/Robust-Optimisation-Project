include("instance.jl")

using JuMP, Gurobi

"""
Solve the separation problem on length with the current values of x
"""
function separation_length(x, data :: Data)
    N = data.N

    model_SPL = Model(Gurobi.Optimizer)
    set_silent(model_SPL)

    @variable(model_SPL, 0 <= delta1[i = 1:N, j = i+1:N] <= data.delta_1_max)

    @constraint(model_SPL, 
        sum(delta1[i, j] for i in 1:N, j in i+1:N) <= data.L
    )

    @objective(model_SPL, Max,
        sum( delta1[i, j] * x[i, j] * (data.l_hat[i] + data.l_hat[j]) for i in 1:N, j in i+1:N )
    )

    optimize!(model_SPL)

    @assert is_solved_and_feasible(model_SPL)

    return objective_value(model_SPL), value.(delta1)
end

"""
Solve the separation problem on length with the current values of y for partition k
"""
function separation_weights(y, k :: Int, data :: Data)
    N = data.N

    model_SPW = Model(Gurobi.Optimizer)
    set_silent(model_SPW)

    @variable(model_SPW, 0 <= delta2[i = 1:N] <= data.delta_2_max[i])

    # @constraint(model_SPW, [i = 1:N], 0 <= delta2[i] <= data.delta_2_max[i])
    @constraint(model_SPW, 
        sum(delta2[i] for i in 1:N) <= data.W
    )

    @objective(model_SPW, Max,
        sum( (1 + delta2[i]) * y[i, k] * data.weights[i] for i in 1:N )
    )

    optimize!(model_SPW)

    @assert is_solved_and_feasible(model_SPW)

    return objective_value(model_SPW), value.(delta2)
    return
end


function cutting_planes_method(data :: Data, ITMAX :: Int = 10, eps :: Float64 = 1e-6)
    N = data.N
    K = data.K

    model = Model(Gurobi.Optimizer)
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
    # @constraint(model, [k in 1:K],
    #     sum(y[i,k] * data.weights[i] for i in 1:N) <= data.B
    # )

    @objective(model, Min,
        sum(x[i,j] * data.edge_lengths[i,j] for i in 1:N, j in i+1:N) + z
    )

    optimality_reached = false
    iter = 0

    while !optimality_reached && iter < ITMAX
        # Re-optimize model
        optimize!(model)
        @assert is_solved_and_feasible(model)
        optimality_reached = true # we will set it back to false as needed
        x_val = value.(x)
        y_val = value.(y)

        # Solve both separation problems, update optimality and add the cuts
        # 1) Length separation
        val, delta1 = separation_length(x_val, data)
        if val > value.(z) + eps 
            optimality_reached = false
            # Add the length cut
            println("Adding a length cut")
            @constraint(model, 
                sum(x[i, j] * delta1[i, j] * (data.l_hat[i] + data.l_hat[j]) for i in 1:N, j in i+1:N) <= z
            )
        end
        # 2) Weight separation
        for k in 1:K
            val, delta2 = separation_weights(y_val, k, data)
            if val > data.B
                optimality_reached = false
                # Add the weight cut for all values of k !
                println("Adding a weight cut")
                @constraint(model, [k = 1:K],
                    sum( (1 + delta2[i]) * y[i, k] * data.weights[i] for i in 1:N) <= data.B
                )
            end
        end
    end

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
    println("Objective value is $(objective_value(model))")
end

data = parse_file("data/10_ulysses_3.tsp");

cutting_planes_method(data)