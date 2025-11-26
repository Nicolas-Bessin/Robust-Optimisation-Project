include("instance.jl")

using JuMP, Gurobi

"""
Solve the separation problem on length with the current values of x
"""
function separation_length(model_SPL, x, data :: Data)
    N = data.N

    delta1 = model_SPL[:delta1]
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
function separation_weights(model_SPW, y, k :: Int, data :: Data)
    N = data.N

    delta2 = model_SPW[:delta2]
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

    # Prepare the separation problems 
    model_SPL = Model(Gurobi.Optimizer)
    set_silent(model_SPL)

    @variable(model_SPL, 0 <= delta1[i = 1:N, j = i+1:N] <= data.delta_1_max)

    @constraint(model_SPL, 
        sum(delta1[i, j] for i in 1:N, j in i+1:N) <= data.L
    )

    models_SPW :: Vector{JuMP.GenericModel{Core.Float64}} = []
    for k in 1:K
        model_SPW_k = Model(Gurobi.Optimizer)
        set_silent(model_SPW_k)

        @variable(model_SPW_k, 0 <= delta2[i = 1:N] <= data.delta_2_max[i])

        # @constraint(model_SPW, [i = 1:N], 0 <= delta2[i] <= data.delta_2_max[i])
        @constraint(model_SPW_k, 
            sum(delta2[i] for i in 1:N) <= data.W
        )
        push!(models_SPW, model_SPW_k)
    end

    optimality_reached = false
    iter = 0
    cutting_planes_count = 0

    while !optimality_reached && iter < ITMAX
        # Re-optimize model
        optimize!(model)
        @assert is_solved_and_feasible(model)
        optimality_reached = true # we will set it back to false as needed
        x_val = value.(x)
        y_val = value.(y)

        # Solve both separation problems, update optimality and add the cuts
        # 1) Length separation
        val, delta1 = separation_length(model_SPL, x_val, data)
        if val > value.(z) + eps 
            optimality_reached = false
            # Add the length cut
            println("Adding a length cut")
            @constraint(model, 
                sum(x[i, j] * delta1[i, j] * (data.l_hat[i] + data.l_hat[j]) for i in 1:N, j in i+1:N) <= z
            )
            cutting_planes_count += 1
        end
        # 2) Weight separation
        for k in 1:K
            val, delta2 = separation_weights(models_SPW[k], y_val, k, data)
            if val > data.B
                optimality_reached = false
                # Add the weight cut for all values of k !
                println("Adding a weight cut")
                @constraint(model, [k = 1:K],
                    sum( (1 + delta2[i]) * y[i, k] * data.weights[i] for i in 1:N) <= data.B
                )
                cutting_planes_count += K
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
    println("Added a total of $cutting_planes_count cutting planes")
end

data = parse_file("data/10_ulysses_3.tsp");

@time cutting_planes_method(data)