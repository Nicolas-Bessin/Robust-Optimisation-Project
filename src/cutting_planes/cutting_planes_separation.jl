include("../data/instance.jl")

using JuMP

"""
Create the initial length separation problem model - without the objective
"""
function cutting_planes_length_separation_model(data :: Data, gurobi_env)
    N = data.N
    model_SPL = Model(() -> Gurobi.Optimizer(gurobi_env))
    set_silent(model_SPL)

    @variable(model_SPL, 0 <= delta1[i = 1:N, j = i+1:N] <= data.delta_1_max)

    @constraint(model_SPL, 
        sum(delta1[i, j] for i in 1:N, j in i+1:N) <= data.L
    )

    return model_SPL
end

"""
Create the initial weights separation problem model - without the objective
"""
function cutting_planes_weights_separation_model(data :: Data, gurobi_env)
    N = data.N
    model_SPW = Model(() -> Gurobi.Optimizer(gurobi_env))
    set_silent(model_SPW)

    @variable(model_SPW, 0 <= delta2[i = 1:N] <= data.delta_2_max[i])

    @constraint(model_SPW, 
        sum(delta2[i] for i in 1:N) <= data.W
    )
    
    return model_SPW
end

"""
Solve the separation problem on length with the current values of x
"""
function separation_length(model_SPL, x, data :: Data, timelimit :: Int = -1)
    N = data.N
    if timelimit > 0
        set_time_limit_sec(model_SPL, timelimit)
    end

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
function separation_weights(model_SPW, y, k :: Int, data :: Data, timelimit :: Int = -1)
    N = data.N
    if timelimit > 0
        set_time_limit_sec(model_SPW, timelimit)
    end

    delta2 = model_SPW[:delta2]
    @objective(model_SPW, Max,
        sum( (1 + delta2[i]) * y[i, k] * data.weights[i] for i in 1:N )
    )

    optimize!(model_SPW)

    @assert is_solved_and_feasible(model_SPW)

    return objective_value(model_SPW), value.(delta2)
    return
end
