include("master.jl")
include("../heuristic/greedy.jl")

function CG_solver_non_compact(
    data :: Data;
    timelimit :: Int,
    itermax :: Int
)

    initial_sol = greedy_init(data)
    master = master_problem(data, initial_sol)

end