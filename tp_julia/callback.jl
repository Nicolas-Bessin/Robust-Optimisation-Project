using JuMP
using Gurobi, GLPK
import Test
import MathOptInterface as MOI


function callbackMain()
    println("------------------------------")
    println("------------------------------")
    
    optimizer = optimizer_with_attributes(
           Gurobi.Optimizer, "Threads" => 1, "Presolve" => 0, "Method" => 1
       );

    m = Model(optimizer)

    # Il est imposé d'utiliser 1 seul thread en Julia avec CPLEX pour
    # utiliser les callbacks
    # MOI.set(m, Gurobi.Threads(), 1)
    # set_attribute(m, "Threads", 1)
    # set_attribute(m, "Presolve", 0)
    # set_attribute(m, "Method", 3)
    # set_silent(m)

    @variable(m, x, Int)
    @constraint(m, x <= 10)
    @objective(m, Max, x)

    lazy_cut_was_called = false
    calls = 0

    function lazy_cut_generic(cb_data)
        calls += 1
        lazy_cut_was_called = true

        status = callback_node_status(cb_data, m)
        if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
            # `callback_value(cb_data, x)` is not integer (to some tolerance).
            # If, for example, your lazy constraint generator requires an
            # integer-feasible primal solution, you can add a `return` here.
            return
        elseif status == MOI.CALLBACK_NODE_STATUS_INTEGER
            # `callback_value(cb_data, x)` is integer (to some tolerance).
        else
            @assert status == MOI.CALLBACK_NODE_STATUS_UNKNOWN
            # `callback_value(cb_data, x)` might be fractional or integer.
        end
        x_val = callback_value(cb_data, x)
        println("Callback was called with x = $x_val")
        if x_val > 1+ 1e-6
            con = @build_constraint(x <= 1)
            println("Adding constraint $con")
            MOI.submit(m, MOI.LazyConstraint(cb_data), con)
        end
        return
    end

    # On précise que le modèle doit utiliser notre fonction de callback 
    set_attribute(m, MOI.LazyConstraintCallback(), lazy_cut_generic)


    # function lazy_cut_gurobi(cb_data, cb_where :: Cint)
    #     calls += 1
    #     lazy_cut_was_called = true

    #     if cb_where != GRB_CB_MIPSOL
    #         return
    #     end

    #     # You can query a callback attribute using GRBcbget
    #     if cb_where == GRB_CB_MIPNODE
    #         resultP = Ref{Cint}()
    #         GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_STATUS, resultP)
    #         if resultP[] != GRB_OPTIMAL
    #             return  # Solution is something other than optimal.
    #         end
    #     end

    #     # Before querying `callback_value`, you must call:
    #     Gurobi.load_callback_variable_primal(cb_data, cb_where)
    #     x_val = callback_value(cb_data, x)

    #     if x_val > 1 + 1e-6
    #         con = @build_constraint(x <= 1)
    #         println("Currently, x = $x_val, Adding a constraint : $con")
    #         MOI.submit(m, MOI.LazyConstraint(cb_data), con)
    #     end
    #     return
    # end

    # # Solver dependent callback
    # MOI.set(m, MOI.RawOptimizerAttribute("LazyConstraints"), 1)
    # MOI.set(m, Gurobi.CallbackFunction(), lazy_cut_gurobi)

    optimize!(m)

    # Prints
    assert_is_solved_and_feasible(m)
    Test.@test lazy_cut_was_called
    Test.@test value(x) <= 1
    println("Optimal solution x = $(value(x))")
    println("Number of callback calls = $calls")
    return

end 

callbackMain()