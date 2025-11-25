using JuMP
using CPLEX

function callbackMain()
    
    m = Model(CPLEX.Optimizer)

    # Il est imposé d'utiliser 1 seul thread en Julia avec CPLEX pour
    # utiliser les callbacks
    MOI.set(m, MOI.NumberOfThreads(), 1)

    @variable(m, x, Int)
    @constraint(m, x <= 10)
    @objective(m, Max, x)

    # Fonction qui sera exécutée par CPLEX à chaque fois qu'une des conditions
    # suivantes est remplie :
    # - 1 solution entière a été trouvée ;
    # - 1 relaxation a été calculée ;
    # - ...
    #
    # Arguments :
    # - context_id : permet de déterminer pour quelle raison le callback à
    # été appelé ;
    # - cb_data  :   permet  d'obtenir   d'autres   informations
    #   (valeur des bornes inférieures et supérieures , meilleure solution
    #   connue, ...)
    function mon_super_callback(cb_data::CPLEX.CallbackContext, context_id::Clong)

        # On teste d'abord si c'est l'obtention d'une solution entière
        # qui a entraîné l'appel du callback
        # (cette fonction isIntegerPoint est définie ci-dessous mais son 
        # contenu n'est pas très important)
        if isIntegerPoint(cb_data, context_id)
            
            # Cette ligne doit être  appelée avant de pouvoir récupérer la
            # solution entière ayant entraîné l'appel du callback
            CPLEX.load_callback_variable_primal(cb_data, context_id)

            # On récupère la valeur de x
            x_val = callback_value(cb_data, x)

            # Si elle est plus grande que 1, on ajoute la contrainte x <= 1
            if x_val >= 2 - 1e-5
                cstr = @build_constraint(x <= 1)
                MOI.submit(m, MOI.LazyConstraint(cb_data), cstr)
                println("Add constraint x <= 1")
            end
        end
    end

    # On précise que le modèle doit utiliser notre fonction de callback
    MOI.set(m, CPLEX.CallbackFunction(), mon_super_callback)
    optimize!(m)
end 

# Fonction  permettant  de  déterminer  si  c'est  l'obtention  d'une
# solution entière qui a entraîné l'appel d'un callback
# (il n'est pas nécessaire d'en comprendre le fonctionnement)
function isIntegerPoint(cb_data::CPLEX.CallbackContext, context_id::Clong)

    # context_id  == CPX_CALLBACKCONTEXT_CANDIDATE si le  callback est
    # appelé dans un des deux cas suivants :
    # cas 1 - une solution entière a été obtenue; ou
    # cas 2 - une relaxation non bornée a été obtenue
    if context_id != CPX_CALLBACKCONTEXT_CANDIDATE
        return false
    end

    # Pour déterminer si on est dans le cas 1 ou 2, on essaie de récupérer la
    # solution entière courante
    ispoint_p = Ref{Cint}()
    ret = CPXcallbackcandidateispoint(cb_data, ispoint_p)

    # S'il n'y a pas de solution entière
    if ret != 0 || ispoint_p[] == 0
        return false
    else
        return true
    end
end
