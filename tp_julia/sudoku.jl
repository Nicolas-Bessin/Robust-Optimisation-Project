using JuMP
using CPLEX
# using GLPK

function sudoku(t::Matrix{Int})

    # Taille de la grille
    n = size(t, 1)

    # Créer le modèle
    m = Model(CPLEX.Optimizer)
    # m = Model(GLPK.Optimizer)

    # Limite le temps d'exécution à 60 secondes
    set_time_limit_sec(m, 60.0)

    ### Variables
    # x[i, j, k] = 1 if cell (i, j) has value k
    @variable(m, x[1:n, 1:n, 1:n], Bin)

    ### Objectif : maximiser la valeur de la case en haut à gauche
    # L'objectif est quelconque car on souhaite simplement trouver une solution faisable
    @objective(m, Max, sum(k * x[1, 1, k] for k in 1:n))

    ### Contraintes
    # Obtenir la taille d'un bloque (3 dans une grille standard)
    blockSize = round.(Int, sqrt(n))

    # ... TODO : Ajouter les contraintes du modèle
    # Remarque : un 0 dans la matrice t, indique que la valeur de la case correspondante de la grille de sudoku n'est pas fixée

    ### Résoudre le problème
    optimize!(m)

    ### Si une solution a été trouvée dans le temps limite
    if primal_status(m) == MOI.FEASIBLE_POINT

        # Afficher la valeur de l'objectif
        objectiveValue = round(Int, JuMP.objective_value(m))
        println("Valeur de l'objectif : ", round(Int, JuMP.objective_value(m)))

        # ... TODO : Afficher la solution
        #
        # Attention : les calculs machine sont approchés et une variable entière z pourrait avoir la valeur 0.9999 au lieu de 1
        # Ainsi si vous voulez tester si z est égale à 1 dans la solution retournée par un solveur, il ne faut pas écrire
        #     if JuMP.value(z) == 1
        # mais plutôt
        #     if JuMP.value(z) > 1 - 1E-4
        # ou
        #     if round(JuMP.value(z)) == 1

        # Si la solution optimale n'est pas obtenue
        if !termination_status(m) == MOI.OPTIMAL

            # Le solveur fournit la meilleure borne supérieure connue sur la solution optimale
            bound = JuMP.objective_bound(m)

            # On peut utiliser cette borne pour savoir à quel point la solution trouvée par le solveur est bonne.
            # Pour cela, on calcule le gap (i.e., l'écart relatif entre l'objectif de la solution trouvée et la borne supérieure).
            gap = 100 * (objectiveValue - bound) / objectiveValue
            println("Gap : ", gap, "%")
        end 
    else
        println("Aucun solution trouvée dans le temps imparti.")
    end   
end

t = [
0 7 0 2 0 3 0 1 0;
3 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 2 0 0;
0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 2;
2 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0]

sudoku(t)
