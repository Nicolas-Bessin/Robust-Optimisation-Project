using Plots

function plot_instance(data :: Data)
    les_x = Vector{Float64}()
    les_y = Vector{Float64}()

    for i in 1:data.N
        x,y = data.coordinates[i]
        push!(les_x,x)
        push!(les_y,y)
    end

    plt = scatter(les_x,les_y, markersize = data.weights)
    display(plt)

    return plt
end