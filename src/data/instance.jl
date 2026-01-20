
struct Data
    # Metadata
    instance_name :: String

    # Size
    N :: Int # Number of vertices
    K :: Int # Number of partition

    # Length-related data
    edge_lengths :: Matrix{Float64} # Static lengths
    L :: Int # Maximum cumulative relative increase on lengths
    l_hat :: Vector{Int}
    delta_1_max :: Int # Should be 3

    # Weight related data
    weights :: Vector{Int} # 
    B :: Int # Weight budget in a given partition
    W :: Int # Maximum cumulative relative increase on weights
    delta_2_max :: Vector{Float64} # Maximum delta_2 per vertex

    # Position data (only used for heuristics)
    coordinates :: Vector{Tuple{Float64, Float64}}
end




function parse_file(filepath :: String) :: Data
    current_section = nothing

    # Data
    n = 0
    K = 0

    coords :: Vector{Tuple{Float64, Float64}} = []
    L = 0.
    l_hat :: Vector{Int} = []
    delta_1_max = 3 # HARDCODED

    weights :: Vector{Int} = []
    B = 0.
    W = 0.
    delta_2_max :: Vector{Float64} = []

    lines = readlines(filepath)

    for line in lines
        line_data = nothing
        if contains(line, "=")
            # New section
            section, line_data = split(line, "=", limit = 2)
            current_section = strip(section)
        else
            # Not a new section, should only happen with coordinates
            line_data = line
        end

        if current_section == "n"
            n = parse(Int, line_data)
        elseif current_section == "L"
            L = parse(Int, line_data)
        elseif current_section == "W"
            W = parse(Int, line_data)
        elseif current_section == "K"
            K = parse(Int, line_data)
        elseif current_section == "B"
            B = parse(Int, line_data)
        elseif current_section == "w_v"
            line_data = strip(line_data)
            # [2:end-1] should eliminate '[' && ']'
            for substr in split(line_data[2:end-1], ",")
                push!(weights, parse(Int, substr))
            end
        elseif current_section == "W_v"
            line_data = strip(line_data)
            # [2:end-1] should eliminate '[' && ']'
            for substr in split(line_data[2:end-1], ",")
                push!(delta_2_max, parse(Float64, substr))
            end
        elseif current_section == "lh"
            line_data = strip(line_data)
            # [2:end-1] should eliminate '[' && ']'
            for substr in split(line_data[2:end-1], ",")
                push!(l_hat, parse(Int, substr))
            end
        elseif current_section == "coordinates"
            entries = split(strip(line_data), " ")
            if length(entries) == 1
                continue
            end
            @assert length(entries) == 3 # Two entries + ';' (or ']')
            x = parse(Float64, entries[1])
            y = parse(Float64, entries[2])
            push!(coords, (x, y))
        end
    end

    # Build the edge length matrix
    edge_lengths = zeros(n, n)
    for i ∈ 1:n
        for j ∈ i+1:n
            ci = coords[i]
            cj = coords[j]
            l = sqrt(
                sum( (ci .- cj).^2 )
            )
            edge_lengths[i, j] = l
            edge_lengths[j, i] = l
        end
    end

    return Data(
        basename(filepath),
        n,
        K,
        edge_lengths,
        L,
        l_hat,
        delta_1_max,
        weights,
        B,
        W,
        delta_2_max,
        coords
    )

end

# data = parse_file("data/202_gr_3.tsp")

# data.l_hat