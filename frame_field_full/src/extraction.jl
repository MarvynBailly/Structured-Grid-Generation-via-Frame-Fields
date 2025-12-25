module Extraction

using ..Types
using Printf

export extract_quad_mesh_isolines, extract_quad_mesh

function extract_quad_mesh_isolines(topo::MeshTopology, u_coords::Vector{Float64}, v_coords::Vector{Float64}; 
                                   grid_spacing=nothing, verbose=false)
    println("\n--- Extracting Quad Mesh via Isoline Tracing ---")
    
    u_min, u_max = minimum(u_coords), maximum(u_coords)
    v_min, v_max = minimum(v_coords), maximum(v_coords)
    u_range = u_max - u_min
    v_range = v_max - v_min
    
    println("U range: [", u_min, ", ", u_max, "], span: ", u_range)
    println("V range: [", v_min, ", ", v_max, "], span: ", v_range)
    
    if grid_spacing === nothing
        target_divs = 8
        u_spacing = u_range / target_divs
        v_spacing = v_range / target_divs
    else
        u_spacing = grid_spacing
        v_spacing = grid_spacing
    end
    
    println("Grid spacing: u=", u_spacing, ", v=", v_spacing)
    
    grid_verts = Dict{Tuple{Int,Int}, Vector{Int}}()
    
    for v_idx in 1:length(u_coords)
        u_grid_idx = round(Int, (u_coords[v_idx] - u_min) / u_spacing)
        v_grid_idx = round(Int, (v_coords[v_idx] - v_min) / v_spacing)
        key = (u_grid_idx, v_grid_idx)
        
        if !haskey(grid_verts, key)
            grid_verts[key] = Int[]
        end
        push!(grid_verts[key], v_idx)
    end
    
    println("Mapped vertices to ", length(grid_verts), " grid cells")
    
    quads = Tuple{Int,Int,Int,Int}[]
    
    all_cells = collect(keys(grid_verts))
    
    for (u_g, v_g) in all_cells
        corners = [
            (u_g, v_g),
            (u_g+1, v_g),
            (u_g+1, v_g+1),
            (u_g, v_g+1)
        ]
        
        if all(haskey(grid_verts, c) && !isempty(grid_verts[c]) for c in corners)
            v1 = grid_verts[corners[1]][1]
            v2 = grid_verts[corners[2]][1]
            v3 = grid_verts[corners[3]][1]
            v4 = grid_verts[corners[4]][1]
            
            if length(Set([v1, v2, v3, v4])) == 4
                push!(quads, (v1, v2, v3, v4))
            end
        end
    end
    
    println("Extracted ", length(quads), " quadrilaterals")
    
    return quads, (u_spacing, v_spacing)
end

function extract_quad_mesh(topo::MeshTopology, u_coords::Vector{Float64}, v_coords::Vector{Float64}; 
                          grid_spacing=nothing, verbose=false)
    return extract_quad_mesh_isolines(topo, u_coords, v_coords; grid_spacing=grid_spacing, verbose=verbose)
end

end