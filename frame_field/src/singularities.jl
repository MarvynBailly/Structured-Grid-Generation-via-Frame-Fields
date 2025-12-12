# Singularity computation for frame fields.
#
# This module implements the computation of singularities in frame fields
# according to the MIQ method (Ray et al. 2008, Bommes et al. 2009).
#
# For interior vertices, the cross field index is:
#     I(v) = Σ_{e ∈ star(v)} p_e / 4
#
# where p_e are the period jumps around the vertex in CCW order.
# Vertices with non-zero index are singularities and will correspond to
# irregular vertices in the final quad mesh.

using GeometryBasics
using LinearAlgebra

"""
    build_mesh_topology(mesh)

Build topology data structures for the mesh.

Returns a named tuple with:
- `edge_to_faces`: Dict mapping edge (v1, v2) -> [face_i, face_j]
- `dual_adj`: Dict mapping face_i -> [(neighbor_face, edge)]
- `vert_to_faces`: Dict mapping vertex -> [faces]
"""
function build_mesh_topology(mesh)
    fs = faces(mesh)
    n_faces = length(fs)
    
    # Build edge to faces mapping
    edge_to_faces = Dict{Tuple{Int,Int}, Vector{Int}}()
    dual_adj = Dict{Int, Vector{Tuple{Int, Tuple{Int,Int}}}}()
    
    for i in 1:n_faces
        dual_adj[i] = []
    end
    
    # Build edge map and dual adjacency
    for (i, face_i) in enumerate(fs)
        for (j, face_j) in enumerate(fs)
            if i >= j
                continue
            end
            
            verts_i = [face_i[1], face_i[2], face_i[3]]
            verts_j = [face_j[1], face_j[2], face_j[3]]
            shared = intersect(verts_i, verts_j)
            
            if length(shared) == 2
                v1, v2 = shared[1], shared[2]
                edge = v1 < v2 ? (v1, v2) : (v2, v1)
                
                if !haskey(edge_to_faces, edge)
                    edge_to_faces[edge] = [i, j]
                    
                    # Add to dual adjacency
                    push!(dual_adj[i], (j, edge))
                    push!(dual_adj[j], (i, edge))
                end
            end
        end
    end
    
    # Build vertex to faces mapping
    vert_to_faces = Dict{Int, Vector{Int}}()
    for (f_idx, face) in enumerate(fs)
        for v_idx in face
            if !haskey(vert_to_faces, v_idx)
                vert_to_faces[v_idx] = Int[]
            end
            push!(vert_to_faces[v_idx], f_idx)
        end
    end
    
    return (edge_to_faces=edge_to_faces, dual_adj=dual_adj, vert_to_faces=vert_to_faces)
end

"""
    get_ordered_1ring(pivot_v, adj_faces, all_faces, dual_adj)

Sorts the faces incident to `pivot_v` in Counter-Clockwise order 
by following the half-edge connectivity. Returns empty array if boundary vertex.
"""
function get_ordered_1ring(pivot_v::Int, adj_faces::Vector{Int}, all_faces, dual_adj)
    if isempty(adj_faces)
        return Int[]
    end
    
    # Start with an arbitrary face
    start_f = adj_faces[1]
    sorted = Int[start_f]
    
    # Set for fast lookup
    candidate_set = Set(adj_faces)
    curr_f = start_f
    
    # Walk CCW until we hit the start or a boundary
    for _ in 1:(length(adj_faces) + 1) # Safety limit
        # Find the "previous" vertex in current triangle relative to pivot
        face_verts = all_faces[curr_f]
        idx_in_tri = findfirst(==(pivot_v), face_verts)
        
        # Get the previous vertex (CCW from pivot)
        prev_idx_in_tri = mod1(idx_in_tri - 1, 3)
        v_prev = face_verts[prev_idx_in_tri]
        
        # Search neighbors for one that shares v_prev
        next_f = -1
        for (nbr, _) in dual_adj[curr_f]
            if nbr in candidate_set && v_prev in all_faces[nbr]
                next_f = nbr
                break
            end
        end
        
        if next_f == -1
            # Boundary reached
            return Int[]
        end
        
        if next_f == start_f
            # Loop closed!
            if length(sorted) != length(adj_faces)
                return Int[]
            end
            return sorted
        end
        
        push!(sorted, next_f)
        curr_f = next_f
    end
    
    return Int[] # Failed to close
end

"""
    find_shared_edge(f1, f2, dual_adj)

Find the edge shared between faces f1 and f2.
"""
function find_shared_edge(f1::Int, f2::Int, dual_adj)
    for (neighbor, edge) in dual_adj[f1]
        if neighbor == f2
            return edge
        end
    end
    return nothing
end

"""
    compute_singularities(mesh, period_jumps)

Compute singularities in the frame field using proper geometric ordering.

# Arguments
- `mesh`: The triangular mesh
- `period_jumps`: Dict mapping edge (v1, v2) to integer period jump value

# Returns
- `singularities`: Dict mapping vertex index to singularity index I(v)
                   Only interior vertices with non-zero index are included.

# Algorithm
For each interior vertex:
1. Sort incident faces in CCW order
2. Traverse the 1-ring, accumulating period jumps with proper orientation
3. Compute index = Σp / 4

Boundary vertices are automatically skipped (their 1-ring doesn't close).
"""
function compute_singularities(mesh, period_jumps::Dict{Tuple{Int,Int}, Int}; debug=false)
    singularities = Dict{Int, Float64}()
    
    # Build mesh topology
    topo = build_mesh_topology(mesh)
    fs = faces(mesh)
    
    # Debug: count non-zero period jumps
    if debug
        non_zero_jumps = count(p -> p != 0, values(period_jumps))
        println("  Total edges with period jumps: $(length(period_jumps))")
        println("  Non-zero period jumps: $non_zero_jumps")
    end
    
    # Iterate over all vertices
    n_interior = 0
    for (v_idx, adj_faces) in topo.vert_to_faces
        # Sort faces in CCW order around this vertex
        sorted_faces = get_ordered_1ring(v_idx, adj_faces, fs, topo.dual_adj)
        
        # If loop doesn't close, it's a boundary vertex -> skip
        if isempty(sorted_faces)
            continue
        end
        
        n_interior += 1
        
        # Sum period jumps around the vertex with proper orientation
        p_sum = 0
        
        for k in 1:length(sorted_faces)
            f_curr = sorted_faces[k]
            f_next = sorted_faces[mod1(k+1, length(sorted_faces))]
            
            # Find shared edge
            edge = find_shared_edge(f_curr, f_next, topo.dual_adj)
            if isnothing(edge)
                if debug
                    println("  Warning: Could not find edge between faces $f_curr and $f_next for vertex $v_idx")
                end
                continue
            end
            
            # Get period jump value for this edge
            p_val = get(period_jumps, edge, 0)
            
            # Determine orientation
            # edge_to_faces stores faces in consistent order
            ef_def = topo.edge_to_faces[edge]
            
            if ef_def == [f_curr, f_next]
                p_sum += p_val
            elseif ef_def == [f_next, f_curr]
                p_sum -= p_val
            else
                if debug
                    println("  Warning: Edge $edge has faces $ef_def but we're transitioning $f_curr -> $f_next")
                end
            end
        end
        
        # Compute singularity index
        # For interior vertices: I(v) = Σp / 4
        if p_sum != 0
            singularities[v_idx] = p_sum / 4.0
        end
    end
    
    if debug
        println("  Interior vertices checked: $n_interior")
        println("  Singularities found: $(length(singularities))")
    end
    
    return singularities
end
"""
    print_singularity_report(mesh, singularities)

Print a formatted report of singularities.

Shows the vertex index, singularity index, and predicted quad mesh valence.
"""
function print_singularity_report(mesh, singularities::Dict{Int, Float64})
    if isempty(singularities)
        println("✓ No singularities found - all vertices will be regular (valence 4) in quad mesh")
        return
    end
    
    # Compute total index
    total_index = sum(values(singularities))
    
    println("\n" * "="^70)
    println("SINGULARITY REPORT")
    println("="^70)
    println("Found $(length(singularities)) singularit$(length(singularities) == 1 ? "y" : "ies"):")
    println("Total index sum: $(round(total_index, digits=6)) (Gauss-Bonnet: χ/4 for closed surfaces)\n")
    
    # Sort by vertex index
    sorted_verts = sort(collect(keys(singularities)))
    
    for vertex_idx in sorted_verts
        index = singularities[vertex_idx]
        # Valence in quad mesh: valence = 4 - 4*index
        # But this is a simplified formula; more precisely:
        # index = +1/4 → valence 3 (one less edge)
        # index = -1/4 → valence 5 (one more edge)
        predicted_valence = 4 - round(Int, 4 * index)
        
        println("  Vertex $vertex_idx:")
        println("    Index:             $(round(index, digits=6))")
        println("    Predicted valence: $predicted_valence")
        
        # Get vertex position for reference
        coords = coordinates(mesh)
        pos = coords[vertex_idx]
        println("    Position:          ($(round(pos[1], digits=3)), $(round(pos[2], digits=3)))")
        println()
    end
    
    println("="^70)
end
