"""
    compute_singularities(mesh, topo, p_values)

Returns a Dict: Vertex Index -> Singularity Index (Float64)
Indices: +0.25 (Valence 3), -0.25 (Valence 5).
"""
function compute_singularities(mesh::Mesh, topo::MeshTopology, p_values::Dict{Int, Int})
    singularities = Dict{Int, Float64}()
    
    # 1. Build unordered adjacency (Vertex -> Faces)
    vert_to_faces = Dict{Int, Vector{Int}}()
    fs = faces(mesh)
    
    for (f_idx, face) in enumerate(fs)
        for v_idx in face
            if !haskey(vert_to_faces, v_idx)
                vert_to_faces[v_idx] = Int[]
            end
            push!(vert_to_faces[v_idx], f_idx)
        end
    end
    
    # 2. Iterate Interior Vertices
    for (v_idx, adj_faces) in vert_to_faces
        # Filter boundary quickly using edge counts if possible, 
        # but the sorting function will detect open loops anyway.
        
        # Sort Faces Geometrically CCW
        sorted_faces = get_ordered_1ring(v_idx, adj_faces, fs, topo)
        
        # If loop is not closed, it's a boundary vertex -> Skip
        if isempty(sorted_faces); continue; end 
        
        # 3. Sum Period Jumps around the loop
        p_sum = 0
        
        for k in 1:length(sorted_faces)
            f_curr = sorted_faces[k]
            f_next = sorted_faces[mod1(k+1, length(sorted_faces))]
            
            # Find shared edge
            e_idx = find_shared_edge(f_curr, f_next, topo)
            if e_idx == -1; continue; end
            
            # Determine Orientation
            # topo.edge_to_faces[e_idx] is always [f_min, f_max] or insertion order
            # We need to check if the transition f_curr -> f_next matches that order
            ef_def = topo.edge_to_faces[e_idx]
            p_val = get(p_values, e_idx, 0)
            
            if ef_def == [f_curr, f_next]
                p_sum += p_val
            elseif ef_def == [f_next, f_curr]
                p_sum -= p_val
            end
        end
        
        # 4. Identify Singularity
        # Valence = 4 - sum(p)
        # Index = sum(p) / 4
        if p_sum != 0
            singularities[v_idx] = p_sum / 4.0
        end
    end
    
    return singularities
end

# --- Robust Helpers ---

"""
    get_ordered_1ring(pivot_v, adj_faces, all_faces, topo)

Sorts the faces incident to `pivot_v` in Counter-Clockwise order 
by following the half-edge connectivity. Returns empty if boundary.
"""
function get_ordered_1ring(pivot_v, adj_faces, all_faces, topo)
    if isempty(adj_faces); return Int[]; end
    
    # Start with an arbitrary face
    start_f = adj_faces[1]
    sorted = Int[start_f]
    
    # We map faces to a Set for fast lookup
    candidate_set = Set(adj_faces)
    
    curr_f = start_f
    
    # Walk CCW until we hit the start or a boundary
    for _ in 1:(length(adj_faces) + 1) # Safety limit
        # Find the "Previous" vertex in the current triangle relative to pivot
        # Triangle: (pivot, next, prev) in CCW order
        # The face sharing (pivot, prev) is the CCW neighbor.
        
        face_verts = all_faces[curr_f]
        idx_in_tri = findfirst(==(pivot_v), face_verts) # 1, 2, or 3
        
        # In a standard CCW triangle (v1, v2, v3):
        # Edge (v1, v2) is outgoing. Edge (v3, v1) is incoming.
        # The neighbor across the "incoming" edge is the next face in the CCW fan.
        # Incoming edge connects: (v_prev, pivot)
        
        prev_idx_in_tri = mod1(idx_in_tri - 1, 3)
        v_prev = face_verts[prev_idx_in_tri]
        
        # Search neighbors for one that shares v_prev
        next_f = -1
        for (nbr, _) in topo.dual_adj[curr_f]
            # Must be in our 1-ring and share the specific edge vertex
            if nbr in candidate_set && v_prev in all_faces[nbr]
                 next_f = nbr
                 break
            end
        end
        
        if next_f == -1
            # Boundary reached (no face on the left)
            return Int[]
        end
        
        if next_f == start_f
            # Loop Closed!
            if length(sorted) != length(adj_faces)
                # We closed a loop but missed some faces? Non-manifold geometry.
                return Int[]
            end
            return sorted
        end
        
        push!(sorted, next_f)
        curr_f = next_f
    end
    
    return Int[] # Failed to close
end

function find_shared_edge(f1, f2, topo)
    for (neighbor, e_idx) in topo.dual_adj[f1]
        if neighbor == f2
            return e_idx
        end
    end
    return -1
end