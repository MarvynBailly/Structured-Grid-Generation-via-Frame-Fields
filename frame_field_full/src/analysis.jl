module Analysis

using LinearAlgebra
using ..Types
using ..Topology
using Printf
using DataStructures

export build_vertex_to_faces, compute_singularities

function build_vertex_to_faces(topo::MeshTopology)
    n_verts = length(topo.vertices)
    v2f = [Int[] for _ in 1:n_verts]
    for (f_idx, tri) in enumerate(topo.faces)
        push!(v2f[tri[1]], f_idx)
        push!(v2f[tri[2]], f_idx)
        push!(v2f[tri[3]], f_idx)
    end
    return v2f
end

function compute_singularities(field::CrossField; verbose=false)
    if verbose
        println("\n--- Computing Singularities (Normalized) ---")
    end
    topo = field.topology
    v2f = build_vertex_to_faces(topo)
    singularities = Tuple{Int, Float64}[]
    
    # Statistics
    boundary_vertices_skipped = 0
    degenerate_vertices_skipped = 0
    
    MIN_ANGLE_THRESHOLD = deg2rad(5.0)
    
    for v_idx in 1:length(topo.vertices)
        faces = v2f[v_idx]
        if isempty(faces) continue end
        
        # --- Sort faces in CCW order ---
        # 1. Pick an arbitrary start face
        start_face = faces[1]
        
        # 2. Traverse CCW ("Next" neighbor is across the edge entering v)
        # 3. Traverse CW ("Prev" neighbor is across the edge leaving v)
        
        deque = Deque{Int}()
        push!(deque, start_face)
        visited = Set{Int}([start_face])
        
        # Go CCW (Forward)
        curr = start_face
        while true
            # Find the 'incoming' edge to v in 'curr'
            tri = topo.faces[curr]
            # Find index of v in tri (1, 2, or 3)
            idx = -1
            if tri[1] == v_idx; idx = 1
            elseif tri[2] == v_idx; idx = 2
            elseif tri[3] == v_idx; idx = 3; end
            
            # In CCW winding (1-2-3), the edge *entering* index 1 is (3,1).
            # The edge *entering* idx is (prev_idx, idx).
            prev_idx = mod1(idx - 1, 3)
            v_prev = tri[prev_idx]
            
            # The shared edge key
            edge_key = minmax(v_idx, v_prev)
            
            # Find the neighbor across this edge
            next_face = -1
            for (n, key) in topo.dual_adj[curr]
                if key == edge_key && n != curr
                    next_face = n
                    break
                end
            end
            
            if next_face != -1 && !(next_face in visited) && (next_face in faces)
                push!(deque, next_face)
                push!(visited, next_face)
                curr = next_face
            else
                break # Hit boundary or looped back
            end
            
            if curr == start_face; break; end # Closed loop
        end
        
        # If not closed, we might have started in the middle of a strip. Go CW (Backward).
        # Check if the CCW search wrapped around to the start (Closed loop condition)
        
        last_face_in_deque = last(deque) # Replaced 'back' with 'last'
        tri = topo.faces[last_face_in_deque]
        idx = findfirst(isequal(v_idx), tri)
        prev_idx = mod1(idx - 1, 3)
        edge_key = minmax(v_idx, tri[prev_idx])
        
        closed_loop = false
        for (n, key) in topo.dual_adj[last_face_in_deque]
            if key == edge_key && n == first(deque) # Replaced 'front' with 'first'
                closed_loop = true
                break
            end
        end

        if !closed_loop
            # Go CW from start to find the other side of the fan
            curr = start_face
            while true
                # Find 'outgoing' edge from v in 'curr' (v -> next vertex)
                tri = topo.faces[curr]
                idx = findfirst(isequal(v_idx), tri)
                next_vert_idx = mod1(idx + 1, 3)
                v_next = tri[next_vert_idx]
                
                edge_key = minmax(v_idx, v_next)
                
                prev_face = -1
                for (n, key) in topo.dual_adj[curr]
                    if key == edge_key && n != curr
                        prev_face = n
                        break
                    end
                end
                
                if prev_face != -1 && !(prev_face in visited) && (prev_face in faces)
                    pushfirst!(deque, prev_face)
                    push!(visited, prev_face)
                    curr = prev_face
                else
                    break
                end
            end
        end
        
        star = collect(deque)
        
        # Check boundary condition again properly
        if !closed_loop
            boundary_vertices_skipped += 1
            continue
        end
        
        # --- End Sort Fix ---
        
        ang_sum = 0.0; k_sum = 0.0; p_sum = 0.0
        n_star = length(star)
        min_angle = Inf
        
        for i in 1:n_star
            curr = star[i]
            next = star[mod1(i+1, n_star)]
            
            tri = topo.faces[curr]
            p = topo.vertices; c = p[v_idx]
            
            if tri[1]==v_idx
                l=p[tri[2]]; r=p[tri[3]]
            elseif tri[2]==v_idx
                l=p[tri[3]]; r=p[tri[1]]
            else
                l=p[tri[1]]; r=p[tri[2]]
            end
            
            u = (l.x-c.x, l.y-c.y, l.z-c.z)
            v = (r.x-c.x, r.y-c.y, r.z-c.z)
            
            angle_at_vertex = acos(clamp((u[1]*v[1]+u[2]*v[2]+u[3]*v[3])/(norm(u)*norm(v)), -1, 1))
            ang_sum += angle_at_vertex
            min_angle = min(min_angle, angle_at_vertex)
            
            # Transport Angle Sum
            k = haskey(field.transport_angles, (curr, next)) ? field.transport_angles[(curr,next)] : -field.transport_angles[(next,curr)]
            k_sum += k
            
            # Period Jump Sum
            edge = minmax(curr, next)
            sign = (curr < next) ? 1.0 : -1.0
            
            if haskey(field.period_jumps, edge)
                p_sum += field.period_jumps[edge] * sign
            end
        end
        
        angle_defect = 2π - ang_sum
        I = (angle_defect + k_sum)/(2π) + p_sum/4.0
        
        # if min_angle < MIN_ANGLE_THRESHOLD
        #     degenerate_vertices_skipped += 1
        #     continue
        # end
        
        if abs(I) > 0.1
            push!(singularities, (v_idx, I))
        end
    end
    
    if verbose
        println("\n=== SINGULARITY STATISTICS ===")
        println("Boundary vertices skipped: $boundary_vertices_skipped")
        println("Degenerate vertices skipped: $degenerate_vertices_skipped")
        println("Total valid interior vertices analyzed: $(length(topo.vertices) - boundary_vertices_skipped - degenerate_vertices_skipped)")
        println("Singularities found: $(length(singularities))")
        for (v_idx, I) in singularities
            @printf("  Vertex %d: Index = %.4f\n", v_idx, I)
        end
        println("================================\n")
    end
    
    return singularities
end

end