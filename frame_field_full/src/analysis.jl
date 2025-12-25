module Analysis

using LinearAlgebra
using ..Types
using ..Topology
using Printf

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
    
    all_indices = Float64[]
    all_angle_defects = Float64[]
    all_k_sums = Float64[]
    all_p_sums = Float64[]
    boundary_vertices_skipped = 0
    degenerate_vertices_skipped = 0
    
    MIN_ANGLE_THRESHOLD = deg2rad(5.0)
    
    for v_idx in 1:length(topo.vertices)
        faces = v2f[v_idx]
        if isempty(faces) continue end
        
        star = Int[]; curr = faces[1]; push!(star, curr)
        while true
            curr = star[end]; found = false
            for (n, _) in topo.dual_adj[curr]
                if !(n in star) && (n in faces)
                    push!(star, n); found = true; break
                end
            end
            if !found || length(star) == length(faces); break; end
        end
        
        is_closed = false
        for (n, _) in topo.dual_adj[star[end]]
            if n == star[1]; is_closed = true; break; end
        end
        if !is_closed
            boundary_vertices_skipped += 1
            continue
        end
        
        ang_sum = 0.0; k_sum = 0.0; p_sum = 0.0
        n_star = length(star)
        min_angle = Inf
        
        for i in 1:n_star
            curr = star[i]; next = star[mod1(i+1, n_star)]
            
            tri = topo.faces[curr]; p = topo.vertices; c = p[v_idx]
            if tri[1]==v_idx; l=p[tri[2]]; r=p[tri[3]]
            elseif tri[2]==v_idx; l=p[tri[3]]; r=p[tri[1]]
            else; l=p[tri[1]]; r=p[tri[2]]; end
            u = (l.x-c.x, l.y-c.y, l.z-c.z); v = (r.x-c.x, r.y-c.y, r.z-c.z)
            angle_at_vertex = acos(clamp((u[1]*v[1]+u[2]*v[2]+u[3]*v[3])/(norm(u)*norm(v)), -1, 1))
            ang_sum += angle_at_vertex
            min_angle = min(min_angle, angle_at_vertex)
            
            k = haskey(field.transport_angles, (curr, next)) ? field.transport_angles[(curr,next)] : -field.transport_angles[(next,curr)]
            k_sum += k
            
            edge = minmax(curr, next)
            sign = (curr < next) ? 1.0 : -1.0
            
            if haskey(field.period_jumps, edge)
                raw_p = field.period_jumps[edge]
                p_sum += raw_p * sign
            end
        end
        
        angle_defect = 2π - ang_sum
        I = (angle_defect + k_sum)/(2π) + p_sum/4.0
        
        is_degenerate = (min_angle < MIN_ANGLE_THRESHOLD)
        if is_degenerate
            degenerate_vertices_skipped += 1
            continue
        end
        
        push!(all_indices, I)
        push!(all_angle_defects, angle_defect)
        push!(all_k_sums, k_sum)
        push!(all_p_sums, p_sum)
        
        if abs(I) > 0.1
            push!(singularities, (v_idx, I))
        end
    end
    
    if verbose
        println("\n=== SINGULARITY STATISTICS ===")
        println("Boundary vertices skipped: $boundary_vertices_skipped")
        println("Degenerate vertices skipped: $degenerate_vertices_skipped")
        println("Total valid interior vertices analyzed: $(length(all_indices))")
        println("Singularities found: $(length(singularities))")
    end
    
    return singularities
end

end