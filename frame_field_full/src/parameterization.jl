module Parameterization

using ..Types
using SparseArrays
using LinearAlgebra

export compute_global_parameterization, smooth_parameterization, compute_parameterization_least_squares

function compute_global_parameterization(field::CrossField, rotations::Vector{Int}, cut_edges::Set{Tuple{Int,Int}})
    println("\n--- Computing Global Parameterization (u, v) ---")
    
    topo = field.topology
    n_verts = length(topo.vertices)
    
    edges = Tuple{Int,Int}[]
    edge_set = Set{Tuple{Int,Int}}()
    
    for tri in topo.faces
        for (v1, v2) in [(tri[1], tri[2]), (tri[2], tri[3]), (tri[3], tri[1])]
            key = minmax(v1, v2)
            if !(key in edge_set) && !(key in cut_edges)
                push!(edge_set, key)
                push!(edges, key)
            end
        end
    end
    
    println("Building least-squares system for $(length(edges)) edges...")
    
    I_rows = Int[]
    J_cols = Int[]
    vals = Float64[]
    b_u = Float64[]
    b_v = Float64[]
    
    equation_idx = 0
    
    for (v1, v2) in edges
        face_idx = 0
        for (f_idx, tri) in enumerate(topo.faces)
            if (v1 in tri) && (v2 in tri)
                face_idx = f_idx
                break
            end
        end
        
        if face_idx == 0
            continue
        end
        
        re = topo.face_ref_edges[face_idx]
        p1 = topo.vertices[re[1]]
        p2 = topo.vertices[re[2]]
        ref_ang = atan(p2.y - p1.y, p2.x - p1.x)
        
        local_theta = field.theta[face_idx]
        r = rotations[face_idx]
        global_angle = ref_ang + local_theta + (r * π / 2.0)
        
        u_dir = (cos(global_angle), sin(global_angle))
        v_dir = (-sin(global_angle), cos(global_angle))
        
        p_v1 = topo.vertices[v1]
        p_v2 = topo.vertices[v2]
        edge_vec = (p_v2.x - p_v1.x, p_v2.y - p_v1.y)
        edge_len = sqrt(edge_vec[1]^2 + edge_vec[2]^2)
        
        if edge_len < 1e-10
            continue
        end
        
        target_du = (edge_vec[1] * u_dir[1] + edge_vec[2] * u_dir[2])
        target_dv = (edge_vec[1] * v_dir[1] + edge_vec[2] * v_dir[2])
        
        equation_idx += 1
        push!(I_rows, equation_idx); push!(J_cols, v1); push!(vals, -1.0)
        push!(I_rows, equation_idx); push!(J_cols, v2); push!(vals, 1.0)
        push!(b_u, target_du)
        push!(b_v, target_dv)
    end
    
    equation_idx += 1
    push!(I_rows, equation_idx); push!(J_cols, 1); push!(vals, 1e3)
    push!(b_u, 0.0)
    push!(b_v, 0.0)
    
    A = sparse(I_rows, J_cols, vals, equation_idx, n_verts)
    
    println("Solving least-squares system ($equation_idx equations, $n_verts unknowns)...")
    
    ATA = A' * A
    ATb_u = A' * b_u
    ATb_v = A' * b_v
    
    u_coords = ATA \ ATb_u
    v_coords = ATA \ ATb_v
    
    println("Computed smooth parameterization for all $n_verts vertices")
    
    return u_coords, v_coords
end

function smooth_parameterization(topo::MeshTopology, u_coords::Vector{Float64}, v_coords::Vector{Float64}, 
                                 cut_edges::Set{Tuple{Int,Int}}; iterations=5)
    println("\n--- Smoothing Parameterization ---")
    
    n_verts = length(u_coords)
    u_smooth = copy(u_coords)
    v_smooth = copy(v_coords)
    
    vert_neighbors = [Int[] for _ in 1:n_verts]
    for face in topo.faces
        edges = [(face[1], face[2]), (face[2], face[3]), (face[3], face[1])]
        for (v1, v2) in edges
            edge_key = minmax(v1, v2)
            if !(edge_key in cut_edges)
                if !(v2 in vert_neighbors[v1])
                    push!(vert_neighbors[v1], v2)
                end
                if !(v1 in vert_neighbors[v2])
                    push!(vert_neighbors[v2], v1)
                end
            end
        end
    end
    
    boundary_verts = Set{Int}()
    edge_counts = Dict{Tuple{Int,Int}, Int}()
    for face in topo.faces
        for i in 1:3
            v1, v2 = face[i], face[mod1(i+1, 3)]
            key = minmax(v1, v2)
            edge_counts[key] = get(edge_counts, key, 0) + 1
        end
    end
    for (edge, count) in edge_counts
        if count == 1
            push!(boundary_verts, edge[1])
            push!(boundary_verts, edge[2])
        end
    end
    
    for iter in 1:iterations
        u_new = copy(u_smooth)
        v_new = copy(v_smooth)
        
        for v_idx in 1:n_verts
            if v_idx in boundary_verts || isempty(vert_neighbors[v_idx])
                continue
            end
            
            neighbors = vert_neighbors[v_idx]
            u_avg = sum(u_smooth[n] for n in neighbors) / length(neighbors)
            v_avg = sum(v_smooth[n] for n in neighbors) / length(neighbors)
            
            weight = 0.5
            u_new[v_idx] = (1 - weight) * u_smooth[v_idx] + weight * u_avg
            v_new[v_idx] = (1 - weight) * v_smooth[v_idx] + weight * v_avg
        end
        
        u_smooth = u_new
        v_smooth = v_new
    end
    
    println("Applied $iterations smoothing iterations")
    return u_smooth, v_smooth
end

function compute_parameterization_least_squares(field, rotations, cut_edges)
    println("\n--- Solving Global Parameterization (Least Squares) ---")
    
    topo = field.topology
    n_verts = length(topo.vertices)
    
    est_edges = n_verts * 3
    I = Int[]; J = Int[]; V = Float64[]
    b_u = zeros(est_edges)
    b_v = zeros(est_edges)
    
    row_idx = 0
    
    cut_set = Set{Tuple{Int,Int}}()
    for (u, v) in cut_edges
        push!(cut_set, minmax(u, v))
    end
    
    seen_edges = Set{Tuple{Int,Int}}()
    
    for i in 1:length(topo.faces)
        tri = topo.faces[i]
        
        re = topo.face_ref_edges[i]
        p1_ref = topo.vertices[re[1]]
        p2_ref = topo.vertices[re[2]]
        ref_ang = atan(p2_ref.y - p1_ref.y, p2_ref.x - p1_ref.x)
        
        global_angle = ref_ang + field.theta[i] + (rotations[i] * π / 2.0)
        
        dir_u = (cos(global_angle), sin(global_angle))
        dir_v = (-sin(global_angle), cos(global_angle))
        
        face_edges = [(tri[1], tri[2]), (tri[2], tri[3]), (tri[3], tri[1])]
        
        for (u, v) in face_edges
            key = minmax(u, v)
            
            if key in cut_set
                continue
            end
            
            if key in seen_edges
                continue
            end
            push!(seen_edges, key)
            
            row_idx += 1
            
            p_u = topo.vertices[u]
            p_v = topo.vertices[v]
            dx = p_v.x - p_u.x
            dy = p_v.y - p_u.y
            
            target_du = dx * dir_u[1] + dy * dir_u[2]
            target_dv = dx * dir_v[1] + dy * dir_v[2]
            
            push!(I, row_idx); push!(J, v); push!(V, 1.0)
            push!(I, row_idx); push!(J, u); push!(V, -1.0)
            
            if row_idx > length(b_u)
                resize!(b_u, length(b_u)*2)
                resize!(b_v, length(b_v)*2)
            end
            b_u[row_idx] = target_du
            b_v[row_idx] = target_dv
        end
    end
    
    row_idx += 1
    push!(I, row_idx); push!(J, 1); push!(V, 1.0)
    
    b_u = b_u[1:row_idx]
    b_v = b_v[1:row_idx]
    
    A = sparse(I, J, V, row_idx, n_verts)
    
    println("  System size: $(size(A)). Solving Normal Equations...")
    
    AtA = A' * A
    Atb_u = A' * b_u
    Atb_v = A' * b_v
    
    AtA_reg = AtA + 1.0e-6 * spdiagm(0 => ones(n_verts))
    
    F = cholesky(AtA_reg)
    
    u_sol = F \ Atb_u
    v_sol = F \ Atb_v
    
    println("  Range U: [$(minimum(u_sol)), $(maximum(v_sol))]")
    println("  Range V: [$(minimum(v_sol)), $(maximum(v_sol))]")
    
    return u_sol, v_sol
end

end