module Parameterization

using ..Types
using ..Topology
using ..Analysis
using DataStructures
using LinearAlgebra
using SparseArrays

export propagate_orientations, compute_parameterization_least_squares, solve_mixed_integer_parameterization, SystemIndexing

# --- 1. Orientation Propagation (Section 5) ---

function get_global_vector_local(topo::MeshTopology, f_idx::Int, theta::Float64, r::Int)
    # 1. Get Reference Frame
    tri = topo.faces[f_idx]
    ref_edge = topo.face_ref_edges[f_idx] # (v_start, v_end)
    
    p1 = topo.vertices[ref_edge[1]]
    p2 = topo.vertices[ref_edge[2]]
    
    # Basis X (e1) = Normalized Reference Edge
    e1 = [p2.x - p1.x, p2.y - p1.y, p2.z - p1.z]
    normalize!(e1)
    
    # Normal (n)
    v1 = topo.vertices[tri[1]]
    v2 = topo.vertices[tri[2]]
    v3 = topo.vertices[tri[3]]
    u = [v2.x - v1.x, v2.y - v1.y, v2.z - v1.z]
    v = [v3.x - v1.x, v3.y - v1.y, v3.z - v1.z]
    n = cross(u, v)
    normalize!(n)
    
    # Basis Y (e2) = n x e1
    e2 = cross(n, e1)
    
    # 2. Total Angle = Optimized Theta + Integer Rotation * 90 degrees
    total_angle = theta + (r * π / 2.0)
    
    # 3. Rotate in Tangent Plane
    v_global = cos(total_angle) .* e1 .+ sin(total_angle) .* e2
    return v_global
end

"""
    propagate_orientations(field::CrossField, cut_graph::Set{Tuple{Int,Int}})

Assigns an integer rotation r_i ∈ {0, 1, 2, 3} to each face i.
Stops propagation at edges defined in cut_graph.
"""
function propagate_orientations(field::CrossField, cut_graph::Set{Tuple{Int,Int}}; 
                                verbose::Bool=true, 
                                callback::Union{Function, Nothing}=nothing)
    if verbose
        println("\n--- Propagating Orientations (Section 5) ---")
    end
    
    topo = field.topology
    n_faces = length(topo.faces)
    rotations = fill(-1, n_faces)
    q = Queue{Int}()
    iter = 0
    
    # Iterate to ensure coverage (handle disconnected components if any)
    for root_face in 1:length(topo.faces)
        if rotations[root_face] != -1
            continue
        end
        
        enqueue!(q, root_face)
        rotations[root_face] = 0
        
        while !isempty(q)
            u = dequeue!(q)
            iter += 1
            
            if !isnothing(callback)
                callback(rotations, iter)
            end
            
            r_u = rotations[u]
            
            # Iterate over neighbors
            for (v, edge_key) in topo.dual_adj[u]
                
                # 1. STOP at Cut Graph
                if edge_key in cut_graph
                    continue
                end
                
                # 2. Check if already visited
                if rotations[v] != -1
                    continue
                end
                
                # 3. Retrieve Period Jump and Transport Angle
                edge_key = minmax(u, v)
                kappa = get(field.transport_angles, edge_key, 0.0)
                
                # Effective kappa for u -> v
                kappa_uv = (u < v) ? kappa : -kappa
                
                # Calculate effective period jump
                # theta_v - theta_u ≈ kappa_uv + p_eff * (pi/2)
                diff = field.theta[v] - field.theta[u]
                p_eff = round(Int, (diff - kappa_uv) / (π/2))
                
                # 4. Calculate neighbor rotation
                rotations[v] = mod(r_u - p_eff, 4)
                
                enqueue!(q, v)
            end
        end
    end
    
    return rotations
end


# --- 2. System Indexing & Assembly ---

struct SystemIndexing
    v_map::Matrix{Int}              # [face_idx, 1..3] -> system_var_idx (1..N)
    n_continuous_vars::Int          # N (number of unique uv pairs)
    int_vars::Dict{Tuple{Int,Int}, Tuple{Int,Int}} # Edge -> (idx_j, idx_k)
    n_total_vars::Int               # Total size of system (2*N + 2*n_cuts)
    phys_to_vars::Vector{Vector{Int}} # Debug/Vis: physical vertex -> [sys_indices...]
end

function index_variables_full(topo::MeshTopology, cut_graph::Set{Tuple{Int,Int}})
    n_faces = length(topo.faces)
    v_map = zeros(Int, n_faces, 3)
    phys_to_vars = [Int[] for _ in 1:length(topo.vertices)]
    
    # 1. Continuous Variables (u,v) - Logical Disk
    visited_face = falses(n_faces)
    current_idx = 0
    queue = Queue{Int}()
    
    # Start at arbitrary root
    root = 1
    enqueue!(queue, root)
    visited_face[root] = true
    
    # Init root vertices
    for k in 1:3
        current_idx += 1
        v_map[root, k] = current_idx
        push!(phys_to_vars[topo.faces[root][k]], current_idx)
    end
    
    while !isempty(queue)
        f_curr = dequeue!(queue)
        for (f_neigh, edge_key) in topo.dual_adj[f_curr]
            if visited_face[f_neigh]; continue; end
            
            is_cut = (edge_key in cut_graph)
            visited_face[f_neigh] = true
            
            # Map vertices
            tri_c = topo.faces[f_curr]
            tri_n = topo.faces[f_neigh]
            shared = intersect(tri_c, tri_n)
            
            for k in 1:3
                v = tri_n[k]
                if v in shared
                    # Find index in current face
                    idx_in_c = findfirst(isequal(v), tri_c)
                    sys_id = v_map[f_curr, idx_in_c]
                    
                    if is_cut
                        # Cut Boundary: Create NEW variable
                        current_idx += 1
                        v_map[f_neigh, k] = current_idx
                        push!(phys_to_vars[v], current_idx)
                    else
                        # Inner Edge: Reuse variable
                        v_map[f_neigh, k] = sys_id
                    end
                else
                    # New vertex (not shared)
                    current_idx += 1
                    v_map[f_neigh, k] = current_idx
                    push!(phys_to_vars[v], current_idx)
                end
            end
            enqueue!(queue, f_neigh)
        end
    end
    
    n_continuous = current_idx
    
    # 2. Integer Variables (j, k) for Cut Edges
    int_vars = Dict{Tuple{Int,Int}, Tuple{Int,Int}}()
    sorted_cuts = sort(collect(cut_graph))
    
    for edge in sorted_cuts
        # Index = 2*N_continuous + 1, etc...
        idx_j = 2 * n_continuous + 2 * length(int_vars) + 1
        idx_k = idx_j + 1
        int_vars[edge] = (idx_j, idx_k)
    end
    
    n_total = 2 * n_continuous + 2 * length(int_vars)
    
    return SystemIndexing(v_map, n_continuous, int_vars, n_total, phys_to_vars)
end

function compute_parameterization_least_squares(field::CrossField, cut_graph::Set{Tuple{Int,Int}}; 
                                              verbose=true,
                                              fixed_jumps=nothing,
                                              fixed_singularities=nothing)
    topo = field.topology
    rotations = propagate_orientations(field, cut_graph; verbose=verbose)
    
    # Index variables
    sys = index_variables_full(topo, cut_graph)
    N = sys.n_continuous_vars
    N_total = sys.n_total_vars
    
    if verbose
        println("Assembling Mixed-Integer System:")
        println("  Continuous Vars (u,v nodes): $(2*N)")
        println("  Integer Vars (jump j,k):     $(2*length(sys.int_vars))")
        println("  Total Matrix Size:           $(N_total) x $(N_total)")
    end
    
    # Helper to add to sparse matrix
    I_idx = Int[]; J_idx = Int[]; V_val = Float64[]
    b = zeros(N_total)
    
    function add_entry!(r, c, v)
        push!(I_idx, r); push!(J_idx, c); push!(V_val, v)
    end
    
    # --- 1. SMOOTHNESS ENERGY (Standard Laplacian) ---
    w_smooth = 1.0
    
    # Calculate Avg Edge Length (h) for scaling
    h = 0.0
    cnt = 0
    for tri in topo.faces
        p1, p2 = topo.vertices[tri[1]], topo.vertices[tri[2]]
        h += norm([p1.x-p2.x, p1.y-p2.y, p1.z-p2.z])
        cnt += 1
    end
    h = (cnt > 0) ? h/cnt : 1.0

    for f in 1:length(topo.faces)
        tri = topo.faces[f]
        ids = sys.v_map[f, :]
        
        # Geometry
        p = [topo.vertices[v] for v in tri]
        e1 = [p[2].x-p[1].x, p[2].y-p[1].y, p[2].z-p[1].z]
        e2 = [p[3].x-p[1].x, p[3].y-p[1].y, p[3].z-p[1].z]
        normal = cross(e1, e2)
        area2 = norm(normal)
        if area2 < 1e-12; continue; end
        area = 0.5 * area2
        
        # Gradients (rotated 90 deg in plane / area2)
        v_opp = [
            [p[3].x-p[2].x, p[3].y-p[2].y, p[3].z-p[2].z],
            [p[1].x-p[3].x, p[1].y-p[3].y, p[1].z-p[3].z],
            [p[2].x-p[1].x, p[2].y-p[1].y, p[2].z-p[1].z]
        ]
        grad_N = [cross(normal, vec) ./ area2 for vec in v_opp]
        
        # Target Vectors
        r = rotations[f]
        u_global = get_global_vector_local(topo, f, field.theta[f], r) .* h
        v_global = get_global_vector_local(topo, f, field.theta[f] + π/2, r) .* h
        
        # Fill Matrix
        # Energy = Area * || sum(u_i gradN_i) - u_target ||^2 + ...
        for i in 1:3
            row_u = ids[i]
            row_v = ids[i] + N
            
            # RHS: b += 2 * Area * (target . gradN_i) * weight
            b[row_u] += 2 * area * dot(u_global, grad_N[i]) * w_smooth
            b[row_v] += 2 * area * dot(v_global, grad_N[i]) * w_smooth
            
            for j in 1:3
                col_u = ids[j]
                col_v = ids[j] + N
                
                val = 2 * area * dot(grad_N[i], grad_N[j]) * w_smooth
                
                add_entry!(row_u, col_u, val)
                add_entry!(row_v, col_v, val)
            end
        end
    end
    
    # --- 2. CUT COMPATIBILITY CONSTRAINTS ---
    # Enforce via Penalty Method
    w_cut = 1e5 * w_smooth 
    
    for (edge, (idx_j, idx_k)) in sys.int_vars
        faces_sharing = find_faces_for_edge(topo, edge)
        if length(faces_sharing) != 2; continue; end
        f_A, f_B = faces_sharing[1], faces_sharing[2]
        
        # Compare Global Vectors to determine integer rotation
        vu_A = get_global_vector_local(topo, f_A, field.theta[f_A], rotations[f_A])
        vv_A = get_global_vector_local(topo, f_A, field.theta[f_A]+π/2, rotations[f_A])
        vu_B = get_global_vector_local(topo, f_B, field.theta[f_B], rotations[f_B])
        
        c_u = dot(vu_B, vu_A) # cos(theta)
        s_u = dot(vu_B, vv_A) # sin(theta)
        
        # Snap to integer
        if abs(c_u) > abs(s_u)
            c_int = sign(c_u); s_int = 0.0
        else
            c_int = 0.0; s_int = sign(s_u)
        end
        
        # Vertices on the edge
        v1, v2 = edge
        
        # Apply constraint to both vertices on the cut edge
        for v_idx in [v1, v2]
            # Find system indices in face A and B
            i_A = findfirst(isequal(v_idx), topo.faces[f_A])
            i_B = findfirst(isequal(v_idx), topo.faces[f_B])
            
            ua_idx = sys.v_map[f_A, i_A]; va_idx = ua_idx + N
            ub_idx = sys.v_map[f_B, i_B]; vb_idx = ub_idx + N
            j_idx = idx_j
            k_idx = idx_k
            
            # U Constraint: u_B - c*u_A - s*v_A - j = 0
            coeffs_u = [(ub_idx, 1.0), (ua_idx, -c_int), (va_idx, -s_int), (j_idx, -1.0)]
            add_squared_residual!(I_idx, J_idx, V_val, coeffs_u, w_cut)
            
            # V Constraint: v_B + s*u_A - c*v_A - k = 0
            coeffs_v = [(vb_idx, 1.0), (ua_idx, s_int), (va_idx, -c_int), (k_idx, -1.0)]
            add_squared_residual!(I_idx, J_idx, V_val, coeffs_v, w_cut)
        end
    end
    
    # 3. Regularization (Fix one vertex to remove translation freedom)
    w_fix = 1e-3 * w_smooth
    add_entry!(1, 1, w_fix)
    add_entry!(1+N, 1+N, w_fix)
    
    # --- 4. FIXED JUMPS (Optional) ---
    if fixed_jumps !== nothing
        w_int_fix = 1e9 * w_smooth
        for (edge, (val_j, val_k)) in fixed_jumps
            if haskey(sys.int_vars, edge)
                idx_j, idx_k = sys.int_vars[edge]
                
                # Add penalty: w * (j - val_j)^2
                add_entry!(idx_j, idx_j, w_int_fix)
                b[idx_j] += w_int_fix * val_j
                
                add_entry!(idx_k, idx_k, w_int_fix)
                b[idx_k] += w_int_fix * val_k
            end
        end
    end
    
    # --- 5. FIXED SINGULARITIES (Optional) ---
    if fixed_singularities !== nothing
        w_sing_fix = 1e9 * w_smooth
        for (v_idx, (val_u, val_v)) in fixed_singularities
            if !isempty(sys.phys_to_vars[v_idx])
                # Use the first instance of the variable for this vertex
                sys_idx = sys.phys_to_vars[v_idx][1]
                idx_u = sys_idx
                idx_v = sys_idx + N
                
                # Add penalty: w * (u - val_u)^2
                add_entry!(idx_u, idx_u, w_sing_fix)
                b[idx_u] += w_sing_fix * val_u
                
                # Add penalty: w * (v - val_v)^2
                add_entry!(idx_v, idx_v, w_sing_fix)
                b[idx_v] += w_sing_fix * val_v
            end
        end
    end

    # Build Sparse Matrix
    M = sparse(I_idx, J_idx, V_val, N_total, N_total)
    
    return M, b, sys, rotations
end

# --- Helpers ---

function find_faces_for_edge(topo, edge)
    faces = Int[]
    v1, v2 = edge
    key = minmax(v1, v2)
    
    for f in 1:length(topo.faces)
        for (_, k) in topo.dual_adj[f]
            if k == key
                push!(faces, f)
            end
        end
        if length(faces) == 2; break; end
    end
    return faces
end

function add_squared_residual!(I, J, V, coeffs, w)
    # Adds w * (sum c_i x_i)^2 to matrix
    for (idx_i, val_i) in coeffs
        for (idx_j, val_j) in coeffs
            push!(I, idx_i)
            push!(J, idx_j)
            push!(V, w * val_i * val_j)
        end
    end
end



"""
    solve_mixed_integer_parameterization(field, cuts; verbose=true)

Fully solves the Global Parameterization problem.
1. Computes the continuous (relaxed) solution.
2. Iteratively snaps integer variables (singularities & cut jumps) to the nearest integers.
3. Returns the final (u, v) coordinates mapped to physical vertices.

Returns:
    u_phys::Vector{Float64}, v_phys::Vector{Float64}
"""
function solve_mixed_integer_parameterization(field::CrossField, cut_graph::Set{Tuple{Int,Int}}; verbose=true)
    topo = field.topology
    
    # 1. Identify Integer Variables
    sings = compute_singularities(field)
    sing_indices = [s[1] for s in sings] # Just vertex indices
    
    if verbose
        println("\n=== Starting Mixed-Integer Parameterization ===")
        println("  Singularities to snap: $(length(sing_indices))")
    end

    # 2. Initial Relaxed Solve
    if verbose; println("  [1/3] Solving continuous system..."); end
    M, b, sys, _ = compute_parameterization_least_squares(field, cut_graph; verbose=false)
    x = M \ b
    
    # 3. Greedy Rounding Loop
    fixed_jumps = Dict{Tuple{Int,Int}, Tuple{Int,Int}}()
    fixed_sings = Dict{Int, Tuple{Int,Int}}()
    
    # Total integers to fix = #Jumps + #Singularities
    jump_keys = collect(keys(sys.int_vars))
    total_to_fix = length(jump_keys) + length(sing_indices)
    
    if verbose; println("  [2/3] Iterative Rounding ($total_to_fix variables)..."); end
    
    while length(fixed_jumps) + length(fixed_sings) < total_to_fix
        best_err = Inf
        best_type = :none
        best_key = nothing
        best_val = (0, 0)
        
        # A. Check Jumps
        for edge in jump_keys
            if haskey(fixed_jumps, edge); continue; end
            
            idx_j, idx_k = sys.int_vars[edge]
            v_j, v_k = x[idx_j], x[idx_k]
            r_j, r_k = round(Int, v_j), round(Int, v_k)
            
            err = max(abs(v_j - r_j), abs(v_k - r_k))
            if err < best_err
                best_err = err; best_type = :jump; best_key = edge; best_val = (r_j, r_k)
            end
        end
        
        # B. Check Singularities
        for v_idx in sing_indices
            if haskey(fixed_sings, v_idx); continue; end
            
            # Use first variable instance of this vertex
            sys_idx = sys.phys_to_vars[v_idx][1]
            v_u, v_v = x[sys_idx], x[sys_idx + sys.n_continuous_vars]
            r_u, r_v = round(Int, v_u), round(Int, v_v)
            
            err = max(abs(v_u - r_u), abs(v_v - r_v))
            if err < best_err
                best_err = err; best_type = :sing; best_key = v_idx; best_val = (r_u, r_v)
            end
        end
        
        # C. Fix the best candidate
        if best_type == :jump
            fixed_jumps[best_key] = best_val
        elseif best_type == :sing
            fixed_sings[best_key] = best_val
        end
        
        # D. Re-solve (Update system with new constraints)
        # Note: We re-assemble M/b. For large meshes, updating M in place is faster, 
        # but re-assembly is safer and easier to implement.
        M_new, b_new, _, _ = compute_parameterization_least_squares(field, cut_graph; 
                                    fixed_jumps=fixed_jumps, 
                                    fixed_singularities=fixed_sings,
                                    verbose=false)
        x = M_new \ b_new
        
        if verbose && (length(fixed_jumps) + length(fixed_sings)) % 5 == 0
             print(".") # Progress dot
        end
    end
    if verbose; println("\n  [3/3] Final solve complete."); end

    # 4. Extract Result (Map back to Physical Vertices)
    # Note: Vertices on cuts have 2 UVs. For simple visualization, we pick the first one.
    u_phys = zeros(length(topo.vertices))
    v_phys = zeros(length(topo.vertices))
    
    for i in 1:length(topo.vertices)
        if isempty(sys.phys_to_vars[i]); continue; end
        
        sys_idx = sys.phys_to_vars[i][1] # Pick primary instance
        u_phys[i] = x[sys_idx]
        v_phys[i] = x[sys_idx + sys.n_continuous_vars]
    end
    
    return u_phys, v_phys
end

end