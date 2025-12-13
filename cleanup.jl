using LinearAlgebra
using SparseArrays
using Printf
using CairoMakie
using DataStructures # You might need to Pkg.add("DataStructures")

# ==============================================================================
# 1. DATA STRUCTURES & TOPOLOGY
# ==============================================================================

struct Point3D
    x::Float64; y::Float64; z::Float64
end

struct MeshTopology
    vertices::Vector{Point3D}
    faces::Vector{Tuple{Int,Int,Int}}
    dual_adj::Vector{Vector{Tuple{Int, Tuple{Int,Int}}}} # face -> [(neighbor, edge_key)]
    face_ref_edges::Vector{Tuple{Int,Int}}
end

mutable struct CrossField
    topology::MeshTopology
    theta::Vector{Float64}
    period_jumps::Dict{Tuple{Int,Int}, Int}
    transport_angles::Dict{Tuple{Int,Int}, Float64}
    constrained_faces::Dict{Int, Float64}
    fixed_edges::Set{Tuple{Int,Int}}
end

# ==============================================================================
# 2. GEOMETRY KERNEL (Same as before)
# ==============================================================================

function read_msh(filename::String)
    lines = readlines(filename)
    vertices = Point3D[]
    faces = Tuple{Int,Int,Int}[]
    tag_map = Dict{Int, Int}()
    i = 1
    while i <= length(lines)
        line = strip(lines[i])
        if line == "\$Nodes"
            i += 1; dims = parse.(Int, split(lines[i])); num_blocks = dims[1]; i += 1
            for _ in 1:num_blocks
                block_head = parse.(Int, split(lines[i])); n_nodes = block_head[4]; i += 1
                tags = Int[]; for _ in 1:n_nodes; push!(tags, parse(Int, lines[i])); i += 1; end
                for k in 1:n_nodes
                    coords = parse.(Float64, split(lines[i]))
                    push!(vertices, Point3D(coords[1], coords[2], coords[3]))
                    tag_map[tags[k]] = length(vertices); i += 1
                end
            end
        elseif line == "\$Elements"
            i += 1; dims = parse.(Int, split(lines[i])); num_blocks = dims[1]; i += 1
            for _ in 1:num_blocks
                block_head = parse.(Int, split(lines[i])); type = block_head[3]; n_elems = block_head[4]; i += 1
                if type == 2
                    for _ in 1:n_elems
                        vals = parse.(Int, split(lines[i]))
                        push!(faces, (tag_map[vals[2]], tag_map[vals[3]], tag_map[vals[4]]))
                        i += 1
                    end
                else; i += n_elems; end
            end
        else; i += 1; end
    end
    return vertices, faces
end

"""
    read_su2(filename)

Parses an SU2 format mesh file (.su2).
Returns vertices and triangle faces.
Note: SU2 uses 0-based indexing, which is converted to 1-based here.
"""
function read_su2(filename::String)
    lines = readlines(filename)
    vertices = Point3D[]
    faces = Tuple{Int,Int,Int}[]
    
    # State tracking
    i = 1
    n_lines = length(lines)
    
    while i <= n_lines
        line = strip(lines[i])
        
        # Skip comments or empty lines
        if isempty(line) || startswith(line, "%")
            i += 1
            continue
        end
        
        # Split by "=" to find keywords
        parts = split(line, "=")
        keyword = strip(parts[1])
        
        if keyword == "NDIME"
            # We assume 2D or 3D is handled by Point3D generic storage
            i += 1
            
        elseif keyword == "NPOIN"
            n_points = parse(Int, strip(parts[2]))
            i += 1
            # Read N lines of vertices
            for _ in 1:n_points
                # Format: x y [z] index (index is usually redundant sequential)
                coords = parse.(Float64, split(lines[i]))
                
                # Handle 2D (x y) or 3D (x y z) cases
                if length(coords) >= 3 && !occursin("index", lines[i]) # Heuristic for 3D
                     push!(vertices, Point3D(coords[1], coords[2], coords[3]))
                else
                     push!(vertices, Point3D(coords[1], coords[2], 0.0))
                end
                i += 1
            end
            
        elseif keyword == "NELEM"
            n_elems = parse(Int, strip(parts[2]))
            i += 1
            # Read N lines of elements
            for _ in 1:n_elems
                # Format: Type v0 v1 v2 ...
                # VTK Type 5 = Triangle, Type 9 = Quad, Type 3 = Line
                vals = parse.(Int, split(lines[i]))
                elem_type = vals[1]
                
                if elem_type == 5 # Triangle
                    # SU2 is 0-based, so add 1
                    n1 = vals[2] + 1
                    n2 = vals[3] + 1
                    n3 = vals[4] + 1
                    push!(faces, (n1, n2, n3))
                end
                # We ignore lines (3) or quads (9) for the pure triangular topology
                i += 1
            end
            
        elseif keyword == "NMARK"
            # Skip boundary markers for now (or implement if you need named boundaries)
            # Usually involves reading "MARKER_TAG= name" and then "MARKER_ELEMS= N"
            n_marks = parse(Int, strip(parts[2]))
            i += 1
            # Fast forward loop to skip markers
            # (A robust parser would handle this, but for simple geometry import we skip)
            # Logic: Read until we find a line that isn't a marker block? 
            # SU2 markers are usually at the end of the file anyway.
            break 
            
        else
            i += 1
        end
    end
    
    return vertices, faces
end

"""
    read_mesh(filename)

Generic mesh reader that dispatches to the correct parser based on extension.
"""
function read_mesh(filename::String)
    if endswith(filename, ".msh")
        return read_msh(filename)
    elseif endswith(filename, ".su2")
        return read_su2(filename)
    else
        error("Unsupported file format: $filename")
    end
end

function build_topology(vertices::Vector{Point3D}, faces::Vector{Tuple{Int,Int,Int}})
    n_faces = length(faces)
    dual_adj = [Tuple{Int, Tuple{Int,Int}}[] for _ in 1:n_faces]
    face_ref_edges = Vector{Tuple{Int,Int}}(undef, n_faces)
    edge_history = Dict{Tuple{Int,Int}, Int}()
    for (i, tri) in enumerate(faces)
        face_ref_edges[i] = (tri[1], tri[2])
        vs = [tri[1], tri[2], tri[3]]
        for k in 1:3
            v1, v2 = vs[k], vs[mod1(k+1, 3)]; key = minmax(v1, v2)
            if haskey(edge_history, key)
                n = edge_history[key]; push!(dual_adj[i], (n, key)); push!(dual_adj[n], (i, key))
            else; edge_history[key] = i; end
        end
    end
    return MeshTopology(vertices, faces, dual_adj, face_ref_edges)
end

function get_edge_vector(topo, v1, v2)
    p1 = topo.vertices[v1]; p2 = topo.vertices[v2]
    return (p2.x - p1.x, p2.y - p1.y, p2.z - p1.z)
end

function compute_transport_angles(topo)
    kappa = Dict{Tuple{Int,Int}, Float64}()
    for i in 1:length(topo.faces)
        for (j, _) in topo.dual_adj[i]
            if i >= j continue end
            vi = topo.faces[i]; vj = topo.faces[j]; shared = intersect(vi, vj)
            u, v = shared[1], shared[2]; ev = get_edge_vector(topo, u, v)
            
            ri = topo.face_ref_edges[i]; rvi = get_edge_vector(topo, ri[1], ri[2])
            ai = atan(ev[2], ev[1]) - atan(rvi[2], rvi[1])
            
            rj = topo.face_ref_edges[j]; rvj = get_edge_vector(topo, rj[1], rj[2])
            aj = atan(ev[2], ev[1]) - atan(rvj[2], rvj[1])
            
            kappa[(i, j)] = mod2pi(aj - ai + π) - π
        end
    end
    return kappa
end

# ==============================================================================
# 3. FAST SOLVER (Local Gauss-Seidel)
# ==============================================================================

function initialize_field(topo, constraints)
    kappa = compute_transport_angles(topo)
    n = length(topo.faces); theta = zeros(n); fixed = Set{Tuple{Int,Int}}()
    p_jumps = Dict{Tuple{Int,Int}, Int}()
    for (f, a) in constraints; theta[f] = a; end
    return CrossField(topo, theta, p_jumps, kappa, constraints, fixed)
end

function build_spanning_tree!(field)
    visited = falses(length(field.topology.faces)); queue = Int[]
    start = isempty(field.constrained_faces) ? 1 : first(keys(field.constrained_faces))
    push!(queue, start); visited[start] = true
    while !isempty(queue)
        u = popfirst!(queue)
        for (v, _) in field.topology.dual_adj[u]
            if !visited[v]
                visited[v] = true; push!(queue, v)
                key = minmax(u, v); push!(field.fixed_edges, key); field.period_jumps[key] = 0
            end
        end
    end
end

# Global solver (still used for initialization)
function solve_global!(field)
    topo = field.topology; n = length(topo.faces)
    I, J, V = Int[], Int[], Float64[]; b = zeros(n)
    
    # Only process FIXED edges for the theta system
    # Free edges do not constrain theta (they absorb error into p)
    for i in 1:n
        for (j, _) in topo.dual_adj[i]
            if i >= j continue end
            edge = (i, j)
            if edge in field.fixed_edges
                k_ij = field.transport_angles[edge]
                p_val = Float64(field.period_jumps[edge])
                # Term: (theta_i - theta_j + shift)^2
                shift = k_ij + p_val * π/2
                # Deriv i: 2(theta_i - theta_j + shift) -> +1 theta_i, -1 theta_j = -shift
                push!(I, i); push!(J, i); push!(V, 1.0)
                push!(I, i); push!(J, j); push!(V, -1.0)
                b[i] += -shift
                
                # Deriv j: -2(theta_i - theta_j + shift) -> -1 theta_i, +1 theta_j = +shift
                push!(I, j); push!(J, j); push!(V, 1.0)
                push!(I, j); push!(J, i); push!(V, -1.0)
                b[j] += shift
            end
        end
    end
    
    penalty = 10000.0
    for (f, ang) in field.constrained_faces
        push!(I, f); push!(J, f); push!(V, penalty); b[f] += penalty * ang
    end
    
    # If system is empty (no fixed edges/constraints), do nothing
    if !isempty(V)
        A = sparse(I, J, V, n, n)
        field.theta = A \ b
    end
end

"""
    local_gauss_seidel!(field, seeds)

Propagates changes starting from `seeds` (face indices) using local updates.
Updates theta until convergence in the local area.
"""
function local_gauss_seidel!(field::CrossField, seeds::Vector{Int})
    queue = Queue{Int}()
    in_queue = Set{Int}()
    
    for s in seeds
        enqueue!(queue, s); push!(in_queue, s)
    end
    
    topo = field.topology
    tol = 1e-7
    max_iter = 100000 # Safety break
    iter = 0
    
    while !isempty(queue) && iter < max_iter
        iter += 1
        u = dequeue!(queue); pop!(in_queue, u)
        
        # Calculate optimal theta for u given neighbors
        # E_local = Sum_neighbors (theta_u - theta_v + shift)^2 + penalty*(theta_u - target)^2
        # dE/dtheta_u = Sum 2(theta_u - theta_v + shift) + 2*penalty*(theta_u - target) = 0
        
        numerator = 0.0
        denominator = 0.0
        
        # 1. Neighbor contributions (Only across FIXED edges)
        for (v, _) in topo.dual_adj[u]
            edge = minmax(u, v)
            if edge in field.fixed_edges
                k_uv = field.transport_angles[edge] # defined for min->max
                p = field.period_jumps[edge]
                
                # Equation is based on theta_min - theta_max + shift = 0
                # shift = k + p*pi/2
                shift = k_uv + p * π/2
                
                denominator += 1.0
                if u < v
                    # term: theta_u - theta_v + shift => theta_u_opt = theta_v - shift
                    numerator += field.theta[v] - shift
                else
                    # term: theta_v - theta_u + shift => -theta_u + ... => theta_u_opt = theta_v + shift
                    numerator += field.theta[v] + shift
                end
            end
        end
        
        # 2. Constraint contribution
        if haskey(field.constrained_faces, u)
            w = 10000.0
            denominator += w
            numerator += w * field.constrained_faces[u]
        end
        
        if denominator > 1e-6
            new_val = numerator / denominator
            diff = abs(new_val - field.theta[u])
            
            if diff > tol
                field.theta[u] = new_val
                # Push neighbors to queue if they are connected via FIXED edges
                for (v, _) in topo.dual_adj[u]
                    edge = minmax(u, v)
                    if edge in field.fixed_edges && !(v in in_queue)
                        enqueue!(queue, v); push!(in_queue, v)
                    end
                end
            end
        end
    end
    @printf("Local Gauss-Seidel converged in %d iterations.\n", iter)
end

"""
Computes the optimal (continuous) p for all free edges based on current theta.
This effectively "solves" the free variables.
"""
function update_free_p_values!(field::CrossField)
    free_p = Tuple{Tuple{Int,Int}, Float64}[]
    
    for i in 1:length(field.topology.faces)
        for (j, _) in field.topology.dual_adj[i]
            if i < j
                edge = (i, j)
                if !(edge in field.fixed_edges)
                    # p_optimal = -(theta_i - theta_j + k_ij) / (pi/2)
                    k = field.transport_angles[edge]
                    val = -(field.theta[i] - field.theta[j] + k) / (π/2)
                    push!(free_p, (edge, val))
                end
            end
        end
    end
    return free_p
end

function solve_miq!(field::CrossField; callback=nothing, debug =false, local_method = false)
    println("--- Starting Fast MIQ Solver ---")
    
    # 1. Init with Global Solve on Tree
    build_spanning_tree!(field)
    solve_global!(field) 
    
    if !isnothing(callback); callback(:start, 0, nothing); end
    
    iter = 0
    while true
        iter += 1
        
        # A. Calculate free p values (lazy update)
        p_candidates = update_free_p_values!(field)
        
        if isempty(p_candidates) break end
        
        # Hook for relax phase (visualize current smooth state)
        if !isnothing(callback); callback(:relax, iter, nothing); end
        
        # B. Greedy Selection
        best_err = Inf; best_idx = -1
        for k in 1:length(p_candidates)
            val = p_candidates[k][2]
            err = abs(val - round(val))
            if err < best_err; best_err = err; best_idx = k; end
        end
        
        # C. Round and Fix
        (edge, p_float) = p_candidates[best_idx]
        p_int = Int(round(p_float))
        debug && @printf("Iter %d: Fixing edge (%d, %d) with p = %d (err=%.4f)\n", iter, edge[1], edge[2], p_int, best_err)
        
        field.period_jumps[edge] = p_int
        push!(field.fixed_edges, edge)
        
        # D. LOCAL GAUSS-SEIDEL UPDATE
        # Only update the two faces adjacent to the new constraint, and propagate
        local_gauss_seidel!(field, [edge[1], edge[2]])
        
        if !isnothing(callback); callback(:round, iter, (edge, p_int)); end
    end
    
    # Final global smooth to be safe
    solve_global!(field)
    if !isnothing(callback); callback(:final, iter, nothing); end
    println("--- Solver Complete ---")
end


using DataStructures # Required for PriorityQueue

"""
    assemble_full_system(field)

Builds the full system A*x = b where x = [theta; p_free].
Returns (A, b, var_map) where var_map maps edge keys to indices in x.
"""
function assemble_full_system(field::CrossField)
    topo = field.topology
    n_faces = length(topo.faces)
    
    free_edges = Tuple{Int,Int}[]
    edge_to_idx = Dict{Tuple{Int,Int}, Int}()
    
    for i in 1:n_faces
        for (j, _) in topo.dual_adj[i]
            if i < j
                edge = (i, j)
                if !(edge in field.fixed_edges)
                    push!(free_edges, edge)
                    edge_to_idx[edge] = n_faces + length(free_edges)
                end
            end
        end
    end
    
    n_vars = n_faces + length(free_edges)
    I = Int[]; J = Int[]; V = Float64[]; b = zeros(n_vars)
    
    for i in 1:n_faces
        # Compute Area of Face i (for weighting)
        # (Simplified: just use 1.0 or heuristic if full geom not avail)
        # Better: Weight by edge length
        
        for (j, _) in topo.dual_adj[i]
            if i >= j continue end
            edge = (i, j)
            k_ij = field.transport_angles[edge]
            idx_i = i; idx_j = j
            idx_p = get(edge_to_idx, edge, 0)
            
            p_val_fixed = (edge in field.fixed_edges) ? Float64(field.period_jumps[edge]) : 0.0
            
            # --- GEOMETRIC WEIGHTING ---
            # Retrieve shared vertices to compute edge length
            verts_i = topo.faces[i]; verts_j = topo.faces[j]
            shared = intersect(verts_i, verts_j)
            u, v = shared[1], shared[2]
            p1 = topo.vertices[u]; p2 = topo.vertices[v]
            len_sq = (p1.x-p2.x)^2 + (p1.y-p2.y)^2 + (p1.z-p2.z)^2
            
            # Weight = 1.0 / Length (penalize twisting on short edges more)
            # This is critical for CFD boundary layers!
            weight = 1.0 / (sqrt(len_sq) + 1e-6)
            # ---------------------------

            indices = [idx_i, idx_j]
            coeffs = [1.0, -1.0]
            
            if idx_p > 0
                push!(indices, idx_p); push!(coeffs, π/2)
            end
            
            rhs_term = -(k_ij + p_val_fixed * π/2)
            
            for r in 1:length(indices)
                row = indices[r]
                b[row] += coeffs[r] * rhs_term * weight # Apply Weight
                for c in 1:length(indices)
                    col = indices[c]
                    push!(I, row); push!(J, col)
                    push!(V, coeffs[r] * coeffs[c] * weight) # Apply Weight
                end
            end
        end
    end
    
    # Increase constraint penalty
    penalty = 1e6 
    for (f, ang) in field.constrained_faces
        push!(I, f); push!(J, f); push!(V, penalty)
        b[f] += penalty * ang
    end
    
    A = sparse(I, J, V, n_vars, n_vars)
    return A, b, edge_to_idx, free_edges
end

"""
    solve_miq_exact!(field; callback=nothing)

Implements the user's Algorithm 1 exactly:
- Unified variable vector x
- Residual-based Local Gauss-Seidel
- Proper dependency queuing
"""
function solve_miq_exact!(field::CrossField; callback=nothing)
    println("--- Starting Exact MIQ Solver (Algorithm 1) ---")
    
    # 1. Initialize Constraints
    build_spanning_tree!(field)
    
    # 2. Build the FULL System ONCE
    # Note: As we fix variables, we won't rebuild A. 
    # We will just strictly enforce the integer value in x and skip updating it in GS.
    A, b, edge_to_idx, free_edges_list = assemble_full_system(field)
    
    n_theta = length(field.topology.faces)
    n_vars = length(b)
    
    # Initial Solve (Continuous Relaxation)
    x = A \ b
    
    # Update field state from x
    field.theta = x[1:n_theta]
    
    # Track which p-variables (indices > n_theta) have been rounded
    # In this logic, 'k' is the number of integer variables
    integer_indices = [edge_to_idx[e] for e in free_edges_list]
    rounded_set = Set{Int}()
    
    iter = 0
    max_k = length(integer_indices)
    
    # Pre-compute column dependencies for A (for fast queue lookups)
    # A is sparse CSC. A.colptr and A.rowval give us the non-zeros.
    # dependency_graph[j] = indices `i` where A[i,j] != 0
    # Because A is symmetric, row `j` also has non-zeros at these indices.
    
    if !isnothing(callback); callback(:start, 0, nothing); end
    
    # === Main Greedy Loop ===
    for step in 1:max_k
        iter += 1
        
        # # 1. Find variable with minimum rounding error
        # best_err = Inf
        # best_var_idx = -1
        
        # for idx in integer_indices
        #     if !(idx in rounded_set)
        #         val = x[idx]
        #         err = abs(val - round(val))
        #         if err < best_err
        #             best_err = err
        #             best_var_idx = idx
        #         end
        #     end
        # end

        # 1. Find variable closest to zero and round it to zero
        best_err = Inf
        best_var_idx = -1
        for idx in integer_indices
            if !(idx in rounded_set)
                val = x[idx]
                err = abs(val)
                if err < best_err
                    best_err = err
                    best_var_idx = idx
                end
            end
        end
        
        if best_var_idx == -1; break; end
        
        # 2. Round and Fix
        rounded_val = round(x[best_var_idx])
        x[best_var_idx] = rounded_val
        push!(rounded_set, best_var_idx)
        
        # Update field for visualization/final result
        # Find which edge this index belongs to
        # (Reverse lookup or just store it. For visualization we update all)
        # We need to lock this value so GS doesn't change it.
        
        # 3. Local Gauss-Seidel Update
        # "Add all variables whose GS update depends on x_i to priority queue"
        # These are rows `m` where A[m, best_var_idx] != 0
        
        # Queue stores unique indices to update
        # Using a Set for uniqueness + a Queue for processing
        Q = Queue{Int}()
        in_Q = Set{Int}()
        
        # Add initial dependencies
        # Iterate non-zeros in column `best_var_idx`
        r_start = A.colptr[best_var_idx]
        r_end = A.colptr[best_var_idx+1] - 1
        
        for k in r_start:r_end
            row_idx = A.rowval[k]
            # Don't add rounded integer variables to queue (they are constants now)
            if !(row_idx in rounded_set) && (row_idx <= n_theta || row_idx in integer_indices)
                if !(row_idx in in_Q)
                    enqueue!(Q, row_idx); push!(in_Q, row_idx)
                end
            end
        end
        
        # GS Parameters
        tau = 1e-4
        gs_iter = 0
        max_gs_iter = 10000 # Safety
        
        while !isempty(Q) && gs_iter < max_gs_iter
            gs_iter += 1
            m = dequeue!(Q); pop!(in_Q, m)
            
            # Compute Residual r_m = b_m - Sum(A_mn * x_n)
            # Efficient sparse row dot product
            # Since A is symmetric CSC, row m is same as col m
            row_start = A.colptr[m]
            row_end = A.colptr[m+1] - 1
            
            Ax_sum = 0.0
            A_mm = 0.0
            
            for k in row_start:row_end
                col_idx = A.rowval[k]
                val = A.nzval[k]
                Ax_sum += val * x[col_idx]
                if col_idx == m
                    A_mm = val
                end
            end
            
            r_m = b[m] - Ax_sum
            
            if abs(r_m) >= tau
                delta = r_m / A_mm
                x[m] += delta
                
                # Push dependencies
                # Rows `n` where A[n, m] != 0
                c_start = A.colptr[m]
                c_end = A.colptr[m+1] - 1
                
                for k in c_start:c_end
                    dep_idx = A.rowval[k]
                    # Only add if not fixed integer
                    is_fixed_int = (dep_idx in rounded_set)
                    if !is_fixed_int
                        if !(dep_idx in in_Q)
                            enqueue!(Q, dep_idx); push!(in_Q, dep_idx)
                        end
                    end
                end
            end
        end
        
        # Visualization Hook
        field.theta = x[1:n_theta]
        # (Optional: Update field.period_jumps from x for viz)
        if !isnothing(callback)
             # Identify edge for viz
             # This is slow O(N), purely for visual debug
             target_edge = (0,0)
             for (e, i) in edge_to_idx; if i == best_var_idx; target_edge = e; break; end; end
             callback(:round, iter, (target_edge, Int(rounded_val)))
        end
    end
    
    # Final cleanup
    # Extract final P values to field struct
    for (edge, idx) in edge_to_idx
        if idx in rounded_set
            field.period_jumps[edge] = Int(round(x[idx]))
        end
    end
    field.theta = x[1:n_theta]
    
    if !isnothing(callback); callback(:final, iter, nothing); end
    println("--- Solver Complete ---")
end

# ==============================================================================
# 5. SINGULARITY EXTRACTION
# ==============================================================================

"""
Helper: Build a lookup of vertex_index -> [face_indices...]
"""
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

"""
Helper: Get the corner angle of a triangle at a specific vertex.
"""
function get_corner_angle(topo::MeshTopology, face_idx::Int, vert_idx::Int)
    tri = topo.faces[face_idx]
    
    # Identify which of the 3 vertices is `vert_idx`
    # Let's say tri is (u, v, w). If vert_idx == u, we want angle at u.
    p1 = topo.vertices[tri[1]]
    p2 = topo.vertices[tri[2]]
    p3 = topo.vertices[tri[3]]
    
    # Map vertex index to A, B, C positions
    # We want angle at A, which is vectors AB and AC
    if tri[1] == vert_idx
        center, left, right = p1, p2, p3
    elseif tri[2] == vert_idx
        center, left, right = p2, p3, p1
    else
        center, left, right = p3, p1, p2
    end
    
    u = (left.x - center.x, left.y - center.y, left.z - center.z)
    v = (right.x - center.x, right.y - center.y, right.z - center.z)
    
    # Angle = acos( dot(u,v) / (|u|*|v|) )
    dot_val = u[1]*v[1] + u[2]*v[2] + u[3]*v[3]
    mag_u = sqrt(u[1]^2 + u[2]^2 + u[3]^2)
    mag_v = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
    
    return acos(clamp(dot_val / (mag_u * mag_v), -1.0, 1.0))
end

"""
    compute_singularities(field)

Computes the index I(v) for every vertex using the star traversal formulas.
Returns a list of (vertex_index, singularity_index).
"""
function compute_singularities(field::CrossField)
    println("\n--- Computing Singularities ---")
    topo = field.topology
    v2f = build_vertex_to_faces(topo)
    singularities = Tuple{Int, Float64}[]
    
    # Iterate over all vertices
    for v_idx in 1:length(topo.vertices)
        incident_faces = v2f[v_idx]
        
        # Skip isolated vertices
        if isempty(incident_faces) continue end
        
        # 1. Sort faces CCW around the vertex to form the "Star"
        # We start with the first face, find the neighbor that shares an edge connected to v, and traverse.
        sorted_faces = Int[]
        
        # Start with arbitrary face
        current_f = incident_faces[1]
        push!(sorted_faces, current_f)
        
        # Walk around
        # (Note: robust handling of boundary vertices requires checking if the loop closes. 
        # This simple looper assumes a closed 1-ring (interior vertex) or walks until it hits a boundary.)
        
        # For simplicity in this snippet, we will assume manifold interior vertices or simple boundaries.
        # We look for the "next" face in the dual adjacency that shares vertex `v_idx`.
        
        found_next = true
        while found_next
            found_next = false
            current_f = sorted_faces[end]
            
            # Look at neighbors of current_f
            for (neighbor, _) in topo.dual_adj[current_f]
                # If we haven't visited it yet AND it is incident to v_idx
                if !(neighbor in sorted_faces) && (neighbor in incident_faces)
                    # Check if they share an edge connected to v_idx
                    # (They must, since they share v_idx and are adjacent faces)
                    push!(sorted_faces, neighbor)
                    found_next = true
                    break
                end
            end
            
            # Stop if we looped back (for closed ring) or hit boundary
            if length(sorted_faces) == length(incident_faces)
                break 
            end
        end

        # If it's a boundary vertex (loop not closed), we usually skip singularity index 
        # or calculate it differently. Let's strictly check for Closed Loop 1-ring.
        # Check if last face is connected to first face
        is_closed_loop = false
        for (n, _) in topo.dual_adj[sorted_faces[end]]
            if n == sorted_faces[1]; is_closed_loop = true; break; end
        end
        
        if !is_closed_loop
            # Boundary vertex: Skip for standard MIQ singularity definition
            continue 
        end
        
        # 2. Compute Quantities over the Star
        angle_sum = 0.0
        kappa_sum = 0.0
        p_sum = 0.0
        
        num_star = length(sorted_faces)
        
        for i in 1:num_star
            f_curr = sorted_faces[i]
            f_next = sorted_faces[mod1(i+1, num_star)]
            
            # A. Angle Defect Component
            angle_sum += get_corner_angle(topo, f_curr, v_idx)
            
            # B. Kappa and P Sums (Edge Traversal)
            # We are crossing the edge between f_curr and f_next
            # We need to ensure we get the sign right: (f_curr -> f_next)
            
            # Retrieve values. Our dicts store canonical keys, so we check orientation.
            # We want Value(curr, next)
            
            # Transport Angle
            if haskey(field.transport_angles, (f_curr, f_next))
                k = field.transport_angles[(f_curr, f_next)]
            else
                k = -field.transport_angles[(f_next, f_curr)]
            end
            kappa_sum += k
            
            # Period Jump
            edge_key = minmax(f_curr, f_next)
            if haskey(field.period_jumps, edge_key)
                # Raw value is for (min, max). 
                raw_p = field.period_jumps[edge_key]
                if f_curr < f_next
                    p_sum += raw_p
                else
                    p_sum += -raw_p
                end
            end
        end
        
        # 3. Calculate Index
        # Ad(v) = 2pi - sum(angles) 
        A_d = 2*π - angle_sum
        
        # I_0(v) = (Ad + Sum Kappa) / 2pi 
        I_0 = (A_d + kappa_sum) / (2*π)
        
        # I(v) = I_0 + Sum p / 4 
        I_v = I_0 + p_sum / 4.0
        
        # 4. Store non-zeros (tolerance for float noise)
        if abs(I_v) > 0.15
            push!(singularities, (v_idx, I_v))
        end
    end
    
    return singularities
end


# ==============================================================================
# 4. VISUALIZATION & MAIN
# ==============================================================================

function get_face_centroid(topo, i)
    tri = topo.faces[i]; p1=topo.vertices[tri[1]]; p2=topo.vertices[tri[2]]; p3=topo.vertices[tri[3]]
    return Point2f((p1.x+p2.x+p3.x)/3, (p1.y+p2.y+p3.y)/3)
end

function get_global_angle(field, i)
    re = field.topology.face_ref_edges[i]; p1=field.topology.vertices[re[1]]; p2=field.topology.vertices[re[2]]
    return atan(p2.y - p1.y, p2.x - p1.x) + field.theta[i]
end

function plot_mesh!(ax, topo)
    for tri in topo.faces
        pts = [topo.vertices[i] for i in tri]; push!(pts, pts[1])
        lines!(ax, [p.x for p in pts], [p.y for p in pts], color=:gray80, linewidth=1)
    end
end

function plot_field!(ax, field; color=:blue, scale::Float64=0.01)
    xs, ys, us, vs = Float64[], Float64[], Float64[], Float64[]
    for i in 1:length(field.topology.faces)
        if haskey(field.constrained_faces, i) continue end
        c = get_face_centroid(field.topology, i); a = get_global_angle(field, i)
        for k in 0:3
            ang = a + k*π/2
            push!(xs, c[1]); push!(ys, c[2]); push!(us, cos(ang)*scale); push!(vs, sin(ang)*scale)
        end
    end
    arrows2d!(ax, xs, ys, us, vs, color=color)
end

function plot_constraints!(ax, field; scale::Float64=0.01)
    for (f, _) in field.constrained_faces
        tri = field.topology.faces[f]
        pts = [Point2f(field.topology.vertices[i].x, field.topology.vertices[i].y) for i in tri]
        poly!(ax, pts, color=(:red, 0.2), strokewidth=2, strokecolor=:red)
    end
    xs, ys, us, vs = Float64[], Float64[], Float64[], Float64[]
    for (f, _) in field.constrained_faces
        c = get_face_centroid(field.topology, f); a = get_global_angle(field, f)
        for k in 0:3
            ang = a + k*π/2
            push!(xs, c[1]); push!(ys, c[2]); push!(us, cos(ang)*scale); push!(vs, sin(ang)*scale)
        end
    end
    arrows2d!(ax, xs, ys, us, vs, color=:red)
end

function run_animated_solver(field, filename)
    fig = Figure(size = (800, 800))
    ax = Axis(fig[1, 1], aspect = DataAspect()); hidedecorations!(ax)
    record(fig, filename; framerate = 5) do io
        function cb(phase, iter, info)
            empty!(ax); plot_mesh!(ax, field.topology); plot_constraints!(ax, field)
            if phase == :start; plot_field!(ax, field, color=:red); recordframe!(io)
            elseif phase == :relax; plot_field!(ax, field, color=:orange); recordframe!(io)
            elseif phase == :round; plot_field!(ax, field, color=:blue); recordframe!(io)
            elseif phase == :final; plot_field!(ax, field, color=:green); for _ in 1:10; recordframe!(io); end
            end
        end
        solve_miq!(field, callback=cb)
    end
    println("Saved $filename")
end


function run_solver(field)
    solve_miq_exact!(field)
end


"""
    plot_singularities!(ax, field; scale=15)
Calculates and plots singularities: Red = Positive (valence < 4), Blue = Negative (valence > 4).
"""
function plot_singularities!(ax, field::CrossField; scale=15)
    # 1. Compute singularities using the function we defined earlier
    sings = compute_singularities(field)
    
    pos_coords = Point2f[]
    neg_coords = Point2f[]
    
    for (v_idx, index_val) in sings
        # Get vertex position
        v_pos = field.topology.vertices[v_idx]
        pt = Point2f(v_pos.x, v_pos.y)
        
        # Sort into buckets
        if index_val > .5
            push!(pos_coords, pt)
        elseif index_val < .5
            push!(neg_coords, pt)
        end
    end
    
    # 2. Plot Positive (Red) - usually Index +1/4 (Valence 3)
    if !isempty(pos_coords)
        scatter!(ax, pos_coords, color=:red, markersize=scale, strokewidth=1, strokecolor=:black)
    end
    
    # 3. Plot Negative (Blue) - usually Index -1/4 (Valence 5)
    if !isempty(neg_coords)
        scatter!(ax, neg_coords, color=:blue, markersize=scale, strokewidth=1, strokecolor=:black)
    end
end

"""
    save_comparison_png(initial_field, smooth_field, filename)
Saves comparison with singularities highlighted on the solution.
"""
function save_comparison_png(initial_field::CrossField, smooth_field::CrossField, filename::String)
    f = Figure(size = (1200, 600))
    
    # Left: Initial State
    ax1 = Axis(f[1, 1], title = "Initial (Unoptimized)", aspect = DataAspect())
    plot_mesh!(ax1, initial_field.topology)
    plot_constraints!(ax1, initial_field)
    plot_field!(ax1, initial_field, color=(:blue, 0.3))
    hidedecorations!(ax1)
    
    # Right: MIQ Solution
    ax2 = Axis(f[1, 2], title = "MIQ Solution", aspect = DataAspect())
    plot_mesh!(ax2, smooth_field.topology)
    plot_constraints!(ax2, smooth_field)
    plot_field!(ax2, smooth_field, color=:blue)
    
    # --- ADDED: Plot Singularities on the solution ---
    plot_singularities!(ax2, smooth_field)
    # -------------------------------------------------
    
    hidedecorations!(ax2)
    
    # Optional: Add a Legend
    Legend(f[1, 3], 
        [MarkerElement(color = :red, marker = :circle), MarkerElement(color = :blue, marker = :circle)],
        ["Pos Index (+1/4)", "Neg Index (-1/4)"],
        "Singularities"
    )
    
    save(filename, f)
    println("Saved comparison with singularities to $filename")
end


"""
    compute_boundary_constraints(topo::MeshTopology)

Identifies all boundary faces and computes the theta value required 
to align the cross field with the boundary edge.
"""
function compute_boundary_constraints(topo::MeshTopology)
    constraints = Dict{Int, Float64}()
    
    # We need to find edges that are NOT in the dual_adj list (i.e., have no neighbor)
    for i in 1:length(topo.faces)
        # Get the 3 edges of the face: (v1,v2), (v2,v3), (v3,v1)
        tri = topo.faces[i]
        face_edges = [
            (tri[1], tri[2]),
            (tri[2], tri[3]),
            (tri[3], tri[1])
        ]
        
        # Check each edge
        for edge in face_edges
            # Canonical key for checking adjacency
            key = minmax(edge[1], edge[2])
            
            # Check if this edge connects to a neighbor
            has_neighbor = false
            for (neighbor, neighbor_edge_key) in topo.dual_adj[i]
                if neighbor_edge_key == key
                    has_neighbor = true
                    break
                end
            end
            
            # If no neighbor, it is a boundary edge
            if !has_neighbor
                # 1. Calculate Global Angle of the Boundary Edge
                # Vector P1 -> P2 (Counter-Clockwise order along boundary)
                p1 = topo.vertices[edge[1]]
                p2 = topo.vertices[edge[2]]
                
                # Tangent vector
                tx = p2.x - p1.x
                ty = p2.y - p1.y
                global_tangent_angle = atan(ty, tx)
                
                # 2. Calculate Reference Edge Angle for this face
                ref_edge = topo.face_ref_edges[i]
                rp1 = topo.vertices[ref_edge[1]]
                rp2 = topo.vertices[ref_edge[2]]
                ref_angle = atan(rp2.y - rp1.y, rp2.x - rp1.x)
                
                # 3. Calculate Constraint
                # We want: Ref_Angle + theta = Global_Tangent_Angle
                # So: theta = Global_Tangent_Angle - Ref_Angle
                theta_constraint = global_tangent_angle - ref_angle
                
                # Normalize to (-pi, pi] for cleanliness
                theta_constraint = mod2pi(theta_constraint + π) - π
                
                constraints[i] = theta_constraint
            end
        end
    end
    
    return constraints
end



function main()
    # filename = "triangulations/disk-radial.msh"
    # filename = "triangulations/disk-radial-fine.msh"
    # filename = "triangulations/hole.msh"
    # filename = "triangulations/300_polygon_sphere_100mm.msh"
    filename = "triangulations/mesh_airfoil_0012.su2"

    if !isfile(filename); println("File not found"); return; end
    
    verts, faces = read_mesh(filename)
    topo = build_topology(verts, faces)
    constraints = Dict(10 => 0.0)
    
    # println("Detecting boundary constraints...")
    # constraints = compute_boundary_constraints(topo)
    # println("Constrained $(length(constraints)) boundary faces.")

    # Run
    field = initialize_field(topo, constraints)
    # run_animated_solver(field, "miq_fast.gif")
    # @time run_solver(field)
    run_solver(field)


    save_comparison_png(initialize_field(topo, constraints), field, "miq_comparison_local.png")
    
    # println("\nTop 5 Angles:")
    # for i in 1:min(5, length(faces))
    #     println("Face $i: $(round(rad2deg(field.theta[i]), digits=2))°")
    # end

    singularities = compute_singularities(field)
    println("Found $(length(singularities)) singularities")
end

main()