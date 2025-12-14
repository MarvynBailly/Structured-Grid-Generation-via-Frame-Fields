using LinearAlgebra
using SparseArrays
using Printf
using CairoMakie
using DataStructures 

# ==============================================================================
# 1. DATA STRUCTURES & TOPOLOGY
# ==============================================================================

struct Point3D
    x::Float64; y::Float64; z::Float64
end

struct MeshTopology
    vertices::Vector{Point3D}
    faces::Vector{Tuple{Int,Int,Int}}
    dual_adj::Vector{Vector{Tuple{Int, Tuple{Int,Int}}}} 
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
# 2. FILE I/O
# ==============================================================================

function read_msh(filename::String)
    lines = readlines(filename)
    vertices = Point3D[]; faces = Tuple{Int,Int,Int}[]; tag_map = Dict{Int, Int}()
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

function read_su2(filename::String)
    lines = readlines(filename)
    vertices = Point3D[]; faces = Tuple{Int,Int,Int}[]
    i = 1
    while i <= length(lines)
        line = strip(lines[i])
        if isempty(line) || startswith(line, "%"); i += 1; continue; end
        parts = split(line, "="); keyword = strip(parts[1])
        if keyword == "NPOIN"
            n_pts = parse(Int, strip(parts[2])); i += 1
            for _ in 1:n_pts
                c = parse.(Float64, split(lines[i]))
                push!(vertices, length(c) >= 3 ? Point3D(c[1], c[2], c[3]) : Point3D(c[1], c[2], 0.0))
                i += 1
            end
        elseif keyword == "NELEM"
            n_el = parse(Int, strip(parts[2])); i += 1
            for _ in 1:n_el
                vals = parse.(Int, split(lines[i]))
                if vals[1] == 5 # Triangle
                    push!(faces, (vals[2]+1, vals[3]+1, vals[4]+1))
                end
                i += 1
            end
        elseif keyword == "NMARK"; break # Stop at markers
        else; i += 1; end
    end
    return vertices, faces
end

function read_mesh(filename::String)
    endswith(filename, ".su2") ? read_su2(filename) : read_msh(filename)
end

# ==============================================================================
# 3. TOPOLOGY & GEOMETRY
# ==============================================================================

function build_topology(vertices, faces)
    n = length(faces)
    dual_adj = [Tuple{Int, Tuple{Int,Int}}[] for _ in 1:n]
    face_ref = Vector{Tuple{Int,Int}}(undef, n)
    edge_hist = Dict{Tuple{Int,Int}, Int}()
    for (i, tri) in enumerate(faces)
        face_ref[i] = (tri[1], tri[2])
        vs = [tri[1], tri[2], tri[3]]
        for k in 1:3
            v1, v2 = vs[k], vs[mod1(k+1, 3)]; key = minmax(v1, v2)
            if haskey(edge_hist, key)
                n = edge_hist[key]; push!(dual_adj[i], (n, key)); push!(dual_adj[n], (i, key))
            else; edge_hist[key] = i; end
        end
    end
    return MeshTopology(vertices, faces, dual_adj, face_ref)
end

function get_edge_vector(topo, v1, v2)
    p1 = topo.vertices[v1]; p2 = topo.vertices[v2]
    return (p2.x - p1.x, p2.y - p1.y, p2.z - p1.z)
end

function compute_transport_angles(topo)
    kappa = Dict{Tuple{Int,Int}, Float64}()
    println("\n--- Computing Transport Angles (kappa) ---")
    for i in 1:length(topo.faces)
        for (j, _) in topo.dual_adj[i]
            if i >= j continue end
            vi = topo.faces[i]; vj = topo.faces[j]; shared = intersect(vi, vj)
            ev = get_edge_vector(topo, shared[1], shared[2])
            
            ri = topo.face_ref_edges[i]; rvi = get_edge_vector(topo, ri[1], ri[2])
            ai = atan(ev[2], ev[1]) - atan(rvi[2], rvi[1])
            
            rj = topo.face_ref_edges[j]; rvj = get_edge_vector(topo, rj[1], rj[2])
            aj = atan(ev[2], ev[1]) - atan(rvj[2], rvj[1])
            
            kappa[(i, j)] = mod2pi(aj - ai + π) - π
            # if i <= 5 || (i <= 20 && abs(kappa[(i,j)]) > 0.1)
            #     @printf("  Edge (%d,%d): kappa = %.4f rad (%.2f deg)\n", i, j, kappa[(i,j)], rad2deg(kappa[(i,j)]))
            # end
        end
    end
    println("Computed $(length(kappa)) transport angles.")
    return kappa
end

function compute_boundary_constraints(topo)
    constraints = Dict{Int, Float64}()
    
    for i in 1:length(topo.faces)
        tri = topo.faces[i]
        # Define the three potential edges of the triangle
        edges = [(tri[1], tri[2]), (tri[2], tri[3]), (tri[3], tri[1])]
        
        for (k, edge) in enumerate(edges)
            key = minmax(edge[1], edge[2])
            has_neighbor = false
            for (_, n_key) in topo.dual_adj[i]
                if n_key == key; has_neighbor = true; break; end
            end
            
            if !has_neighbor
                # --- Boundary Edge Found ---
                
                # 1. Identify the Third Vertex (Internal Point)
                # The triangle has vertices tri[1], tri[2], tri[3].
                # If the edge uses two of them, the third one is the remaining index.
                # In your edges array:
                # k=1 uses (1,2) -> third is 3
                # k=2 uses (2,3) -> third is 1
                # k=3 uses (3,1) -> third is 2
                third_idx = if k == 1; tri[3]; elseif k == 2; tri[1]; else; tri[2]; end
                p3 = topo.vertices[third_idx]
                
                # 2. Get Boundary Edge Coordinates
                p1 = topo.vertices[edge[1]]
                p2 = topo.vertices[edge[2]]
                
                # 3. Compute Tangent Components
                dx = p2.x - p1.x
                dy = p2.y - p1.y
                global_tan = atan(dy, dx)
                
                # 4. Apply Direction Logic based on Third Point Y
                # If third point is BELOW y=0 (Bottom Surface), force Right (+X)
                # If third point is ABOVE y=0 (Top Surface), force Left (-X)
                
                if p3.y < 0.0
                    # Bottom Surface: We want +X direction
                    if dx < 0.0
                        global_tan += π
                    end
                elseif p3.y > 0.0
                    # Top Surface: We want -X direction (assuming specific winding desire)
                    # Based on your previous snippet: "if cy > 0.01 ... if dx > 0.0 global_tan += π"
                    # This implies you want the top surface to flow Left (-X).
                    if dx > 0.0
                        global_tan += π
                    end
                end
                
                # Note: If p3.y is exactly 0.0, we do nothing (preserve original topology)

                # 5. Compute Constraint relative to Face Reference Frame
                re = topo.face_ref_edges[i]
                rp1 = topo.vertices[re[1]]
                rp2 = topo.vertices[re[2]]
                ref_ang = atan(rp2.y - rp1.y, rp2.x - rp1.x)
                
                val = mod2pi(global_tan - ref_ang + π) - π
                constraints[i] = val
            end
        end
    end
    return constraints
end
# ==============================================================================
# 4. UNIFIED GREEDY SOLVER
# ==============================================================================

function initialize_field(topo, constraints)
    kappa = compute_transport_angles(topo)
    n = length(topo.faces); theta = zeros(n); fixed = Set{Tuple{Int,Int}}()
    p_jumps = Dict{Tuple{Int,Int}, Int}()
    for (f, a) in constraints; theta[f] = a; end
    return CrossField(topo, theta, p_jumps, kappa, constraints, fixed)
end

function build_spanning_tree!(field; verbose=false)
    visited = falses(length(field.topology.faces))
    queue = Int[]
    
    if verbose
        println("\n--- Building Spanning Forest ---")
    end
    
    # Initialize queue with ALL constrained faces (creates a forest, not just a tree)
    n_roots = 0
    if !isempty(field.constrained_faces)
        for f in keys(field.constrained_faces)
            push!(queue, f)
            visited[f] = true
            n_roots += 1
        end
        if verbose
            println("Starting from $n_roots constrained faces (roots)")
        end
    else
        # Fallback: if no constraints, start from face 1
        push!(queue, 1)
        visited[1] = true
        n_roots = 1
        if verbose
            println("No constrained faces - starting from face 1")
        end
    end
    
    while !isempty(queue)
        u = popfirst!(queue)
        for (v, _) in field.topology.dual_adj[u]
            if !visited[v]
                visited[v] = true
                push!(queue, v)
                key = minmax(u, v)
                push!(field.fixed_edges, key)
                field.period_jumps[key] = 0
            end
        end
    end
end

"""
Builds the unified system (theta and p) with GEOMETRIC WEIGHTING.
"""
function assemble_system_weighted(field)
    topo = field.topology
    n_faces = length(topo.faces)
    
    # 1. Map Variables
    theta_map = Dict{Int,Int}()
    p_map = Dict{Tuple{Int,Int}, Int}()
    curr_idx = 0
    
    # Map Theta (only unconstrained faces are variables)
    for f in 1:n_faces
        if !haskey(field.constrained_faces, f)
            curr_idx += 1; theta_map[f] = curr_idx
        end
    end
    
    # Map P (only free edges)
    for i in 1:n_faces
        for (j, _) in topo.dual_adj[i]
            if i < j
                edge = (i, j)
                # Skip if fixed or if between two constrained faces (optimization)
                is_constrained_edge = haskey(field.constrained_faces, i) && haskey(field.constrained_faces, j)
                if !(edge in field.fixed_edges) && !is_constrained_edge
                    curr_idx += 1; p_map[edge] = curr_idx
                end
            end
        end
    end
    
    n_vars = curr_idx
    I = Int[]; J = Int[]; V = Float64[]; b = zeros(n_vars)
    
    function add!(r, c, val)
        push!(I, r); push!(J, c); push!(V, val)
    end
    
    # 2. Build Energy
    for i in 1:n_faces
        for (j, _) in topo.dual_adj[i]
            if i >= j continue end
            edge = (i, j)
            kappa = field.transport_angles[edge]
            
            idx_i = get(theta_map, i, 0)
            idx_j = get(theta_map, j, 0)
            idx_p = get(p_map, edge, 0)
            
            # Constraint Values
            val_i = get(field.constrained_faces, i, 0.0)
            val_j = get(field.constrained_faces, j, 0.0)
            
            # P Value logic
            val_p = 0.0
            if idx_p == 0
                if edge in field.fixed_edges
                    val_p = Float64(field.period_jumps[edge])
                elseif haskey(field.constrained_faces, i) && haskey(field.constrained_faces, j)
                    # Round explicitly between two constraints
                    term = (2.0/π) * (val_j - val_i - kappa)
                    val_p = round(term)
                    # We store this back to field for consistency
                    field.period_jumps[edge] = Int(val_p)
                end
            end
            
            # --- WEIGHTING FIX ---
            # Compute edge length
            vi = topo.faces[i]; vj = topo.faces[j]; shared = intersect(vi, vj)
            p1 = topo.vertices[shared[1]]; p2 = topo.vertices[shared[2]]
            len_sq = (p1.x-p2.x)^2 + (p1.y-p2.y)^2 + (p1.z-p2.z)^2
            weight = 1.0 #/ (sqrt(len_sq) + 1e-6)
            # ---------------------
            
            # RHS
            # Equation: 1*th_i - 1*th_j + pi/2*p = -kappa
            rhs_const = -kappa
            if idx_i == 0; rhs_const -= val_i; end
            if idx_j == 0; rhs_const += val_j; end
            if idx_p == 0; rhs_const -= (π/2) * val_p; end
            
            # Add derivatives with Weight scaling
            # Row i
            if idx_i > 0
                add!(idx_i, idx_i, 2.0 * weight)
                if idx_j > 0; add!(idx_i, idx_j, -2.0 * weight); end
                if idx_p > 0; add!(idx_i, idx_p, π * weight); end
                b[idx_i] += 2.0 * rhs_const * weight
            end
            
            # Row j
            if idx_j > 0
                if idx_i > 0; add!(idx_j, idx_i, -2.0 * weight); end
                add!(idx_j, idx_j, 2.0 * weight)
                if idx_p > 0; add!(idx_j, idx_p, -π * weight); end
                b[idx_j] += -2.0 * rhs_const * weight
            end
            
            # Row p
            if idx_p > 0
                if idx_i > 0; add!(idx_p, idx_i, π * weight); end
                if idx_j > 0; add!(idx_p, idx_j, -π * weight); end
                add!(idx_p, idx_p, (π^2/2.0) * weight)
                b[idx_p] += π * rhs_const * weight
            end
        end
    end
    
    A = sparse(I, J, V, n_vars, n_vars)
    return A, b, theta_map, p_map
end

function local_gauss_seidel_queue!(A, b, x, fixed_mask, start_node, diag_A)
    q = Queue{Int}(); in_q = Set{Int}()
    rows = rowvals(A); vals = nonzeros(A)
    
    # Add neighbors of start node
    for idx in nzrange(A, start_node)
        n = rows[idx]
        if !fixed_mask[n] && !(n in in_q)
            enqueue!(q, n); push!(in_q, n)
        end
    end
    
    iter = 0; max_iter = 50000; tol = 1e-5
    while !isempty(q) && iter < max_iter
        iter += 1
        k = dequeue!(q); delete!(in_q, k)
        
        # Residual
        Ax_k = 0.0
        for idx in nzrange(A, k)
            Ax_k += vals[idx] * x[rows[idx]]
        end
        r_k = b[k] - Ax_k
        
        if abs(r_k) > tol
            x[k] += r_k / diag_A[k]
            # Propagate
            for idx in nzrange(A, k)
                n = rows[idx]
                if !fixed_mask[n] && !(n in in_q)
                    enqueue!(q, n); push!(in_q, n)
                end
            end
        end
    end
    return iter
end

function solve_greedy!(field; verbose=false, animate=false, frame_dir="animation_frames", frame_interval=10)
    if verbose
        println("--- Starting Unified MIQ Solver ---")
    end
    build_spanning_tree!(field; verbose=verbose)
    
    # Create animation directory if needed
    if animate
        mkpath(frame_dir)
        println("Saving animation frames to: $frame_dir")
    end
    
    # 1. Build System
    if verbose
        println("Assembling system...")
    end
    A, b, theta_map, p_map = assemble_system_weighted(field)
    n_vars = length(b)
    
    # 2. Initial Solve
    if verbose
        println("Initial solve...")
    end
    x = A \ b
    
    # Extract initial state to field
    for (f, idx) in theta_map
        field.theta[f] = x[idx]
    end
    for (edge, idx) in p_map
        field.period_jumps[edge] = Int(round(x[idx]))
    end
    
    # Save initial frame
    if animate
        frame_path = joinpath(frame_dir, "frame_0000.png")
        plot_frame(field, frame_path, "Initial State (Continuous Solution)")
        println("  Saved initial frame")
    end
    
    # 3. Greedy Rounding
    int_indices = collect(values(p_map))
    fixed_mask = falses(n_vars)
    diag_A = [A[i,i] for i in 1:n_vars]
    
    num_p = length(int_indices)
    if verbose
        println("Greedy rounding $num_p integer variables...")
    end
    
    count = 0
    frame_count = 1
    while true
        best_idx = -1; min_err = Inf; best_val = 0.0
        
        for idx in int_indices
            if !fixed_mask[idx]
                v = x[idx]; r = round(v); err = abs(v-r)
                if err < min_err
                    min_err = err; best_idx = idx; best_val = r
                end
            end
        end
        
        if best_idx == -1; break; end
        
        # Round
        x[best_idx] = best_val
        fixed_mask[best_idx] = true
        count += 1
        
        # Relax
        local_gauss_seidel_queue!(A, b, x, fixed_mask, best_idx, diag_A)
        
        # Save animation frame
        if animate && (count % frame_interval == 0 || count == num_p)
            # Extract current state to field
            for (f, idx) in theta_map
                field.theta[f] = x[idx]
            end
            for (edge, idx) in p_map
                field.period_jumps[edge] = Int(round(x[idx]))
            end
            
            frame_path = joinpath(frame_dir, @sprintf("frame_%04d.png", frame_count))
            progress_pct = round(100 * count / num_p, digits=1)
            title = @sprintf("Greedy Rounding: %d/%d (%.1f%%) - Error: %.4f", count, num_p, progress_pct, min_err)
            plot_frame(field, frame_path, title)
            
            if verbose
                println("  Saved frame $frame_count (iteration $count)")
            end
            frame_count += 1
        end
        
        if verbose && count % 50 == 0
            @printf("  Rounded %d/%d (err: %.4f)\n", count, num_p, min_err)
        end
    end
    
    # 4. Extract Final Results
    for (f, idx) in theta_map
        field.theta[f] = x[idx]
    end
    for (edge, idx) in p_map
        field.period_jumps[edge] = Int(round(x[idx]))
    end
    
    if verbose
        println("--- Solver Complete ---")
        if animate
            println("Animation: Saved $frame_count frames total")
        end
    end
end

# ==============================================================================
# 5. SINGULARITIES & VIZ
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

function compute_singularities(field::CrossField; verbose=false)
    if verbose
        println("\n--- Computing Singularities (Normalized) ---")
    end
    topo = field.topology
    v2f = build_vertex_to_faces(topo)
    singularities = Tuple{Int, Float64}[]
    
    # Statistics
    all_indices = Float64[]
    all_angle_defects = Float64[]
    all_k_sums = Float64[]
    all_p_sums = Float64[]
    boundary_vertices_skipped = 0
    degenerate_vertices_skipped = 0
    
    # Thresholds
    MIN_ANGLE_THRESHOLD = deg2rad(5.0)  # Skip vertices with angles < 5 degrees
    MIN_ANGLE_SUM_THRESHOLD = deg2rad(120.0)  # Skip if total angle sum < 120 degrees
    
    for v_idx in 1:length(topo.vertices)
        faces = v2f[v_idx]
        if isempty(faces) continue end
        
        # Star Traversal (same as before)
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
        
        # Compute Index
        ang_sum = 0.0; k_sum = 0.0; p_sum = 0.0
        n_star = length(star)
        min_angle = Inf
        
        # println("\n=== Vertex $v_idx (star size: $n_star) ===")
        
        for i in 1:n_star
            curr = star[i]; next = star[mod1(i+1, n_star)]
            
            # Angle Defect calculation (same as before)
            tri = topo.faces[curr]; p = topo.vertices; c = p[v_idx]
            if tri[1]==v_idx; l=p[tri[2]]; r=p[tri[3]]
            elseif tri[2]==v_idx; l=p[tri[3]]; r=p[tri[1]]
            else; l=p[tri[1]]; r=p[tri[2]]; end
            u = (l.x-c.x, l.y-c.y, l.z-c.z); v = (r.x-c.x, r.y-c.y, r.z-c.z)
            angle_at_vertex = acos(clamp((u[1]*v[1]+u[2]*v[2]+u[3]*v[3])/(norm(u)*norm(v)), -1, 1))
            ang_sum += angle_at_vertex
            min_angle = min(min_angle, angle_at_vertex)
            
            # @printf("  Face %d: angle = %.4f rad (%.2f deg)\n", curr, angle_at_vertex, rad2deg(angle_at_vertex))
            
            # Transport
            k = haskey(field.transport_angles, (curr, next)) ? field.transport_angles[(curr,next)] : -field.transport_angles[(next,curr)]
            k_sum += k
            # @printf("  Edge (%d->%d): kappa = %.4f rad (%.2f deg)\n", curr, next, k, rad2deg(k))
            
            # Period Jump (NORMALIZED)
            edge = minmax(curr, next)
            sign = (curr < next) ? 1.0 : -1.0
            
            if haskey(field.period_jumps, edge)
                raw_p = field.period_jumps[edge]
                
                # For cross fields (4-RoSy), period jumps represent multiples of π/2
                # No normalization needed - the raw value is what we need
                # The division by 4 in the index formula handles the periodicity
                # @printf("  Edge %s: raw_p = %d (sign = %.1f, contribution = %.2f)\n", 
                #         edge, raw_p, sign, raw_p * sign)
                p_sum += raw_p * sign
            else
                println("  Edge $edge: no period jump")
            end
        end
        
        angle_defect = 0#2π - ang_sum
        # For cross fields: index = (angle_defect)/(2π) + (kappa_sum)/(2π) + (period_sum)/4
        # The period_sum/4 accounts for the 4-fold symmetry (π/2 periodicity)
        I = (angle_defect + k_sum)/(2π) + p_sum/4.0
        
        # Round to nearest 1/4 for cross fields
        # I_rounded = round(I * 4) / 4
        # @printf("  INDEX rounded to nearest 1/4: %.4f (%.2f/4)\n", I_rounded, I_rounded * 4)
        
        # Check for degenerate geometry
        is_degenerate = (min_angle < MIN_ANGLE_THRESHOLD) #|| (ang_sum < MIN_ANGLE_SUM_THRESHOLD)
        
        if is_degenerate
            # println("  ⚠️  DEGENERATE VERTEX SKIPPED (nearly-zero angles indicate bad mesh quality)")
            # degenerate_vertices_skipped += 1
            #             @printf("  Angle sum: %.4f rad (%.2f deg)\n", ang_sum, rad2deg(ang_sum))
            # @printf("  Angle defect: %.4f rad (%.2f deg)\n", angle_defect, rad2deg(angle_defect))
            # @printf("  Kappa sum: %.4f rad (%.2f deg)\n", k_sum, rad2deg(k_sum))
            # @printf("  Period sum (raw): %.4f\n", p_sum)
            # @printf("  Min angle: %.4f rad (%.2f deg)\n", min_angle, rad2deg(min_angle))
            # @printf("  INDEX = (%.4f + %.4f)/(2π) + %.4f/4 = %.4f\n", angle_defect, k_sum, p_sum, I)
            # @printf("  INDEX breakdown: geom=(%.4f) + period=(%.4f) = %.4f\n", 
            #     (angle_defect + k_sum)/(2π), p_sum/4.0, I)
            continue
        end
        
        push!(all_indices, I)
        push!(all_angle_defects, angle_defect)
        push!(all_k_sums, k_sum)
        push!(all_p_sums, p_sum)
        
        # Filter noise
        if abs(I) > 0.15
            push!(singularities, (v_idx, I))
            println("  ⚠️  SINGULARITY DETECTED!")
            @printf("  Angle sum: %.4f rad (%.2f deg)\n", ang_sum, rad2deg(ang_sum))
            @printf("  Angle defect: %.4f rad (%.2f deg)\n", angle_defect, rad2deg(angle_defect))
            @printf("  Kappa sum: %.4f rad (%.2f deg)\n", k_sum, rad2deg(k_sum))
            @printf("  Period sum (raw): %.4f\n", p_sum)
            @printf("  Min angle: %.4f rad (%.2f deg)\n", min_angle, rad2deg(min_angle))
            @printf("  INDEX = (%.4f + %.4f)/(2π) + %.4f/4 = %.4f\n", angle_defect, k_sum, p_sum, I)
            @printf("  INDEX breakdown: geom=(%.4f) + period=(%.4f) = %.4f\n", 
                (angle_defect + k_sum)/(2π), p_sum/4.0, I)
       
        end
    end
    
    # Print statistics
    if verbose
        println("\n=== SINGULARITY STATISTICS ===")
        println("Boundary vertices skipped: $boundary_vertices_skipped")
        println("Degenerate vertices skipped: $degenerate_vertices_skipped")
        println("Total valid interior vertices analyzed: $(length(all_indices))")
        println("Singularities found: $(length(singularities))")
        if degenerate_vertices_skipped > 0
            println("\n⚠️  WARNING: $degenerate_vertices_skipped degenerate vertices detected!")
            println("This indicates poor mesh quality (nearly flat/needle triangles).")
            println("Consider remeshing with better quality constraints.")
        end
    end
    if verbose && !isempty(all_indices)
        @printf("Index range: [%.4f, %.4f]\n", minimum(all_indices), maximum(all_indices))
        @printf("Index mean: %.4f, std: %.4f\n", sum(all_indices)/length(all_indices), 
                sqrt(sum((all_indices .- sum(all_indices)/length(all_indices)).^2)/length(all_indices)))
        @printf("Angle defect mean: %.4f rad (%.2f deg)\n", 
                sum(all_angle_defects)/length(all_angle_defects), 
                rad2deg(sum(all_angle_defects)/length(all_angle_defects)))
        @printf("Kappa sum mean: %.4f rad (%.2f deg)\n", 
                sum(all_k_sums)/length(all_k_sums), 
                rad2deg(sum(all_k_sums)/length(all_k_sums)))
        @printf("Period sum mean: %.4f\n", sum(all_p_sums)/length(all_p_sums))
    end
    
    # Compute Euler characteristic and verify topological constraint
    if verbose
        println("\n=== TOPOLOGICAL VERIFICATION ===")
    end
    n_vertices = length(topo.vertices)
    n_faces = length(topo.faces)
    
    # Count edges (each edge shared by 2 faces, boundary edges by 1)
    edge_set = Set{Tuple{Int,Int}}()
    for face in topo.faces
        push!(edge_set, minmax(face[1], face[2]))
        push!(edge_set, minmax(face[2], face[3]))
        push!(edge_set, minmax(face[3], face[1]))
    end
    n_edges = length(edge_set)
    
    euler_char = n_vertices - n_edges + n_faces
    
    # Determine genus and expected index sum
    # For orientable surfaces: χ = 2 - 2g (g = genus)
    # Expected total index for cross fields (4-RoSy): (χ/2)
    genus = (2 - euler_char) / 2
    expected_index_sum = euler_char / 2.0
    
    if verbose
        @printf("Mesh topology: V=%d, E=%d, F=%d\n", n_vertices, n_edges, n_faces)
        @printf("Euler characteristic χ = V - E + F = %d\n", euler_char)
        @printf("Genus g = %.1f ", genus)
        if abs(genus - round(genus)) < 0.01
            g_int = Int(round(genus))
            if g_int == 0; println("(topological sphere)")
            elseif g_int == 1; println("(topological torus)")
            else; println("(genus-$g_int surface)")
            end
        else
            println("(non-integer genus - mesh may have boundary or be non-orientable)")
        end
        
        @printf("\nExpected singularity index sum (for cross field): %.4f\n", expected_index_sum)
    end
    
    if !isempty(singularities)
        actual_index_sum = sum(I for (_, I) in singularities)
        
        if verbose
            @printf("Actual singularity index sum: %.4f\n", actual_index_sum)
            @printf("Difference: %.4f\n", abs(actual_index_sum - expected_index_sum))
            
            if abs(actual_index_sum - expected_index_sum) < 0.3
                println("✓ Singularity sum matches topological constraint!")
            else
                println("⚠️  WARNING: Singularity sum does not match expected value!")
                println("   This could indicate:")
                println("   - Missing singularities (filtered as degenerate)")
                println("   - Incorrect period jumps")
                println("   - Mesh topology issues")
            end
        end
        
        # Show individual singularities
        if verbose
            println("\nDetected singularities:")
        end
        high_valence_count = 0
        for (v_idx, I) in singularities
            p = topo.vertices[v_idx]
            I_rounded = round(I * 4) / 4
            k = Int(round(I_rounded * 4))  # Index as k/4
            valence = 4 - k  # Quad mesh valence
            
            if verbose
                @printf("  Vertex %d at (%.4f, %.4f):\n", v_idx, p.x, p.y)
                @printf("    Index = %.4f ≈ %d/4 = %.4f\n", I, k, I_rounded)
                @printf("    Quad mesh valence = 4 - (%d) = %d", k, valence)
            end
            
            if valence < 3 || valence > 5
                if verbose
                    print(" ⚠️  HIGH VALENCE!")
                end
                high_valence_count += 1
            end
            if verbose
                println()
            end
        end
        
        if verbose && high_valence_count > 0
            println("\n⚠️  WARNING: $high_valence_count high-valence singularities detected!")
            println("Ideal quad meshes should only have valence-3 and valence-5 vertices.")
            println("High valence vertices indicate:")
            println("  - Poor triangle mesh quality (check for degenerate triangles)")
            println("  - Solver producing incorrect period jumps")
            println("  - Possible issues with boundary constraints")
            println("\nSuggestions:")
            println("  1. Try a different mesh with better quality")
            println("  2. Verify boundary constraints are aligned with mesh edges")
            println("  3. Check if $degenerate_vertices_skipped degenerate vertices are hiding the real singularities")
        end
    elseif verbose
        @printf("No singularities detected, but expected sum: %.4f\n", expected_index_sum)
        if abs(expected_index_sum) > 0.1
            println("⚠️  This suggests singularities are missing (possibly filtered as degenerate)")
        end
    end
    
    return singularities
end

function plot_frame(field, filename, title_text; cut_edges=nothing, show_singularities=false)
    """Plot a single frame for animation (simplified for speed)"""
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], aspect=DataAspect(), title=title_text)
    
    # Mesh (lighter)
    for tri in field.topology.faces
        pts = [field.topology.vertices[i] for i in tri]; push!(pts, pts[1])
        lines!(ax, [p.x for p in pts], [p.y for p in pts], color=(:gray, 0.2), linewidth=0.3)
    end
    
    # Constrained faces
    for (face_idx, _) in field.constrained_faces
        tri = field.topology.faces[face_idx]
        pts = [field.topology.vertices[i] for i in tri]; push!(pts, pts[1])
        lines!(ax, [p.x for p in pts], [p.y for p in pts], color=:orange, linewidth=1.5)
    end
    
    # Cross field vectors
    n_faces = length(field.topology.faces)
    sample_rate = max(1, n_faces ÷ 500)  # More aggressive sampling for animation
    
    xs, ys, us, vs = Float64[], Float64[], Float64[], Float64[]
    for i in 1:n_faces
        if i % sample_rate != 0; continue; end
        
        re = field.topology.face_ref_edges[i]
        p1 = field.topology.vertices[re[1]]; p2 = field.topology.vertices[re[2]]
        ref_ang = atan(p2.y - p1.y, p2.x - p1.x)
        ang = ref_ang + field.theta[i]
        
        tri = field.topology.faces[i]
        c = (field.topology.vertices[tri[1]].x + field.topology.vertices[tri[2]].x + field.topology.vertices[tri[3]].x)/3
        cy = (field.topology.vertices[tri[1]].y + field.topology.vertices[tri[2]].y + field.topology.vertices[tri[3]].y)/3
        
        push!(xs, c); push!(ys, cy)
        push!(us, cos(ang)*0.03); push!(vs, sin(ang)*0.03)
    end
    
    arrows2d!(ax, xs, ys, us, vs, color=:blue)
    
    save(filename, fig)
    return nothing
end

function plot_results(field, filename; verbose=false, cut_edges=nothing)
    fig = Figure(size=(1000, 800))
    ax = Axis(fig[1, 1], aspect=DataAspect(), title="MIQ Solution (Constrained Faces & Singularities)")
    
    # Mesh
    for tri in field.topology.faces
        pts = [field.topology.vertices[i] for i in tri]; push!(pts, pts[1])
        lines!(ax, [p.x for p in pts], [p.y for p in pts], color=:gray90, linewidth=0.5)
    end
    
    # Determine sampling rate based on mesh size
    n_faces = length(field.topology.faces)
    sample_rate = if n_faces < 1000
        1  # Show all crosses
    elseif n_faces < 10000
        2  # Show every 2nd
    elseif n_faces < 50000
        5  # Show every 5th
    elseif n_faces < 200000
        10  # Show every 10th
    else
        500  # Show every 500th for very large meshes
    end
    
    if verbose && sample_rate > 1
        println("Plotting every $(sample_rate)th cross field (mesh has $n_faces faces)")
    end
    
    # Visualize spanning tree edges (fixed edges with p=0)
    for edge in field.fixed_edges
        # Find the two faces that share this edge
        i, j = edge
        # Get centroids of both faces
        tri_i = field.topology.faces[i]
        cx_i = (field.topology.vertices[tri_i[1]].x + field.topology.vertices[tri_i[2]].x + field.topology.vertices[tri_i[3]].x)/3
        cy_i = (field.topology.vertices[tri_i[1]].y + field.topology.vertices[tri_i[2]].y + field.topology.vertices[tri_i[3]].y)/3
        
        tri_j = field.topology.faces[j]
        cx_j = (field.topology.vertices[tri_j[1]].x + field.topology.vertices[tri_j[2]].x + field.topology.vertices[tri_j[3]].x)/3
        cy_j = (field.topology.vertices[tri_j[1]].y + field.topology.vertices[tri_j[2]].y + field.topology.vertices[tri_j[3]].y)/3
        
        lines!(ax, [cx_i, cx_j], [cy_i, cy_j], color=:green, linewidth=2, alpha=0.6)
    end
    
    # Highlight constrained faces
    for (face_idx, _) in field.constrained_faces
        tri = field.topology.faces[face_idx]
        pts = [field.topology.vertices[i] for i in tri]; push!(pts, pts[1])
        lines!(ax, [p.x for p in pts], [p.y for p in pts], color=:orange, linewidth=2)
        
        # Fill constrained faces with semi-transparent color
        poly!(ax, [Point2f(p.x, p.y) for p in pts[1:3]], color=(:orange, 0.2))
    end
    
    # Field - Main direction (first cross) and secondary directions
    xs_main, ys_main, us_main, vs_main = Float64[], Float64[], Float64[], Float64[]
    xs_sec, ys_sec, us_sec, vs_sec = Float64[], Float64[], Float64[], Float64[]
    
    for i in 1:length(field.topology.faces)
        # Skip faces based on sampling rate
        if i % sample_rate != 0
            continue
        end
        
        re = field.topology.face_ref_edges[i]
        p1 = field.topology.vertices[re[1]]; p2 = field.topology.vertices[re[2]]
        ref_ang = atan(p2.y - p1.y, p2.x - p1.x)
        ang = ref_ang + field.theta[i]
        
        tri = field.topology.faces[i]
        c = (field.topology.vertices[tri[1]].x + field.topology.vertices[tri[2]].x + field.topology.vertices[tri[3]].x)/3
        cy = (field.topology.vertices[tri[1]].y + field.topology.vertices[tri[2]].y + field.topology.vertices[tri[3]].y)/3
        
        # Main direction (k=0)
        push!(xs_main, c); push!(ys_main, cy); push!(us_main, cos(ang)*0.03); push!(vs_main, sin(ang)*0.03)
        
        # Secondary directions (k=1,2,3)
        for k in 1:3
            a = ang + k*π/2
            push!(xs_sec, c); push!(ys_sec, cy); push!(us_sec, cos(a)*0.03); push!(vs_sec, sin(a)*0.03)
        end
    end
    # Draw secondary directions first (lighter color, in background)
    arrows2d!(ax, xs_sec, ys_sec, us_sec, vs_sec, color=:blue)#, arrowsize=0, lengthscale=1.0)
    # Draw main direction on top (darker, more prominent)
    arrows2d!(ax, xs_main, ys_main, us_main, vs_main, color=:red)#, arrowsize=0, lengthscale=1.0, linewidth=2)
    
    # Singularities
    sings = compute_singularities(field; verbose=verbose)
    pos = Point2f[]; neg = Point2f[]
    for (v, I) in sings
        p = field.topology.vertices[v]
        pt = Point2f(p.x, p.y)
        if I > 0; push!(pos, pt); else; push!(neg, pt); end
    end
    scatter!(ax, pos, color=:red, markersize=15, strokecolor=:black, strokewidth=1, label="+ Index")
    scatter!(ax, neg, color=:cyan, markersize=15, strokecolor=:black, strokewidth=1, label="- Index")
    

    cx, cy = Float64[], Float64[]
    for (u, v) in cut_edges
        p1, p2 = field.topology.vertices[u], field.topology.vertices[v]
        push!(cx, p1.x, p2.x, NaN)
        push!(cy, p1.y, p2.y, NaN)
    end
    lines!(ax, cx, cy, color=:red, linewidth=5.0)


    # Add mesh cut edges to legend
    if cut_edges !== nothing && !isempty(cut_edges)
        # Dummy plot for legend entry
        lines!(ax, [NaN], [NaN], color=:red, linewidth=5.0, label="Mesh Cut Edges")
    end

    # Add constrained faces to legend
    if !isempty(field.constrained_faces)
        # Dummy plot for legend entry
        lines!(ax, [NaN], [NaN], color=:orange, linewidth=2, label="Constrained Faces")
    end
    
    # Add spanning tree to legend
    if !isempty(field.fixed_edges)
        lines!(ax, [NaN], [NaN], color=:green, linewidth=2, alpha=0.6, label="Spanning Tree")
    end
    
    if !isempty(pos) || !isempty(neg) || !isempty(field.constrained_faces) || !isempty(field.fixed_edges)
        axislegend(ax)
    end
    


    save(filename, fig)
    if verbose
        println("Saved visualization to $filename")
        println("Found $(length(sings)) singularities.")
    end
end


# ==============================================================================
# 6. ANIMATION UTILITIES
# ==============================================================================

# Animation functions removed - use create_gif.py Python script instead
# See ANIMATION_README.md for instructions

# ==============================================================================
# 7. MESH CUTTING (TOPOLOGICAL DISK)
# ==============================================================================

function compute_cut_graph(topo::MeshTopology)
    println("\n--- Computing Cut Graph (Dual Spanning Tree) ---")
    n_faces = length(topo.faces)
    
    # 1. Build Dual Spanning Tree
    # We use BFS to grow a tree connecting all faces.
    # Any edge we CROSS is a "Tree Edge".
    # Any edge we DO NOT CROSS (but connects two visited faces) is a "Cut Edge".
    
    visited_face = falses(n_faces)
    tree_edges = Set{Tuple{Int,Int}}() # Edges we crossed
    
    # Pick a random start face (or Face 1)
    queue = Int[1]
    visited_face[1] = true
    
    while !isempty(queue)
        f_curr = popfirst!(queue)
        
        # Check all neighbors
        for (f_neigh, edge_key) in topo.dual_adj[f_curr]
            if !visited_face[f_neigh]
                visited_face[f_neigh] = true
                push!(tree_edges, edge_key)
                push!(queue, f_neigh)
            end
        end
    end
    
    # 2. Identify Initial Cut Edges (Primal of non-tree edges)
    # These are all internal edges NOT in the spanning tree.
    raw_cut_edges = Set{Tuple{Int,Int}}()
    
    for i in 1:n_faces
        for (j, key) in topo.dual_adj[i]
            # Only consider internal edges (i < j to avoid duplicates)
            if i < j
                if !(key in tree_edges)
                    push!(raw_cut_edges, key)
                end
            end
        end
    end
    
    println("Initial Cut Graph size: $(length(raw_cut_edges)) edges")
    
    # 3. Prune Open Paths (The "Iterative Reduction" step)
    # A cut edge is "useless" if it doesn't help open a loop.
    # Visually, these are "dead ends" in the cut graph (valence 1).
    # We remove them unless they touch the mesh boundary.
    
    # Identify mesh boundary vertices
    boundary_verts = Set{Int}()
    edge_counts = Dict{Tuple{Int,Int}, Int}()
    
    # Simple boundary check: edges with only 1 neighbor
    all_half_edges = Dict{Tuple{Int,Int}, Int}()
    for tri in topo.faces
        e1, e2, e3 = minmax(tri[1],tri[2]), minmax(tri[2],tri[3]), minmax(tri[3],tri[1])
        all_half_edges[e1] = get(all_half_edges, e1, 0) + 1
        all_half_edges[e2] = get(all_half_edges, e2, 0) + 1
        all_half_edges[e3] = get(all_half_edges, e3, 0) + 1
    end
    for (e, count) in all_half_edges
        if count == 1
            push!(boundary_verts, e[1])
            push!(boundary_verts, e[2])
        end
    end
    
    # Iterative Pruning
    final_cuts = copy(raw_cut_edges)
    changed = true
    
    while changed
        changed = false
        
        # Build adjacency of the CUT GRAPH
        cut_adj = Dict{Int, Vector{Int}}()
        for (u, v) in final_cuts
            if !haskey(cut_adj, u); cut_adj[u] = Int[]; end
            if !haskey(cut_adj, v); cut_adj[v] = Int[]; end
            push!(cut_adj[u], v)
            push!(cut_adj[v], u)
        end
        
        # Find leaves (valence 1) that are NOT on boundary
        to_remove = Set{Tuple{Int,Int}}()
        for (v, neighbors) in cut_adj
            if length(neighbors) == 1 && !(v in boundary_verts)
                # This is a dead end inside the mesh. Prune it.
                neighbor = neighbors[1]
                edge = minmax(v, neighbor)
                push!(to_remove, edge)
                changed = true
            end
        end
        
        # Apply removal
        if !isempty(to_remove)
            setdiff!(final_cuts, to_remove)
            # println("  Pruned $(length(to_remove)) dead-end edges...")
        end
    end
    
    println("Final Cut Graph size: $(length(final_cuts)) edges")
    return final_cuts
end



function propagate_orientations(field::CrossField, cut_edges::Set{Tuple{Int,Int}})
    topo = field.topology
    n_faces = length(topo.faces)
    
    # r[f] stores the integer rotation (in 90-degree steps) for face f
    # such that the parameter orientation matches the global frame.
    r = zeros(Int, n_faces)
    visited = falses(n_faces)
    
    # BFS Queue
    q = Queue{Int}()
    
    # Start from a root face (e.g., face 1)
    root = 1
    enqueue!(q, root)
    visited[root] = true
    r[root] = 0 # Canonical base orientation
    
    println("\n--- Propagating Global Orientation ---")
    
    while !isempty(q)
        u = dequeue!(q)
        
        # Check all neighbors
        for (v, edge_key) in topo.dual_adj[u]
            # CRITICAL: Do not propagate orientation across cut edges!
            # The cut edges act as boundaries for the parameterization domain.
            if edge_key in cut_edges
                continue
            end
            
            if !visited[v]
                # Retrieve the cross field period jump p_uv
                # Note: field.period_jumps is stored with key minmax(u, v)
                p = get(field.period_jumps, edge_key, 0)
                
                # Apply compatibility rule derived from: theta_j - theta_i approx p * pi/2
                # We want: (theta_j + r_j) - (theta_i + r_i) = 0
                # Result: r_j = r_i - p (if u < v)
                if u < v
                    r[v] = r[u] - p
                else
                    r[v] = r[u] + p
                end
                
                visited[v] = true
                enqueue!(q, v)
            end
        end
    end
    
    # Compute the transition rotations for the cut edges
    # These are the 'k' values in the transition function: (u', v') = Rot^k (u, v)
    cut_rotations = Dict{Tuple{Int,Int}, Int}()
    
    for edge in cut_edges
        u, v = (edge[1] < edge[2]) ? (edge[1], edge[2]) : (edge[2], edge[1])
        
        # Ensure we reached both sides (sanity check for connected mesh)
        if visited[u] && visited[v]
            p = get(field.period_jumps, edge, 0)
            
            # The rotation mismatch across the cut
            # k = r_v - (r_u - p)  =>  k = r_v - r_u + p
            k = r[v] - r[u] + p
            
            # Normalize k to [-1, 2] usually, or just keep integer
            cut_rotations[(u, v)] = k
        else
            println("Warning: Cut edge $edge connects unvisited faces (disconnected component?)")
        end
    end
    
    println("Propagated orientation to $(count(visited)) faces.")
    println("Computed rotations for $(length(cut_rotations)) cut edges.")
    
    return r, cut_rotations
end


function plot_smooth_global_field(field::CrossField, rotations::Vector{Int}, filename::String; cut_edges=nothing, verbose=false)
    # Create Figure
    fig = Figure(size=(1000, 800))
    ax = Axis(fig[1, 1], aspect=DataAspect(), title="Globally Smoothed 'u' Direction (Seams at Cuts)")
    
    # 1. Plot Mesh Wireframe (faint background)
    for tri in field.topology.faces
        pts = [field.topology.vertices[i] for i in tri]; push!(pts, pts[1])
        lines!(ax, [p.x for p in pts], [p.y for p in pts], color=(:gray, 0.3), linewidth=0.5)
    end

    # 2. Plot Cut Edges (Red Lines)
    if cut_edges !== nothing
        cx, cy = Float64[], Float64[]
        for (u, v) in cut_edges
            p1, p2 = field.topology.vertices[u], field.topology.vertices[v]
            push!(cx, p1.x, p2.x, NaN)
            push!(cy, p1.y, p2.y, NaN)
        end
        lines!(ax, cx, cy, color=:red, linewidth=4.0, label="Cut Graph")
    end

    # 3. Plot Smooth Vectors
    # We plot ONE vector per face: The local cross direction rotated by the integer rotation
    xs, ys, us, vs = Float64[], Float64[], Float64[], Float64[]
    
    # Sampling for large meshes
    n_faces = length(field.topology.faces)
    sample_rate = n_faces > 10000 ? (n_faces ÷ 2000) : 1
    
    for i in 1:n_faces
        if i % sample_rate != 0; continue; end

        # A. Get Geometry of Reference Edge
        re = field.topology.face_ref_edges[i]
        p1 = field.topology.vertices[re[1]]
        p2 = field.topology.vertices[re[2]]
        ref_ang = atan(p2.y - p1.y, p2.x - p1.x)
        
        # B. Get Local Theta
        local_theta = field.theta[i]
        
        # C. Get Integer Rotation (from propagation)
        # r is the number of 90-degree turns
        r = rotations[i]
        
        # D. Compute Global Smooth Angle
        # Angle = Reference + Local_Theta + (Rotation * pi/2)
        global_angle = ref_ang + local_theta + (r * π / 2.0)
        
        # E. Centroid
        tri = field.topology.faces[i]
        c = (field.topology.vertices[tri[1]].x + field.topology.vertices[tri[2]].x + field.topology.vertices[tri[3]].x)/3
        cy = (field.topology.vertices[tri[1]].y + field.topology.vertices[tri[2]].y + field.topology.vertices[tri[3]].y)/3
        
        push!(xs, c); push!(ys, cy)
        push!(us, cos(global_angle) * 0.04) # Scale length as needed
        push!(vs, sin(global_angle) * 0.04)
    end

    arrows2d!(ax, xs, ys, us, vs, color=:blue, label="Global 'u' Field")
    
    # 4. Singularities (Optional context)
    sings = compute_singularities(field)
    pos = Point2f[]; neg = Point2f[]
    for (v, I) in sings
        p = field.topology.vertices[v]
        if I > 0; push!(pos, Point2f(p.x, p.y)); else; push!(neg, Point2f(p.x, p.y)); end
    end
    if !isempty(pos); scatter!(ax, pos, color=:red, markersize=12, label="Singularity (+)"); end
    if !isempty(neg); scatter!(ax, neg, color=:cyan, markersize=12, label="Singularity (-)"); end

    axislegend(ax)
    save(filename, fig)
    if verbose; println("Saved smooth field plot to $filename"); end
end

# ==============================================================================
# 8. GLOBAL PARAMETERIZATION (U, V EXTRACTION)
# ==============================================================================

"""
Compute global (u, v) parameterization coordinates for each triangle vertex.
Uses Poisson integration to recover the scalar fields from the cross field.
"""
function compute_global_parameterization(field::CrossField, rotations::Vector{Int}, cut_edges::Set{Tuple{Int,Int}})
    println("\n--- Computing Global Parameterization (u, v) ---")
    
    topo = field.topology
    n_verts = length(topo.vertices)
    n_faces = length(topo.faces)
    
    # Build vertex-to-face adjacency
    vert_to_faces = [Int[] for _ in 1:n_verts]
    for (f_idx, tri) in enumerate(topo.faces)
        push!(vert_to_faces[tri[1]], f_idx)
        push!(vert_to_faces[tri[2]], f_idx)
        push!(vert_to_faces[tri[3]], f_idx)
    end
    
    # Build edge list (excluding cut edges)
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
    
    # Build least-squares system: minimize ||∇u - target_u||² + ||∇v - target_v||²
    # For each edge, we want the gradient to match the cross field direction
    
    I_rows = Int[]
    J_cols = Int[]
    vals = Float64[]
    b_u = Float64[]
    b_v = Float64[]
    
    equation_idx = 0
    
    for (v1, v2) in edges
        # Find a face containing this edge
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
        
        # Get global field direction at this face
        re = topo.face_ref_edges[face_idx]
        p1 = topo.vertices[re[1]]
        p2 = topo.vertices[re[2]]
        ref_ang = atan(p2.y - p1.y, p2.x - p1.x)
        
        local_theta = field.theta[face_idx]
        r = rotations[face_idx]
        global_angle = ref_ang + local_theta + (r * π / 2.0)
        
        # Target directions
        u_dir = (cos(global_angle), sin(global_angle))
        v_dir = (-sin(global_angle), cos(global_angle))
        
        # Edge vector in 3D space
        p_v1 = topo.vertices[v1]
        p_v2 = topo.vertices[v2]
        edge_vec = (p_v2.x - p_v1.x, p_v2.y - p_v1.y)
        edge_len = sqrt(edge_vec[1]^2 + edge_vec[2]^2)
        
        if edge_len < 1e-10
            continue
        end
        
        # Project edge onto u and v
        target_du = (edge_vec[1] * u_dir[1] + edge_vec[2] * u_dir[2])
        target_dv = (edge_vec[1] * v_dir[1] + edge_vec[2] * v_dir[2])
        
        # Add equation: u[v2] - u[v1] ≈ target_du
        equation_idx += 1
        push!(I_rows, equation_idx); push!(J_cols, v1); push!(vals, -1.0)
        push!(I_rows, equation_idx); push!(J_cols, v2); push!(vals, 1.0)
        push!(b_u, target_du)
        push!(b_v, target_dv)
    end
    
    # Add constraint: fix first vertex
    equation_idx += 1
    push!(I_rows, equation_idx); push!(J_cols, 1); push!(vals, 1e3)  # Strong constraint
    push!(b_u, 0.0)
    push!(b_v, 0.0)
    
    # Solve least-squares system
    A = sparse(I_rows, J_cols, vals, equation_idx, n_verts)
    
    println("Solving least-squares system ($equation_idx equations, $n_verts unknowns)...")
    
    # Solve using normal equations: A'A x = A'b
    ATA = A' * A
    ATb_u = A' * b_u
    ATb_v = A' * b_v
    
    u_coords = ATA \ ATb_u
    v_coords = ATA \ ATb_v
    
    println("Computed smooth parameterization for all $n_verts vertices")
    
    return u_coords, v_coords
end

"""
Smooth the parameterization using Poisson smoothing to reduce distortion.
"""
function smooth_parameterization(topo::MeshTopology, u_coords::Vector{Float64}, v_coords::Vector{Float64}, 
                                 cut_edges::Set{Tuple{Int,Int}}; iterations=5)
    println("\n--- Smoothing Parameterization ---")
    
    n_verts = length(u_coords)
    u_smooth = copy(u_coords)
    v_smooth = copy(v_coords)
    
    # Build vertex adjacency (skip cut edges)
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
    
    # Identify boundary vertices (don't smooth these)
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
    
    # Iterative smoothing
    for iter in 1:iterations
        u_new = copy(u_smooth)
        v_new = copy(v_smooth)
        
        for v_idx in 1:n_verts
            if v_idx in boundary_verts || isempty(vert_neighbors[v_idx])
                continue
            end
            
            # Average with neighbors
            neighbors = vert_neighbors[v_idx]
            u_avg = sum(u_smooth[n] for n in neighbors) / length(neighbors)
            v_avg = sum(v_smooth[n] for n in neighbors) / length(neighbors)
            
            # Weighted average (0.5 = moderate smoothing)
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

"""
Extract quad mesh by tracing isolines through the triangle mesh.
This creates new vertices at isoline intersections.
"""
function extract_quad_mesh_isolines(topo::MeshTopology, u_coords::Vector{Float64}, v_coords::Vector{Float64}; 
                                   grid_spacing=nothing, verbose=false)
    println("\n--- Extracting Quad Mesh via Isoline Tracing ---")
    
    u_min, u_max = minimum(u_coords), maximum(u_coords)
    v_min, v_max = minimum(v_coords), maximum(v_coords)
    u_range = u_max - u_min
    v_range = v_max - v_min
    
    println("U range: [$u_min, $u_max], span: $u_range")
    println("V range: [$v_min, $v_max], span: $v_range")
    
    # Use SEPARATE grid spacings for u and v based on their ranges
    if grid_spacing === nothing
        # Target number of divisions
        target_divs = 8
        u_spacing = u_range / target_divs
        v_spacing = v_range / target_divs
    else
        u_spacing = grid_spacing
        v_spacing = grid_spacing
    end
    
    println("Grid spacing: u=$u_spacing, v=$v_spacing")
    
    # Map vertices to grid cells with separate u and v indices
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
    
    println("Mapped vertices to $(length(grid_verts)) grid cells")
    
    # Build quads from 2x2 grid cells
    quads = Tuple{Int,Int,Int,Int}[]
    
    # Get all grid cell coordinates
    all_cells = collect(keys(grid_verts))
    
    for (u_g, v_g) in all_cells
        # Check if we have all 4 corners for a quad
        corners = [
            (u_g, v_g),
            (u_g+1, v_g),
            (u_g+1, v_g+1),
            (u_g, v_g+1)
        ]
        
        if all(haskey(grid_verts, c) && !isempty(grid_verts[c]) for c in corners)
            # Use the vertex closest to the grid cell center
            v1 = grid_verts[corners[1]][1]
            v2 = grid_verts[corners[2]][1]
            v3 = grid_verts[corners[3]][1]
            v4 = grid_verts[corners[4]][1]
            
            # Quality check: all 4 vertices must be distinct
            if length(Set([v1, v2, v3, v4])) == 4
                # Additional quality check: quad shouldn't be too distorted
                # (we could check aspect ratio, angles, etc.)
                push!(quads, (v1, v2, v3, v4))
            end
        end
    end
    
    println("Extracted $(length(quads)) quadrilaterals")
    
    return quads, (u_spacing, v_spacing)
end

"""
Simple vertex snapping method (original approach).
"""
function extract_quad_mesh(topo::MeshTopology, u_coords::Vector{Float64}, v_coords::Vector{Float64}; 
                          grid_spacing=nothing, verbose=false)
    return extract_quad_mesh_isolines(topo, u_coords, v_coords; grid_spacing=grid_spacing, verbose=verbose)
end

"""  
Plot the extracted quad mesh with quality analysis.
"""
function plot_quad_mesh(topo::MeshTopology, u_coords::Vector{Float64}, v_coords::Vector{Float64}, 
                        quads::Vector{Tuple{Int,Int,Int,Int}}, filename::String; verbose=false)
    println("\n--- Plotting Quad Mesh ---")
    
    fig = Figure(size=(1600, 500))
    
    # Left: UV parameter space
    ax1 = Axis(fig[1, 1], aspect=DataAspect(), title="UV Parameter Space")
    
    # Plot quads in parameter space
    for quad in quads
        v1, v2, v3, v4 = quad
        us = [u_coords[v1], u_coords[v2], u_coords[v3], u_coords[v4], u_coords[v1]]
        vs = [v_coords[v1], v_coords[v2], v_coords[v3], v_coords[v4], v_coords[v1]]
        lines!(ax1, us, vs, color=:blue, linewidth=1.5)
    end
    
    # Scatter all vertices
    scatter!(ax1, u_coords, v_coords, color=:red, markersize=4)
    
    # Middle: Original mesh with parameterization
    ax2 = Axis(fig[1, 2], aspect=DataAspect(), title="Parameterization on Mesh")
    
    # Plot triangles faintly
    for tri in topo.faces
        pts = [topo.vertices[i] for i in tri]
        xs = [pts[1].x, pts[2].x, pts[3].x, pts[1].x]
        ys = [pts[1].y, pts[2].y, pts[3].y, pts[1].y]
        lines!(ax2, xs, ys, color=(:gray, 0.15), linewidth=0.3)
    end
    
    # Scatter vertices colored by u-coordinate
    xs = [v.x for v in topo.vertices]
    ys = [v.y for v in topo.vertices]
    scatter!(ax2, xs, ys, color=u_coords, colormap=:viridis, markersize=4)
    
    # Right: Extracted quad mesh
    ax3 = Axis(fig[1, 3], aspect=DataAspect(), title="Quad Mesh ($(length(quads)) quads)")
    
    # Plot original mesh very faintly
    for tri in topo.faces
        pts = [topo.vertices[i] for i in tri]
        xs = [pts[1].x, pts[2].x, pts[3].x, pts[1].x]
        ys = [pts[1].y, pts[2].y, pts[3].y, pts[1].y]
        lines!(ax3, xs, ys, color=(:gray, 0.08), linewidth=0.2)
    end
    
    # Compute quad quality (aspect ratio) and color code
    quad_qualities = Float64[]
    
    for quad in quads
        v1, v2, v3, v4 = quad
        p1, p2, p3, p4 = topo.vertices[v1], topo.vertices[v2], topo.vertices[v3], topo.vertices[v4]
        
        # Compute edge lengths
        e1 = sqrt((p2.x-p1.x)^2 + (p2.y-p1.y)^2)
        e2 = sqrt((p3.x-p2.x)^2 + (p3.y-p2.y)^2)
        e3 = sqrt((p4.x-p3.x)^2 + (p4.y-p3.y)^2)
        e4 = sqrt((p1.x-p4.x)^2 + (p1.y-p4.y)^2)
        
        # Quality metric: min/max edge length ratio
        min_edge = min(e1, e2, e3, e4)
        max_edge = max(e1, e2, e3, e4)
        quality = max_edge > 1e-10 ? min_edge / max_edge : 0.0
        
        push!(quad_qualities, quality)
    end
    
    # Plot quads colored by quality
    for (q_idx, quad) in enumerate(quads)
        v1, v2, v3, v4 = quad
        pts = [topo.vertices[v1], topo.vertices[v2], topo.vertices[v3], topo.vertices[v4]]
        xs = [pts[1].x, pts[2].x, pts[3].x, pts[4].x, pts[1].x]
        ys = [pts[1].y, pts[2].y, pts[3].y, pts[4].y, pts[1].y]
        
        # Color by quality (green = good, red = bad)
        quality = quad_qualities[q_idx]
        color = quality > 0.5 ? :green : (quality > 0.3 ? :orange : :red)
        
        lines!(ax3, xs, ys, color=color, linewidth=2)
    end
    
    # Add colorbar for reference
    Colorbar(fig[1, 4], limits=(minimum(u_coords), maximum(u_coords)), 
             colormap=:viridis, label="u-coordinate")
    
    save(filename, fig)
    
    if verbose
        println("Saved quad mesh visualization to $filename")
        if !isempty(quad_qualities)
            @printf("Quad quality - min: %.3f, avg: %.3f, max: %.3f\n", 
                    minimum(quad_qualities), sum(quad_qualities)/length(quad_qualities), maximum(quad_qualities))
        end
    end
end# ==============================================================================
# 6. MAIN
# ==============================================================================

function compute_parameterization_least_squares(field, rotations, cut_edges)
    println("\n--- Solving Global Parameterization (Least Squares) ---")
    
    topo = field.topology
    n_verts = length(topo.vertices)
    
    # 1. Build the Linear System
    # We want to minimize: sum over edges [ (u_j - u_i - target_du)^2 ]
    # This is an overdetermined system Ax = b.
    # We solve the Normal Equations: (A'A)x = A'b
    
    # Estimate number of edges for allocation (Euler: E approx 3V)
    est_edges = n_verts * 3
    I = Int[]; J = Int[]; V = Float64[]
    b_u = zeros(est_edges)
    b_v = zeros(est_edges)
    
    row_idx = 0
    
    # Helper to check if edge is a cut
    cut_set = Set{Tuple{Int,Int}}()
    for (u, v) in cut_edges
        push!(cut_set, minmax(u, v))
    end
    
    # Iterate over all VALID edges (those NOT in the cut graph)
    # We can iterate faces and look at edges, but we need unique edges.
    seen_edges = Set{Tuple{Int,Int}}()
    
    for i in 1:length(topo.faces)
        tri = topo.faces[i]
        
        # Get Global Angle for this face
        # Angle = Reference + Local_Theta + (Rotation * pi/2)
        re = topo.face_ref_edges[i]
        p1_ref = topo.vertices[re[1]]
        p2_ref = topo.vertices[re[2]]
        ref_ang = atan(p2_ref.y - p1_ref.y, p2_ref.x - p1_ref.x)
        
        global_angle = ref_ang + field.theta[i] + (rotations[i] * π / 2.0)
        
        # Basis vectors for this face
        dir_u = (cos(global_angle), sin(global_angle))
        dir_v = (-sin(global_angle), cos(global_angle))
        
        # Process 3 edges of the triangle
        face_edges = [(tri[1], tri[2]), (tri[2], tri[3]), (tri[3], tri[1])]
        
        for (u, v) in face_edges
            key = minmax(u, v)
            
            # Skip cut edges (they are boundaries in the parameter domain)
            if key in cut_set
                continue
            end
            
            # Optimization: Only add each edge once. 
            # Note: For weighted Laplacian, we might want to average the target_du from both faces.
            # Here we just take the first face's opinion for simplicity.
            if key in seen_edges
                continue
            end
            push!(seen_edges, key)
            
            row_idx += 1
            
            # Geometry
            p_u = topo.vertices[u]
            p_v = topo.vertices[v]
            dx = p_v.x - p_u.x
            dy = p_v.y - p_u.y
            
            # Target differences
            target_du = dx * dir_u[1] + dy * dir_u[2]
            target_dv = dx * dir_v[1] + dy * dir_v[2]
            
            # Equation: u_target - u_source = target
            # Add +1 at column v, -1 at column u
            push!(I, row_idx); push!(J, v); push!(V, 1.0)
            push!(I, row_idx); push!(J, u); push!(V, -1.0)
            
            # Ensure proper sign direction based on edge traversal
            # We computed target based on v - u
            if row_idx > length(b_u)
                resize!(b_u, length(b_u)*2)
                resize!(b_v, length(b_v)*2)
            end
            b_u[row_idx] = target_du
            b_v[row_idx] = target_dv
        end
    end
    
    # 2. Add Soft Constraints (Pin one vertex to prevent floating)
    # Pin vertex 1 to (0,0)
    row_idx += 1
    push!(I, row_idx); push!(J, 1); push!(V, 1.0) 
    # (b is effectively 0.0 for this row)
    
    # Resize b to actual count
    b_u = b_u[1:row_idx]
    b_v = b_v[1:row_idx]
    
    # 3. Solve System
    # Construct Sparse Matrix A (Equations x Vertices)
    A = sparse(I, J, V, row_idx, n_verts)
    
    println("  System size: $(size(A)). Solving Normal Equations...")
    
    # Solve (A'A)x = A'b
    # In Julia, backslash on non-square sparse matrices usually does QR or Least Squares automatically.
    # But explicitly forming normal equations is often faster for 2D Laplacians.
    AtA = A' * A
    Atb_u = A' * b_u
    Atb_v = A' * b_v
    
    # Add regularization for stability (create sparse diagonal matrix)
    AtA_reg = AtA + 1.0e-6 * spdiagm(0 => ones(n_verts))
    
    # Factorize once
    F = cholesky(AtA_reg)
    
    u_sol = F \ Atb_u
    v_sol = F \ Atb_v
    
    println("  Range U: [$(minimum(u_sol)), $(maximum(u_sol))]")
    println("  Range V: [$(minimum(v_sol)), $(maximum(v_sol))]")
    
    return u_sol, v_sol
end


function main()
    # Replace with your file
    # filename = "triangulations/mesh_airfoil_dae11.su2" 
    filename = "triangulations/regular-square-8x8.msh"
    # filename = "triangulations/disk-radial-fine.msh"
    # filename = "triangulations/crmhl_test.su2"
    
    if !isfile(filename)
        println("File not found: $filename")
        return
    end
    
    println("Reading mesh...")
    verts, faces = read_mesh(filename)
    topo = build_topology(verts, faces)
    
    println("Detecting boundaries...")
    constraints = compute_boundary_constraints(topo)
    # println("Constrained $(length(constraints)) boundary faces.")
    # constraints = Dict{Int, Float64}(1 => 0.0)  # Example: Fix face 1 to angle 0.0
    
    # Set verbosity and animation parameters
    verbose = false
    animate = false  # Set to true to generate animation (disable for now to focus on quad mesh)
    frame_dir = "output/animation_frames"
    frame_interval = 5  # Save a frame every N iterations (adjust based on mesh size)

    field = initialize_field(topo, constraints)
    
    # Solve with animation
    println("Solving cross field...")
    solve_greedy!(field; verbose=true, animate=animate, frame_dir=frame_dir, frame_interval=frame_interval)
    println("Cross field computed successfully.")


    # compute singularities
    println("Computing singularities...")
    sings = compute_singularities(field; verbose=verbose)
    println("Found $(length(sings)) singularities.")

    # Compute Cut Graph
    cut_edges = compute_cut_graph(topo)

    # 1. Propagate Orientation
    face_rotations, cut_transitions = propagate_orientations(field, cut_edges)
    
    println("Plotting globally smooth field...")
    plot_smooth_global_field(field, face_rotations, "global_smooth_field.png"; 
                             cut_edges=cut_edges, verbose=true)
    
    # # 2. Compute Global Parameterization (u, v)
    # u_coords, v_coords = compute_global_parameterization(field, face_rotations, cut_edges)
    
    # # 2b. Smooth the parameterization to reduce distortion
    # u_smooth, v_smooth = smooth_parameterization(topo, u_coords, v_coords, cut_edges; iterations=10)
    
    # # 3. Extract Quad Mesh (using smoothed coordinates)
    # quads, grid_spacing = extract_quad_mesh(topo, u_smooth, v_smooth; verbose=true)

    # 2. Compute Global Parameterization (Use Least Squares now!)
    u_coords, v_coords = compute_parameterization_least_squares(field, face_rotations, cut_edges)
    u_smooth = u_coords
    v_smooth = v_coords
    # 3. Extract Quad Mesh
    quads, grid_spacing = extract_quad_mesh(topo, u_coords, v_coords; verbose=true)

    # Note: Using smoothed parameterization produces better quad quality
    println("  Using smoothed parameterization for better quality")
    
    # 4. Visualize Quad Mesh
    println("Plotting quad mesh...")
    plot_quad_mesh(topo, u_smooth, v_smooth, quads, "quad_mesh_result.png"; verbose=true)
    
    # Visualize
    println("Generating visualization...")
    plot_results(field, "miq_solution.png"; verbose=verbose, cut_edges=cut_edges)
    println("Complete! Visualization saved to miq_solution.png")
    
    # Create animation GIF
    if animate
        println("\n--- Creating Animation ---")
        println("Method 1: Trying Python script...")
        
        # Try Python script first
        try
            python_path = "C:\\Users\\admin\\AppData\\Local\\Programs\\Python\\Python311\\python.exe"
            run(`$python_path create_gif.py`)
            println("✓ GIF created successfully with Python!")
        catch e
            println("Python method failed. Frames are available in: $frame_dir")
            println("\nTo create GIF manually:")
            println("  Option 1 (ImageMagick): magick -delay 20 -loop 0 \"$frame_dir\\frame_*.png\" greedy_solver_animation.gif")
            println("  Option 2 (FFmpeg): ffmpeg -framerate 5 -pattern_type glob -i \"$frame_dir\\frame_*.png\" animation.mp4")
            println("  Option 3 (Python): python create_gif.py")
            println("\nSee ANIMATION_README.md for more details.")
        end
    end
end

main()