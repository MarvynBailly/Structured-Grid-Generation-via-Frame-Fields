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

function compute_boundary_constraints(topo; ccw_boundary=true)
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
                third_idx = if k == 1; tri[3]; elseif k == 2; tri[1]; else; tri[2]; end
                p3 = topo.vertices[third_idx]
                
                # 2. Get Boundary Edge Coordinates
                p1 = topo.vertices[edge[1]]
                p2 = topo.vertices[edge[2]]
                
                # 3. Compute Edge Vector and Tangent
                dx = p2.x - p1.x
                dy = p2.y - p1.y
                global_tan = atan(dy, dx)
                
                # 4. Apply Direction Logic
                if ccw_boundary
                    # CCW Mode: Check if third vertex is on the left of the edge
                    # Cross product: (p2 - p1) × (p3 - p1)
                    cross_z = dx * (p3.y - p1.y) - dy * (p3.x - p1.x)
                    
                    # If cross_z < 0, the third vertex is on the right (clockwise)
                    # We need to flip the tangent by π to make it go CCW
                    if cross_z < 0.0
                        global_tan += π
                    end
                else
                    # Original Y-based Logic
                    if p3.y < 0.0
                        # Bottom Surface: Force +X direction
                        if dx < 0.0
                            global_tan += π
                        end
                    elseif p3.y > 0.0
                        # Top Surface: Force -X direction
                        if dx > 0.0
                            global_tan += π
                        end
                    end
                end
                
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
        is_degenerate =(min_angle < MIN_ANGLE_THRESHOLD) #|| (ang_sum < MIN_ANGLE_SUM_THRESHOLD)
        
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

function create_gif_native(frame_dir::String, output_file::String; fps=10)
    """
    Create an animated GIF from frames using Julia's FileIO and ImageIO packages.
    
    Arguments:
    - frame_dir: Directory containing frame_####.png files
    - output_file: Output GIF filename
    - fps: Frames per second (default: 10)
    """
    println("\n--- Creating Animated GIF (Native Julia) ---")
    
    # Check if frames exist
    frame_files = filter(f -> endswith(f, ".png"), readdir(frame_dir))
    if isempty(frame_files)
        println("ERROR: No PNG frames found in $frame_dir")
        return false
    end
    
    sort!(frame_files)  # Ensure frames are in order
    println("Found $(length(frame_files)) frames")
    
    try
        # Try to use FileIO for GIF creation
        println("Loading frames...")
        images = []
        for (i, fname) in enumerate(frame_files)
            fpath = joinpath(frame_dir, fname)
            img = load(fpath)
            push!(images, img)
            if i % 10 == 0
                println("  Loaded $i/$(length(frame_files)) frames")
            end
        end
        
        println("Creating GIF with $fps fps...")
        save(output_file, images, fps=fps)
        
        println("✓ Successfully created: $output_file")
        
        # Report file size
        file_size = filesize(output_file) / (1024 * 1024)  # MB
        @printf("  File size: %.2f MB\n", file_size)
        
        return true
    catch e
        println("ERROR: Failed to create GIF with native method")
        println("Error: $e")
        println("\nAlternative: Use ImageMagick manually:")
        println("  Install from: https://imagemagick.org/")
        println("  Then run: magick -delay $(round(Int, 100/fps)) -loop 0 \"$frame_dir\\frame_*.png\" $output_file")
        return false
    end
end

function create_gif(frame_dir::String, output_file::String; fps=10, loop=0)
    """
    Create an animated GIF from frames in frame_dir.
    Tries ImageMagick first, falls back to native Julia method.
    
    Arguments:
    - frame_dir: Directory containing frame_####.png files
    - output_file: Output GIF filename
    - fps: Frames per second (default: 10)
    - loop: Number of loops (0 = infinite, default: 0) - only for ImageMagick
    """
    println("\n--- Creating Animated GIF ---")
    
    # Check if frames exist
    frames = filter(f -> endswith(f, ".png"), readdir(frame_dir))
    if isempty(frames)
        println("ERROR: No PNG frames found in $frame_dir")
        return false
    end
    
    sort!(frames)  # Ensure frames are in order
    println("Found $(length(frames)) frames")
    
    # Try ImageMagick first
    try
        # Calculate delay for ImageMagick (in 1/100th of a second)
        delay = round(Int, 100 / fps)
        
        # Build ImageMagick command
        frame_pattern = joinpath(frame_dir, "frame_*.png")
        cmd = `magick -delay $delay -loop $loop $frame_pattern $output_file`
        
        println("Trying ImageMagick...")
        run(cmd)
        println("✓ Successfully created with ImageMagick: $output_file")
        
        # Report file size
        file_size = filesize(output_file) / (1024 * 1024)  # MB
        @printf("  File size: %.2f MB\n", file_size)
        
        return true
    catch e
        println("ImageMagick not available, trying native Julia method...")
        return create_gif_native(frame_dir, output_file; fps=fps)
    end
end

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
    
    # r[f] stores the integer rotation (in 90-degree steps) for face f.
    # r is an integer "winding number" that flattens the connection.
    r = zeros(Int, n_faces)
    visited = falses(n_faces)
    
    # Queue for BFS
    q = Queue{Int}()
    
    println("\n--- Propagating Global Orientation ---")
    
    # Outer loop: Handle multiple connected components
    for start_node in 1:n_faces
        if visited[start_node]
            continue
        end
        
        # Start a new component
        enqueue!(q, start_node)
        visited[start_node] = true
        r[start_node] = 0 # Canonical base orientation for this component
        
        while !isempty(q)
            u = dequeue!(q)
            
            # Check all neighbors
            for (v, edge_key) in topo.dual_adj[u]
                # CRITICAL: Stop at cut edges!
                # These edges act as the "seams" of our UV map.
                if edge_key in cut_edges
                    continue
                end
                
                if !visited[v]
                    # Retrieve solver's period jump 'p' (defined for min < max)
                    p = get(field.period_jumps, edge_key, 0)
                    
                    # --- DERIVATION OF SIGN LOGIC ---
                    # The solver ensures: theta_j - theta_i approx p * (pi/2)
                    # We want global angles Theta to match: Theta_j approx Theta_i
                    # Where Theta = theta + r * (pi/2)
                    #
                    # Substitution:
                    # (theta_j + r_j*pi/2) - (theta_i + r_i*pi/2) = 0
                    # (theta_j - theta_i) + (r_j - r_i)*pi/2 = 0
                    # p*pi/2 + (r_j - r_i)*pi/2 = 0
                    # p + r_j - r_i = 0  =>  r_j = r_i - p
                    #
                    # This derivation assumes direction i -> j matches the solver's 'p' direction.
                    # The solver stores 'p' for the edge (min, max).
                    
                    if u < v
                        # Moving in direction of stored p: r_v = r_u - p
                        r[v] = r[u] - p
                    else
                        # Moving against direction of stored p: r_v = r_u + p
                        r[v] = r[u] + p
                    end
                    
                    visited[v] = true
                    enqueue!(q, v)
                end
            end
        end
    end
    
    # Compute the transitions (k) across the cut edges
    # These are needed if you want to stitch the UV map later, 
    # though the Least Squares solver implicitly handles boundaries by not adding equations for them.
    cut_rotations = Dict{Tuple{Int,Int}, Int}()
    
    for edge in cut_edges
        u, v = (edge[1] < edge[2]) ? (edge[1], edge[2]) : (edge[2], edge[1])
        
        if visited[u] && visited[v]
            p = get(field.period_jumps, edge, 0)
            
            # Calculate the rotation mismatch k
            # If the field was perfectly continuous (no singularity), k would be 0.
            # k != 0 indicates a singularity is enclosed by the loop formed by this cut.
            k = r[v] - r[u] + p
            cut_rotations[edge] = k
        end
    end
    
    println("  Propagated orientation to $(count(visited)) / $n_faces faces.")
    if count(visited) < n_faces
        println("  WARNING: Some faces were unreachable (bad cut graph?).")
    end
    
    return r, cut_rotations
end


function plot_smooth_global_field(field::CrossField, rotations::Vector{Int}, filename::String; cut_edges=nothing, verbose=false)
    # Create Figure with single panel
    fig = Figure(size=(1000, 800))
    ax = Axis(fig[1, 1], aspect=DataAspect(), title="Globally Smoothed Cross Field (u and v)")
    
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

    # 3. Compute Smooth Vectors for both u and v
    xs, ys, us_u, vs_u, us_v, vs_v = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
    
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
        
        # D. Compute Global Smooth Angles
        # u direction: Reference + Local_Theta + (Rotation * pi/2)
        global_angle_u = ref_ang + local_theta + (r * π / 2.0)
        # v direction: perpendicular to u (add pi/2)
        global_angle_v = global_angle_u + π / 2.0
        
        # E. Centroid
        tri = field.topology.faces[i]
        c = (field.topology.vertices[tri[1]].x + field.topology.vertices[tri[2]].x + field.topology.vertices[tri[3]].x)/3
        cy = (field.topology.vertices[tri[1]].y + field.topology.vertices[tri[2]].y + field.topology.vertices[tri[3]].y)/3
        
        push!(xs, c); push!(ys, cy)
        push!(us_u, cos(global_angle_u) * 0.04) # Scale length as needed
        push!(vs_u, sin(global_angle_u) * 0.04)
        push!(us_v, cos(global_angle_v) * 0.04)
        push!(vs_v, sin(global_angle_v) * 0.04)
    end

    # Plot both u and v arrows on same axis
    arrows2d!(ax, xs, ys, us_u, vs_u, color=:blue, label="Global 'u' Field")
    arrows2d!(ax, xs, ys, us_v, vs_v, color=:green, label="Global 'v' Field")
    
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
# 6. MAIN
# ==============================================================================

# ==============================================================================
# 8. PARAMETERIZATION & EXTRACTION
# ==============================================================================

"""
Solve for global (u, v) parameterization using Least Squares.
Minimize energy: sum_edges ||(u_j - u_i) - target_du||^2
"""
function compute_parameterization_least_squares(field::CrossField, rotations::Vector{Int}, cut_edges::Set{Tuple{Int,Int}})
    println("\n--- Solving Global Parameterization (Boundary Locked) ---")
    
    topo = field.topology
    n_verts = length(topo.vertices)
    n_faces = length(topo.faces)

    # 1. Identify Boundary Edges
    # An edge is a boundary if it only has one incident face in dual_adj
    boundary_edges = Set{Tuple{Int,Int}}()
    edge_counts = Dict{Tuple{Int,Int}, Int}()
    
    for tri in topo.faces
        edges = [(tri[1], tri[2]), (tri[2], tri[3]), (tri[3], tri[1])]
        for (u, v) in edges
            key = minmax(u, v)
            edge_counts[key] = get(edge_counts, key, 0) + 1
        end
    end
    for (key, count) in edge_counts
        if count == 1
            push!(boundary_edges, key)
        end
    end
    println("  identified $(length(boundary_edges)) boundary edges to lock.")

    # 2. Build Linear System
    I, J, V = Int[], Int[], Float64[]
    b_u, b_v = Float64[], Float64[]
    row_idx = 0
    
    seen_edges = Set{Tuple{Int,Int}}()

    # Pre-calculate face basis
    face_basis_u = Vector{Tuple{Float64, Float64}}(undef, n_faces)
    face_basis_v = Vector{Tuple{Float64, Float64}}(undef, n_faces)

    for i in 1:n_faces
        re = field.topology.face_ref_edges[i]
        p1 = field.topology.vertices[re[1]]; p2 = field.topology.vertices[re[2]]
        ref_ang = atan(p2.y - p1.y, p2.x - p1.x)
        global_angle = ref_ang + field.theta[i] + (rotations[i] * π / 2.0)
        face_basis_u[i] = (cos(global_angle), sin(global_angle))
        face_basis_v[i] = (-sin(global_angle), cos(global_angle))
    end

    for i in 1:n_faces
        tri = topo.faces[i]
        edges = [(tri[1], tri[2]), (tri[2], tri[3]), (tri[3], tri[1])]

        for (u, v) in edges
            key = minmax(u, v)
            if key in seen_edges; continue; end
            push!(seen_edges, key)
            if key in cut_edges; continue; end
            
            row_idx += 1
            
            # Geometry
            p_u_vert = topo.vertices[u]; p_v_vert = topo.vertices[v]
            dx = p_v_vert.x - p_u_vert.x; dy = p_v_vert.y - p_u_vert.y

            # Project edge onto local field
            dir_u, dir_v = face_basis_u[i], face_basis_v[i]
            target_du = dx * dir_u[1] + dy * dir_u[2]
            target_dv = dx * dir_v[1] + dy * dir_v[2]

            # --- BOUNDARY LOCKING LOGIC ---
            # If this is a boundary edge, we force it to align strictly 
            # with EITHER u or v, whichever it is closer to.
            weight_u, weight_v = 1.0, 1.0

            if key in boundary_edges
                # Check alignment
                # Dot product of Edge Vector (dx, dy) with U-basis (dir_u)
                # Normalize geometry for fair comparison
                len = sqrt(dx^2 + dy^2)
                dot_u = abs(target_du / len)
                dot_v = abs(target_dv / len)

                STIFFNESS = 1000.0  # How strongly to force the alignment
                
                if dot_u > dot_v
                    # Edge is mostly parallel to U. 
                    # This means it should NOT change in V. Force dV = 0.
                    target_dv = 0.0
                    weight_v = STIFFNESS
                else
                    # Edge is mostly parallel to V.
                    # This means it should NOT change in U. Force dU = 0.
                    target_du = 0.0
                    weight_u = STIFFNESS
                end
            end
            # ------------------------------

            # Add Equations with Weights
            push!(I, row_idx); push!(J, v); push!(V, 1.0 * weight_u)
            push!(I, row_idx); push!(J, u); push!(V, -1.0 * weight_u)
            push!(b_u, target_du * weight_u)
            
            # Note: We use the same row index for b_v, but we need separate matrix lines 
            # if we solved them as one big 2N x 2N system. 
            # Since we solve U and V separately with the same topology, 
            # we need distinct V values for the V-system if weights differ.
            # However, standard sparse solvers take A and b. 
            # To handle different weights for U and V efficiently in this simple setup,
            # we will just accept the higher weight constraint into both or solve purely coupled.
            
            # SIMPLIFICATION for separate solves:
            # We construct two separate A matrices effectively.
            # Actually, let's just push to b_v here and handle the matrix split below.
            push!(b_v, target_dv * weight_v)
        end
    end

    # NOTE: Because we applied different weights to U and V equations based on alignment,
    # We can no longer use a single matrix 'A' for both. We must build Au and Av.
    
    # Re-loop to build separate matrices (Cleanest way)
    # (Or just store weights separately and construct sparse twice)
    
    println("  Building stiffened systems...")
    I_u, J_u, V_u = Int[], Int[], Float64[]
    I_v, J_v, V_v = Int[], Int[], Float64[]
    
    # Anchor vertex (e.g. 1) to remove translation drift
    # Force vertex 1 to (0,0)
    # Note: For airfoils, pinning a vertex on the LEADING EDGE is better than random.
    # But vertex 1 is usually okay.
    anchor_idx = 1
    
    # Re-iterate to fill matrices correctly
    row = 0
    seen_edges = Set{Tuple{Int,Int}}() # reset
    for i in 1:n_faces
        tri = topo.faces[i]
        edges = [(tri[1], tri[2]), (tri[2], tri[3]), (tri[3], tri[1])]
        for (u, v) in edges
            key = minmax(u, v)
            if key in seen_edges || key in cut_edges; continue; end
            push!(seen_edges, key)
            row += 1
            
            # ... (re-calc weights from logic above) ...
            p_u_vert = topo.vertices[u]; p_v_vert = topo.vertices[v]
            dx = p_v_vert.x - p_u_vert.x; dy = p_v_vert.y - p_u_vert.y
            re = field.topology.face_ref_edges[i]
            p1 = field.topology.vertices[re[1]]; p2 = field.topology.vertices[re[2]]
            ref_ang = atan(p2.y - p1.y, p2.x - p1.x)
            global_angle = ref_ang + field.theta[i] + (rotations[i] * π / 2.0)
            dir_u = (cos(global_angle), sin(global_angle))
            dir_v = (-sin(global_angle), cos(global_angle))

            target_du = dx * dir_u[1] + dy * dir_u[2]
            target_dv = dx * dir_v[1] + dy * dir_v[2]
            
            w_u, w_v = 1.0, 1.0
            if key in boundary_edges
                len = sqrt(dx^2 + dy^2)
                dot_u = abs(target_du / len)
                dot_v = abs(target_dv / len)
                STIFFNESS = 1000.0
                if dot_u > dot_v; target_dv = 0.0; w_v = STIFFNESS; else; target_du = 0.0; w_u = STIFFNESS; end
            end

            push!(I_u, row); push!(J_u, v); push!(V_u, w_u); push!(I_u, row); push!(J_u, u); push!(V_u, -w_u)
            b_u[row] = target_du * w_u # update b with weight
            
            push!(I_v, row); push!(J_v, v); push!(V_v, w_v); push!(I_v, row); push!(J_v, u); push!(V_v, -w_v)
            b_v[row] = target_dv * w_v # update b with weight
        end
    end
    
    # Add anchor
    row += 1
    push!(I_u, row); push!(J_u, anchor_idx); push!(V_u, 1.0); push!(b_u, 0.0)
    push!(I_v, row); push!(J_v, anchor_idx); push!(V_v, 1.0); push!(b_v, 0.0)

    # Solve U
    Au = sparse(I_u, J_u, V_u, row, n_verts)
    u_sol = Au \ b_u
    
    # Solve V
    Av = sparse(I_v, J_v, V_v, row, n_verts)
    v_sol = Av \ b_v
    
    println("  Range U: [$(minimum(u_sol)), $(maximum(u_sol))]")
    println("  Range V: [$(minimum(v_sol)), $(maximum(v_sol))]")
    return u_sol, v_sol
end

"""
Extract the quad mesh by tracing isolines of the scalar fields u and v.
"""
function extract_quad_mesh(topo::MeshTopology, u_coords::Vector{Float64}, v_coords::Vector{Float64}, grid_size_x::Float64, grid_size_y::Float64)
    println("\n--- Extracting Quad Mesh ---")
    quads = Vector{Tuple{Point3D, Point3D, Point3D, Point3D}}()
    
    # Normalize coords by grid size
    u_scaled = u_coords ./ grid_size_x
    v_scaled = v_coords ./ grid_size_y

    # This is a visualization-only extraction (finding integer crossings)
    # For a topological extraction (generating a new MeshTopology), you need to cut triangles.
    # Here we generate line segments for visualization.
    
    u_lines = Vector{Tuple{Point3D, Point3D}}()
    v_lines = Vector{Tuple{Point3D, Point3D}}()

    for i in 1:length(topo.faces)
        tri = topo.faces[i]
        ids = [tri[1], tri[2], tri[3]]
        
        # Check for U integer crossings
        u_vals = [u_scaled[id] for id in ids]
        v_vals = [v_scaled[id] for id in ids]
        
        # Find min/max integer range in this triangle
        u_min, u_max = floor(Int, minimum(u_vals)), ceil(Int, maximum(u_vals))
        v_min, v_max = floor(Int, minimum(v_vals)), ceil(Int, maximum(v_vals))

        # Trace U-isolines (where u = constant integer)
        for u_iso in u_min:u_max
            # Find intersection points on edges
            pts = Point3D[]
            for k in 1:3
                p1_idx, p2_idx = ids[k], ids[mod1(k+1,3)]
                u1, u2 = u_vals[k], u_vals[mod1(k+1,3)]
                
                # Check if isoline crosses edge
                if (u1 <= u_iso && u2 >= u_iso) || (u1 >= u_iso && u2 <= u_iso)
                    if abs(u2 - u1) > 1e-6
                        t = (u_iso - u1) / (u2 - u1)
                        p1, p2 = topo.vertices[p1_idx], topo.vertices[p2_idx]
                        px = p1.x + t * (p2.x - p1.x)
                        py = p1.y + t * (p2.y - p1.y)
                        push!(pts, Point3D(px, py, 0.0))
                    end
                end
            end
            if length(pts) == 2
                push!(u_lines, (pts[1], pts[2]))
            end
        end

        # Trace V-isolines (where v = constant integer)
        for v_iso in v_min:v_max
            pts = Point3D[]
            for k in 1:3
                p1_idx, p2_idx = ids[k], ids[mod1(k+1,3)]
                v1, v2 = v_vals[k], v_vals[mod1(k+1,3)]
                
                if (v1 <= v_iso && v2 >= v_iso) || (v1 >= v_iso && v2 <= v_iso)
                    if abs(v2 - v1) > 1e-6
                        t = (v_iso - v1) / (v2 - v1)
                        p1, p2 = topo.vertices[p1_idx], topo.vertices[p2_idx]
                        px = p1.x + t * (p2.x - p1.x)
                        py = p1.y + t * (p2.y - p1.y)
                        push!(pts, Point3D(px, py, 0.0))
                    end
                end
            end
            if length(pts) == 2
                push!(v_lines, (pts[1], pts[2]))
            end
        end
    end
    
    return u_lines, v_lines
end

function plot_quad_extraction(field, u_lines, v_lines, filename)
    fig = Figure(size=(1000, 800))
    ax = Axis(fig[1, 1], aspect=DataAspect(), title="Extracted Quad Grid")
    
    # Plot background mesh (faint)
    for tri in field.topology.faces
        pts = [field.topology.vertices[i] for i in tri]; push!(pts, pts[1])
        lines!(ax, [p.x for p in pts], [p.y for p in pts], color=(:gray, 0.1), linewidth=0.2)
    end

    # Plot U lines (Blue)
    for seg in u_lines
        lines!(ax, [seg[1].x, seg[2].x], [seg[1].y, seg[2].y], color=:blue, linewidth=1.5)
    end

    # Plot V lines (Red)
    for seg in v_lines
        lines!(ax, [seg[1].x, seg[2].x], [seg[1].y, seg[2].y], color=:red, linewidth=1.5)
    end

    save(filename, fig)
    println("Saved quad extraction to $filename")
end

function main()
    # Replace with your file
    # filename = "triangulations/mesh_airfoil_dae11.su2" 
    # filename = "triangulations/regular-square-8x8.msh"
    # filename = "triangulations/disk-radial-fine.msh"
    # filename = "triangulations/crmhl_test.su2"
    # filename = "triangulations/two_element.msh"
    filename = "triangulations/back-step.msh"
    
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
    # constraints = Dict{Int, Float64}(1 => 0.0)  # Dummy constraint on face 1

    # Set verbosity and animation parameters
    verbose = false
    animate = false  # Set to true to generate animation
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
    # Visualize
    println("Generating visualization...")
    plot_results(field, "miq_solution.png"; verbose=verbose, cut_edges=cut_edges)
    println("Complete! Visualization saved to miq_solution.png")
u_sol, v_sol = compute_parameterization_least_squares(field, face_rotations, cut_edges)

    # 9. Extract and Visualize Grid
    # Define grid sizing (target edge length). 
    # Calculate approximate bounding box diagonal to guess a good size if unknown
    min_x = minimum(v.x for v in topo.vertices)
    max_x = maximum(v.x for v in topo.vertices)
    grid_size_x = (max_x - min_x) / 100.0 # Aim for ~30 quads across width

    min_y = minimum(v.y for v in topo.vertices)
    max_y = maximum(v.y for v in topo.vertices)
    grid_size_y = (max_y - min_y) / 100
    
    println("Extracting grid with size: ($grid_size_x, $grid_size_y)")
    u_lines, v_lines = extract_quad_mesh(topo, u_sol, v_sol, grid_size_x, grid_size_y)
    
    plot_quad_extraction(field, u_lines, v_lines, "quad_extraction.png")
    
    println("Process Complete.")
end



main()