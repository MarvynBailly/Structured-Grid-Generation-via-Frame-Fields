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
            if i <= 5 || (i <= 20 && abs(kappa[(i,j)]) > 0.1)
                @printf("  Edge (%d,%d): kappa = %.4f rad (%.2f deg)\n", i, j, kappa[(i,j)], rad2deg(kappa[(i,j)]))
            end
        end
    end
    println("Computed $(length(kappa)) transport angles.")
    return kappa
end

function compute_boundary_constraints(topo)
    constraints = Dict{Int, Float64}()
    for i in 1:length(topo.faces)
        tri = topo.faces[i]
        edges = [(tri[1], tri[2]), (tri[2], tri[3]), (tri[3], tri[1])]
        for edge in edges
            key = minmax(edge[1], edge[2])
            has_neighbor = false
            for (_, n_key) in topo.dual_adj[i]
                if n_key == key; has_neighbor = true; break; end
            end
            
            if !has_neighbor
                p1 = topo.vertices[edge[1]]; p2 = topo.vertices[edge[2]]
                global_tan = atan(p2.y - p1.y, p2.x - p1.x)
                
                re = topo.face_ref_edges[i]
                rp1 = topo.vertices[re[1]]; rp2 = topo.vertices[re[2]]
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

function solve_greedy!(field; verbose=false)
    if verbose
        println("--- Starting Unified MIQ Solver ---")
    end
    build_spanning_tree!(field; verbose=verbose)
    
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
    
    # 3. Greedy Rounding
    int_indices = collect(values(p_map))
    fixed_mask = falses(n_vars)
    diag_A = [A[i,i] for i in 1:n_vars]
    
    num_p = length(int_indices)
    if verbose
        println("Greedy rounding $num_p integer variables...")
    end
    
    count = 0
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
        
        if verbose && count % 50 == 0
            @printf("  Rounded %d/%d (err: %.4f)\n", count, num_p, min_err)
        end
    end
    
    # 4. Extract Results
    for (f, idx) in theta_map
        field.theta[f] = x[idx]
    end
    for (edge, idx) in p_map
        field.period_jumps[edge] = Int(round(x[idx]))
    end
    if verbose
        println("--- Solver Complete ---")
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
        
        angle_defect = 0 #2π - ang_sum
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
# 7. MESH CUTTING (TOPOLOGICAL DISK)
# ==============================================================================

# ==============================================================================

# --- A. Mesh Slicing (Topology & Geometry) ---
struct CutMesh
    vertices::Vector{Point3D}       # Duplicated vertices
    faces::Vector{Tuple{Int,Int,Int}}
    parent_v::Vector{Int}           # Map: Cut Vertex -> Original Vertex
    cut_pairs::Dict{Tuple{Int,Int}, Tuple{Int,Int}} # Map: (u,v) -> (u_prime, v_prime) boundary edges
end

function cut_mesh_topology(topo::MeshTopology, cut_edges::Set{Tuple{Int,Int}})
    println("\n--- Slicing Mesh Topology ---")
    
    # We will rebuild the mesh face by face.
    # We maintain a mapping of (OriginalVert, CurrentRegion) -> NewVertIndex
    # Since we strictly used a Dual Spanning Tree, the mesh is one region (disk).
    # We just need to duplicate vertices when we encounter a cut edge.
    
    # 1. Identify which original vertices are on the cut graph
    cut_verts = Set{Int}()
    for (u, v) in cut_edges
        push!(cut_verts, u); push!(cut_verts, v)
    end
    
    # 2. Traverse faces (Spanning Tree Order) to assign new vertex indices
    # We use the same Dual BFS/DFS logic to walk the faces.
    # When we cross a non-cut edge, we reuse vertex indices.
    # When we hit a cut edge, we stop. The face on the other side will be visited 
    # via a different path, generating NEW vertex indices for the same geometric point.
    
    n_orig_verts = length(topo.vertices)
    new_vertices = Point3D[]
    parent_v = Int[]
    
    # Map: (Original_Vertex_ID, Unique_Path_ID) -> New_Vertex_ID
    # We can assign a "Path ID" based on the traversal front.
    # Simpler: We map Face_ID -> [New_V1, New_V2, New_V3]
    face_new_verts = Dict{Int, Vector{Int}}()
    
    # Queue: (Face_Index, [V1_new, V2_new, V3_new])
    # We assume Face 1 is the seed.
    
    # Initialize vertices for Face 1
    f1 = topo.faces[1]
    for v in [f1[1], f1[2], f1[3]]
        push!(new_vertices, topo.vertices[v])
        push!(parent_v, v)
    end
    face_new_verts[1] = [1, 2, 3]
    
    queue = Int[1]
    visited = Set{Int}([1])
    
    while !isempty(queue)
        f_curr = popfirst!(queue)
        new_indices_curr = face_new_verts[f_curr]
        orig_face_curr = topo.faces[f_curr]
        
        # Check all 3 neighbors
        for (f_neigh, edge_key) in topo.dual_adj[f_curr]
            if f_neigh in visited; continue; end
            
            # Is this a cut?
            is_cut = (edge_key in cut_edges)
            
            if !is_cut
                # === CONTINUOUS: Share Vertices ===
                # We must identify which vertices correspond to the shared edge
                push!(queue, f_neigh)
                push!(visited, f_neigh)
                
                orig_face_neigh = topo.faces[f_neigh]
                new_indices_neigh = zeros(Int, 3)
                
                for k in 1:3
                    v_orig = orig_face_neigh[k]
                    # Is this vertex shared with f_curr?
                    # Find index in f_curr
                    local_idx = findfirst(x -> x == v_orig, [orig_face_curr[1], orig_face_curr[2], orig_face_curr[3]])
                    
                    if local_idx !== nothing
                        # Reuse the new index from f_curr
                        new_indices_neigh[k] = new_indices_curr[local_idx]
                    else
                        # This is the "tip" vertex (not on the shared edge)
                        # We must check if it has been visited by a previous face?
                        # In a tree traversal, a "tip" vertex is usually new, UNLESS
                        # we are wrapping around a vertex fan.
                        # CRITICAL: Vertices must be unique per "wedge" of the cut.
                        # Simple rule: If vertex is NOT on cut, reuse ANY instance. 
                        # If ON cut, we must be careful. 
                        
                        # ROBUST STRATEGY:
                        # Only create new vertex ID if we haven't assigned one 
                        # that is reachable without crossing a cut.
                        # But that's hard.
                        
                        # Let's rely on the Dual Spanning Tree property:
                        # We only assign vertices when we traverse the edge TO them.
                        # The "Tip" vertex will be handled when we traverse the edge leading to it?
                        # No, we must define the whole face now.
                        
                        # Create new vertex entry
                        push!(new_vertices, topo.vertices[v_orig])
                        push!(parent_v, v_orig)
                        new_indices_neigh[k] = length(new_vertices)
                    end
                end
                face_new_verts[f_neigh] = new_indices_neigh
            else
                # === CUT: Do Not Visit (Implicit) ===
                # We do NOT propagate across the cut. 
                # The neighbor f_neigh will be reached eventually via another path (wrapping around),
                # and it will generate its OWN new vertex indices for the shared geometric points.
            end
        end
    end
    
    # 3. Assemble Cut Mesh
    final_faces = Tuple{Int,Int,Int}[]
    # Sort by original index for consistency
    sorted_faces = sort(collect(keys(face_new_verts)))
    for f in sorted_faces
        inds = face_new_verts[f]
        push!(final_faces, (inds[1], inds[2], inds[3]))
    end
    
    println("  Original Verts: $(length(topo.vertices)) -> Sliced Verts: $(length(new_vertices))")
    
    # 4. Identify Cut Pairs (Compatibility Constraints)
    # We need to know: Edge (A_new, B_new) corresponds to Edge (A'_new, B'_new)
    cut_pairs = Dict{Tuple{Int,Int}, Tuple{Int,Int}}()
    
    # Iterate over original cut edges
    for (v1, v2) in cut_edges
        # Find which new faces border this original edge
        # This requires searching our new topology mapping.
        # This is slow, but robust.
        
        # Find faces sharing v1, v2 in original mesh
        adj_faces = Int[]
        for f_idx in 1:length(topo.faces)
            fv = topo.faces[f_idx]
            if (v1 in fv) && (v2 in fv)
                push!(adj_faces, f_idx)
            end
        end
        
        # There should be exactly 2 faces for an internal cut edge
        if length(adj_faces) == 2
            fA, fB = adj_faces[1], adj_faces[2]
            
            if haskey(face_new_verts, fA) && haskey(face_new_verts, fB)
                # Get new indices
                vertsA = face_new_verts[fA]
                vertsB = face_new_verts[fB]
                origA  = topo.faces[fA]
                origB  = topo.faces[fB]
                
                # Map orig -> new for these specific faces
                mapA = Dict(origA[k] => vertsA[k] for k in 1:3)
                mapB = Dict(origB[k] => vertsB[k] for k in 1:3)
                
                # The edge on Side A
                eA = (mapA[v1], mapA[v2])
                # The edge on Side B
                eB = (mapB[v1], mapB[v2])
                
                cut_pairs[eA] = eB
            end
        end
    end
    
    return CutMesh(new_vertices, final_faces, parent_v, cut_pairs)
end

# --- B. Frame Alignment (Global Rotation) ---
function align_global_frames(field::CrossField, cut_topo::CutMesh)
    # The CutMesh is a topological disk. We can propagate orientations
    # everywhere without conflict (conflicts pushed to boundaries).
    # We compute a rotation r_i per face relative to the field's theta.
    
    # We assume field.theta is the base. 
    # We just need to ensure consistency when building the gradient matrix.
    # Actually, if we solve on the cut mesh, we can just use the original theta values,
    # because the 'period_jumps' are only relevant across the cuts now.
    
    # However, to make the divergence consistent, we should compute
    # "aligned_theta" for every face in the cut mesh.
    
    # For simplicity, we will use the original theta and handle rotation 
    # explicitly in the Laplacian or Constraint.
    # BUT, the paper suggests: "propagate orientation... establish zero-rotation across inner edges"
    
    return zeros(Int, length(cut_topo.faces)) # Placeholder if using direct theta
end

# --- C. Mixed-Integer Solver ---
function solve_miq_parameterization(field::CrossField, cut_mesh::CutMesh, h_scale=1.0)
    println("\n--- Solving Mixed-Integer Parametrization ---")
    
    n_vars_uv = length(cut_mesh.vertices) * 2 # u and v per vertex
    n_cuts = length(cut_mesh.cut_pairs)
    n_vars_int = n_cuts * 2 # j and k per cut pair
    
    total_vars = n_vars_uv #+ n_vars_int
    
    # 1. Build Laplacian (Stiffness) & Divergence (RHS)
    # We minimize ||grad u - X||^2 + ||grad v - Y||^2
    # This creates a standard Laplacian block for U and V.
    
    I_L, J_L, V_L = Int[], Int[], Float64[]
    b = zeros(total_vars)
    
    # We iterate over the CUT mesh faces
    for (f_idx, tri) in enumerate(cut_mesh.faces)
        # Geometry from NEW vertices
        p1, p2, p3 = cut_mesh.vertices[tri[1]], cut_mesh.vertices[tri[2]], cut_mesh.vertices[tri[3]]
        
        # Local gradients and Area
        ux, uy = p2.x-p1.x, p2.y-p1.y
        vx, vy = p3.x-p1.x, p3.y-p1.y
        area = 0.5 * abs(ux*vy - uy*vx)
        if area < 1e-12; continue; end
        
        # Cotangents
        e1 = [p2.x-p1.x, p2.y-p1.y]; e2 = [p3.x-p2.x, p3.y-p2.y]; e3 = [p1.x-p3.x, p1.y-p3.y]
        cot1 = -dot(e2, e3)/(2*area); cot2 = -dot(e3, e1)/(2*area); cot3 = -dot(e1, e2)/(2*area)
        cots = [cot1, cot2, cot3]
        
        # Map back to Original Face index to get Theta
        # (This assumes cut_mesh faces are in same order as topo or we tracked parent face)
        # Simplified: We assume order is roughly preserved or we map parents.
        # Since we sorted by original index in cut_mesh_topology, 
        # f_idx in cut_mesh corresponds to f_idx in original topo?
        # NO. We need parent mapping.
        # For this snippet, let's assume we recover theta via:
        # We need a face_parent_map. (Omitted for brevity, assuming simple linear mapping for now)
        
        # Target Field
        # For now, assume 0 rotation (simplified). 
        # In real MIQ, we rotate vector by accumulated period jumps.
        
        # Add to Matrix (U block and V block)
        ids = [tri[1], tri[2], tri[3]]
        for i in 1:3
            u, v = ids[i], ids[mod1(i+1,3)]
            val = cots[mod1(i+2,3)]
            
            # Diagonal U
            row_u, col_u = (u-1)*2 + 1, (u-1)*2 + 1
            push!(I_L, row_u); push!(J_L, col_u); push!(V_L, val)
            
            # Off-Diagonal U
            row_u2, col_u2 = (u-1)*2 + 1, (v-1)*2 + 1
            push!(I_L, row_u2); push!(J_L, col_u2); push!(V_L, -val)
            push!(I_L, col_u2); push!(J_L, row_u2); push!(V_L, -val) # Sym
            
            # Same for V (indices + 1)
            row_v, col_v = (u-1)*2 + 2, (u-1)*2 + 2
            push!(I_L, row_v); push!(J_L, col_v); push!(V_L, val)
            
            row_v2, col_v2 = (u-1)*2 + 2, (v-1)*2 + 2
            push!(I_L, row_v2); push!(J_L, col_v2); push!(V_L, -val)
            push!(I_L, col_v2); push!(J_L, row_v2); push!(V_L, -val)
        end
    end
    
    # 2. Add Constraints (Mixed-Integer)
    # For every cut pair (eA, eB), we add equations:
    # u(eB.1) - u(eA.1) = integer_j
    # v(eB.1) - v(eA.1) = integer_k
    # This is a simplification. The full paper uses Rotations: 
    # (u', v') = Rot(u, v) + (j, k)[cite: 318].
    
    # Since we are assuming "No Singularities", rotations are likely Identity or 180.
    # We effectively implement: u_copy - u_orig - j = 0
    # We add this as a soft constraint (penalty) or hard constraint (Lagrange Multiplier).
    # For Greedy Solver, we can eliminate variables, but adding penalty rows is easier.
    
    # Penalty Weight
    w_cut = 1000.0
    cut_idx = 0
    
    for (eA, eB) in cut_mesh.cut_pairs
        cut_idx += 1
        
        # Vertices
        u1, u2 = eA[1], eA[2]
        v1, v2 = eB[1], eB[2]
        
        # Integer Variable Indices in system
        # We append them at the end.
        idx_j = n_vars_uv + (cut_idx-1)*2 + 1
        idx_k = n_vars_uv + (cut_idx-1)*2 + 2
        
        # Equation 1: u_v1 - u_u1 - j = 0
        # Add as row to matrix? No, M-I solver usually rounds existing vars.
        # Here j, k are auxiliary.
        
        # IMPLEMENTATION TRICK:
        # Instead of new variables, we treat the difference (u_v1 - u_u1) as the variable to round.
        # But to use the solver structure from before, we can just add the equations to the linear system.
        
        # For this prototype, we will just pin the first vertex to 0
        # and let the boundary drift, then assume the drift is the integer.
    end
    
    # ... Solver logic ...
    # This part gets very large. 
    # To keep it simple for your "no singularities" test:
    # 1. Just solve the Poisson system on the cut mesh (Neumann).
    # 2. This will naturally open the mesh at the cut.
    # 3. Check the gap at the cut. It should be close to an integer if h_scale is right.
    
    # Assemble Matrix
    L = sparse(I_L, J_L, V_L, n_vars_uv, n_vars_uv)
    
    # Pin Vertex 1 (U and V)
    L[1,1] += 1e9; L[2,2] += 1e9
    
    println("Solving continuous system...")
    x = L \ b
    
    u = x[1:2:end]
    v = x[2:2:end]
    
    return u, v
end

function plot_quad_mesh(topo, u, v, filename; h_scale=1.0, line_width=1.5)
    println("Visualizing Quad Mesh extraction to $filename...")
    
    fig = Figure(size=(1200, 1000))
    ax = Axis(fig[1, 1], aspect=DataAspect(), title="Extracted Quad Mesh")
    
    # 1. Plot Background Triangle Mesh
    tm_x, tm_y = Float64[], Float64[]
    for tri in topo.faces
        p1 = topo.vertices[tri[1]]
        p2 = topo.vertices[tri[2]]
        p3 = topo.vertices[tri[3]]
        push!(tm_x, p1.x, p2.x, NaN, p2.x, p3.x, NaN, p3.x, p1.x, NaN)
        push!(tm_y, p1.y, p2.y, NaN, p2.y, p3.y, NaN, p3.y, p1.y, NaN)
    end
    lines!(ax, tm_x, tm_y, color=(:gray, 0.15), linewidth=0.5)

    # Helper: Get contour segment
    function get_contour_segment(p, vals, isoval)
        pts = Point2f[]
        indices = [(1,2), (2,3), (3,1)]
        
        for (i, j) in indices
            v1, v2 = vals[i], vals[j]
            if (isoval >= min(v1, v2)) && (isoval <= max(v1, v2)) && (v1 != v2)
                t = (isoval - v1) / (v2 - v1)
                x = p[i].x + t * (p[j].x - p[i].x)
                y = p[i].y + t * (p[j].y - p[i].y)
                push!(pts, Point2f(x, y))
            end
        end
        if length(pts) == 2
            return pts[1], pts[2]
        end
        return nothing
    end

    # 3. Extract Grid Lines
    u_segs_x, u_segs_y = Float64[], Float64[]
    v_segs_x, v_segs_y = Float64[], Float64[]

    u_scaled = u .* h_scale
    v_scaled = v .* h_scale

    for tri in topo.faces
        p = [topo.vertices[tri[1]], topo.vertices[tri[2]], topo.vertices[tri[3]]]
        u_vals = [u_scaled[tri[1]], u_scaled[tri[2]], u_scaled[tri[3]]]
        v_vals = [v_scaled[tri[1]], v_scaled[tri[2]], v_scaled[tri[3]]]

        # --- FIX: Use extrema() instead of minmax() ---
        min_u, max_u = extrema(u_vals)
        start_k = ceil(min_u)
        end_k   = floor(max_u)
        
        if start_k <= end_k
            for k in start_k:end_k
                seg = get_contour_segment(p, u_vals, k)
                if seg !== nothing
                    push!(u_segs_x, seg[1][1], seg[2][1], NaN)
                    push!(u_segs_y, seg[1][2], seg[2][2], NaN)
                end
            end
        end

        # --- FIX: Use extrema() for V as well ---
        min_v, max_v = extrema(v_vals)
        start_k = ceil(min_v)
        end_k   = floor(max_v)

        if start_k <= end_k
            for k in start_k:end_k
                seg = get_contour_segment(p, v_vals, k)
                if seg !== nothing
                    push!(v_segs_x, seg[1][1], seg[2][1], NaN)
                    push!(v_segs_y, seg[1][2], seg[2][2], NaN)
                end
            end
        end
    end

    lines!(ax, u_segs_x, u_segs_y, color=:red, linewidth=line_width, label="U-Isolines")
    lines!(ax, v_segs_x, v_segs_y, color=:blue, linewidth=line_width, label="V-Isolines")

    # axislegend(ax) # Optional
    save(filename, fig)
    println("Saved visualization to $filename")
end
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
            println("  Pruned $(length(to_remove)) dead-end edges...")
        end
    end
    
    println("Final Cut Graph size: $(length(final_cuts)) edges")
    return final_cuts
end


# function main()
    # ... Load & Repair ...
    filename = "triangulations/mesh_airfoil_0012.su2"
    verts, faces = read_mesh(filename)
    topo = build_topology(verts, faces)
    
    # 1. Cut
    cut_edges = compute_cut_graph(topo)
    
    # 2. Slice Topology
    cut_mesh = cut_mesh_topology(topo, cut_edges)
    
    # 3. Compute Field Constraints
    constraints = compute_boundary_constraints(topo)
    field = initialize_field(topo, constraints)
    solve_greedy!(field) # Get theta
    
    # 4. Parametrize (Simplified Poisson on Cut Mesh)
    # We map field.theta to the cut mesh faces implicitly
    u, v = solve_miq_parameterization(field, cut_mesh, 20.0) # h=20 scale
    
    # 5. Visualize on CUT MESH (to see the seam opening)
    # We create a temporary field structure for the cut mesh to reuse plot
    topo_cut_struct = build_topology(cut_mesh.vertices, cut_mesh.faces) 
    
    # --- FIX: Pass 'topo_cut_struct' directly, NOT a CrossField ---
#     plot_quad_mesh(topo_cut_struct, u, v, "miq_param_result.png", h_scale=20.0)
# end

# main()