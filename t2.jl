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
            weight = 1.0 / (sqrt(len_sq) + 1e-6)
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

function solve_greedy!(field; verbose=true)
    println("--- Starting Unified MIQ Solver ---")
    build_spanning_tree!(field)
    
    # 1. Build System
    println("Assembling system...")
    A, b, theta_map, p_map = assemble_system_weighted(field)
    n_vars = length(b)
    
    # 2. Initial Solve
    println("Initial solve...")
    x = A \ b
    
    # 3. Greedy Rounding
    int_indices = collect(values(p_map))
    fixed_mask = falses(n_vars)
    diag_A = [A[i,i] for i in 1:n_vars]
    
    num_p = length(int_indices)
    println("Greedy rounding $num_p integer variables...")
    
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
    println("--- Solver Complete ---")
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

function compute_singularities(field::CrossField)
    println("\n--- Computing Singularities (Normalized) ---")
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
        
        println("\n=== Vertex $v_idx (star size: $n_star) ===")
        
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
            
            @printf("  Face %d: angle = %.4f rad (%.2f deg)\n", curr, angle_at_vertex, rad2deg(angle_at_vertex))
            
            # Transport
            k = haskey(field.transport_angles, (curr, next)) ? field.transport_angles[(curr,next)] : -field.transport_angles[(next,curr)]
            k_sum += k
            @printf("  Edge (%d->%d): kappa = %.4f rad (%.2f deg)\n", curr, next, k, rad2deg(k))
            
            # Period Jump (NORMALIZED)
            edge = minmax(curr, next)
            sign = (curr < next) ? 1.0 : -1.0
            
            if haskey(field.period_jumps, edge)
                raw_p = field.period_jumps[edge]
                
                # --- CHANGE IS HERE ---
                # Remove windings of 4 (360 degrees)
                # This maps ..., -4, 0, 4, ... -> 0
                # And ..., -5, -1, 3, ... -> -1
                norm_p = raw_p - 4 * round(raw_p / 4)
                @printf("  Edge %s: raw_p = %d -> norm_p = %d (sign = %.1f, contribution = %.2f)\n", 
                        edge, raw_p, Int(norm_p), sign, norm_p * sign)
                p_sum += norm_p * sign
            else
                println("  Edge $edge: no period jump")
            end
        end
        
        angle_defect = 2π - ang_sum
        I = (angle_defect + k_sum)/(2π) + p_sum/4.0
        
        @printf("  Angle sum: %.4f rad (%.2f deg)\n", ang_sum, rad2deg(ang_sum))
        @printf("  Angle defect: %.4f rad (%.2f deg)\n", angle_defect, rad2deg(angle_defect))
        @printf("  Kappa sum: %.4f rad (%.2f deg)\n", k_sum, rad2deg(k_sum))
        @printf("  Period sum: %.4f\n", p_sum)
        @printf("  Min angle: %.4f rad (%.2f deg)\n", min_angle, rad2deg(min_angle))
        @printf("  INDEX = (%.4f + %.4f)/(2π) + %.4f/4 = %.4f\n", angle_defect, k_sum, p_sum, I)
        
        # Check for degenerate geometry
        is_degenerate = (min_angle < MIN_ANGLE_THRESHOLD) || (ang_sum < MIN_ANGLE_SUM_THRESHOLD)
        
        if is_degenerate
            println("  ⚠️  DEGENERATE VERTEX SKIPPED (nearly-zero angles indicate bad mesh quality)")
            degenerate_vertices_skipped += 1
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
        end
    end
    
    # Print statistics
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
    if !isempty(all_indices)
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
    println("\n=== TOPOLOGICAL VERIFICATION ===")
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
    @printf("Mesh topology: V=%d, E=%d, F=%d\n", n_vertices, n_edges, n_faces)
    @printf("Euler characteristic χ = V - E + F = %d\n", euler_char)
    
    # Determine genus and expected index sum
    # For orientable surfaces: χ = 2 - 2g (g = genus)
    # Expected total index for cross fields (4-RoSy): (χ/2)
    genus = (2 - euler_char) / 2
    expected_index_sum = euler_char / 2.0
    
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
    
    if !isempty(singularities)
        actual_index_sum = sum(I for (_, I) in singularities)
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
        
        # Show individual singularities
        println("\nDetected singularities:")
        for (v_idx, I) in singularities
            p = topo.vertices[v_idx]
            @printf("  Vertex %d at (%.4f, %.4f): index = %.4f\n", v_idx, p.x, p.y, I)
        end
    else
        @printf("No singularities detected, but expected sum: %.4f\n", expected_index_sum)
        if abs(expected_index_sum) > 0.1
            println("⚠️  This suggests singularities are missing (possibly filtered as degenerate)")
        end
    end
    
    return singularities
end

function plot_results(field, filename)
    fig = Figure(size=(1000, 800))
    ax = Axis(fig[1, 1], aspect=DataAspect(), title="MIQ Solution (Constrained Faces & Singularities)")
    
    # Mesh
    for tri in field.topology.faces
        pts = [field.topology.vertices[i] for i in tri]; push!(pts, pts[1])
        lines!(ax, [p.x for p in pts], [p.y for p in pts], color=:gray90, linewidth=0.5)
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
    sings = compute_singularities(field)
    pos = Point2f[]; neg = Point2f[]
    for (v, I) in sings
        p = field.topology.vertices[v]
        pt = Point2f(p.x, p.y)
        if I > 0; push!(pos, pt); else; push!(neg, pt); end
    end
    scatter!(ax, pos, color=:red, markersize=15, strokecolor=:black, strokewidth=1, label="+ Index")
    scatter!(ax, neg, color=:cyan, markersize=15, strokecolor=:black, strokewidth=1, label="- Index")
    
    # Add constrained faces to legend
    if !isempty(field.constrained_faces)
        # Dummy plot for legend entry
        lines!(ax, [NaN], [NaN], color=:orange, linewidth=2, label="Constrained Faces")
    end
    
    if !isempty(pos) || !isempty(neg) || !isempty(field.constrained_faces)
        axislegend(ax)
    end
    
    save(filename, fig)
    println("Saved visualization to $filename")
    println("Found $(length(sings)) singularities.")
end

# ==============================================================================
# 6. MAIN
# ==============================================================================

function main()
    # Replace with your file
    filename = "triangulations/mesh_airfoil_dae11.su2" 
    
    if !isfile(filename)
        println("File not found: $filename")
        return
    end
    
    println("Reading mesh...")
    verts, faces = read_mesh(filename)
    topo = build_topology(verts, faces)
    
    println("Detecting boundaries...")
    constraints = compute_boundary_constraints(topo)
    println("Constrained $(length(constraints)) boundary faces.")
    
    field = initialize_field(topo, constraints)
    
    # Solve
    solve_greedy!(field)
        
    # Visualize
    plot_results(field, "miq_solution.png")
end

main()