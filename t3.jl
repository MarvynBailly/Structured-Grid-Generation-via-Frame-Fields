using LinearAlgebra
using SparseArrays
using Printf
using CairoMakie
using DataStructures 

# --- Data Structures ---
struct Point3D
    x::Float64; y::Float64; z::Float64
end

struct MeshTopology
    vertices::Vector{Point3D}
    faces::Vector{Tuple{Int,Int,Int}}
    dual_adj::Vector{Vector{Tuple{Int, Tuple{Int,Int}}}} 
    face_ref_edges::Vector{Tuple{Int,Int}}
    vert_to_faces::Vector{Vector{Int}} 
end

mutable struct CrossField
    topology::MeshTopology
    theta::Vector{Float64}
    period_jumps::Dict{Tuple{Int,Int}, Int}
    transport_angles::Dict{Tuple{Int,Int}, Float64}
    constrained_faces::Dict{Int, Float64}
    fixed_edges::Set{Tuple{Int,Int}}
end

# --- File I/O ---
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
                # Handle 2D or 3D lines
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
        elseif keyword == "NMARK"; break 
        else; i += 1; end
    end
    return vertices, faces
end

# --- Robust Mesh Repair (The Fix) ---
function safe_repair_mesh(vertices, faces)
    println("\n--- Starting Safe Mesh Repair ---")
    println("Initial: $(length(vertices)) vertices, $(length(faces)) faces")

    valid_faces = Tuple{Int,Int,Int}[]
    seen_faces = Set{Tuple{Int,Int,Int}}()
    
    skipped_area = 0
    skipped_dupe = 0
    
    for tri in faces
        # Sort indices to detect duplicates regardless of winding
        v_sorted = sort([tri[1], tri[2], tri[3]])
        key = (v_sorted[1], v_sorted[2], v_sorted[3])
        
        # Check A: Degenerate Indices (e.g. 1-1-2)
        if key[1] == key[2] || key[2] == key[3]; continue; end
        
        # Check B: Duplicate Face
        if key in seen_faces; skipped_dupe += 1; continue; end
        push!(seen_faces, key)
        
        # Check C: 3D Geometric Area (Robust against Z-plane issues)
        p1, p2, p3 = vertices[tri[1]], vertices[tri[2]], vertices[tri[3]]
        
        # Vector u = p2 - p1, v = p3 - p1
        ux, uy, uz = p2.x - p1.x, p2.y - p1.y, p2.z - p1.z
        vx, vy, vz = p3.x - p1.x, p3.y - p1.y, p3.z - p1.z
        
        # Cross Product
        cx = uy*vz - uz*vy; cy = uz*vx - ux*vz; cz = ux*vy - uy*vx
        area = 0.5 * sqrt(cx^2 + cy^2 + cz^2)
        
        if area > 1e-12
            push!(valid_faces, tri)
        else
            skipped_area += 1
        end
    end
    
    println("Removed $skipped_area sliver faces (Zero Area)")
    println("Removed $skipped_dupe duplicate faces")
    
    # Compact Vertices
    used_verts = fill(false, length(vertices))
    for f in valid_faces; used_verts[f[1]]=true; used_verts[f[2]]=true; used_verts[f[3]]=true; end
    
    old_to_new = zeros(Int, length(vertices))
    new_verts = Point3D[]
    for (i, is_used) in enumerate(used_verts)
        if is_used
            push!(new_verts, vertices[i])
            old_to_new[i] = length(new_verts)
        end
    end
    
    final_faces = Tuple{Int,Int,Int}[]
    for f in valid_faces
        push!(final_faces, (old_to_new[f[1]], old_to_new[f[2]], old_to_new[f[3]]))
    end
    
    println("Final: $(length(new_verts)) vertices, $(length(final_faces)) faces")
    return new_verts, final_faces
end

# --- Topology Construction ---
function build_topology(vertices, faces)
    n = length(faces)
    dual_adj = [Tuple{Int, Tuple{Int,Int}}[] for _ in 1:n]
    face_ref = Vector{Tuple{Int,Int}}(undef, n)
    vert_to_faces = [Int[] for _ in 1:length(vertices)]
    edge_hist = Dict{Tuple{Int,Int}, Int}()

    for (i, tri) in enumerate(faces)
        face_ref[i] = (tri[1], tri[2])
        for v in [tri[1], tri[2], tri[3]]; push!(vert_to_faces[v], i); end

        for k in 1:3
            v1, v2 = [tri[1], tri[2], tri[3]][k], [tri[1], tri[2], tri[3]][mod1(k+1, 3)]
            key = minmax(v1, v2)
            if haskey(edge_hist, key)
                n_neighbor = edge_hist[key]
                push!(dual_adj[i], (n_neighbor, key)); push!(dual_adj[n_neighbor], (i, key))
            else
                edge_hist[key] = i
            end
        end
    end
    return MeshTopology(vertices, faces, dual_adj, face_ref, vert_to_faces)
end

# --- MIQ Functions (Transport, Constraints, Solve) ---
function get_edge_vector(topo, v1, v2)
    p1, p2 = topo.vertices[v1], topo.vertices[v2]
    return (p2.x - p1.x, p2.y - p1.y, p2.z - p1.z)
end

function compute_transport_angles(topo)
    kappa = Dict{Tuple{Int,Int}, Float64}()
    for i in 1:length(topo.faces)
        for (j, _) in topo.dual_adj[i]
            if i >= j continue end
            vi, vj = topo.faces[i], topo.faces[j]; shared = intersect(vi, vj)
            ev = get_edge_vector(topo, shared[1], shared[2])
            
            ri = topo.face_ref_edges[i]; rvi = get_edge_vector(topo, ri[1], ri[2])
            ai = atan(ev[2], ev[1]) - atan(rvi[2], rvi[1])
            
            rj = topo.face_ref_edges[j]; rvj = get_edge_vector(topo, rj[1], rj[2])
            aj = atan(ev[2], ev[1]) - atan(rvj[2], rvj[1])
            
            kappa[(i, j)] = mod2pi(aj - ai + π) - π
        end
    end
    return kappa
end

function compute_boundary_constraints(topo)
    constraints = Dict{Int, Float64}()
    for i in 1:length(topo.faces)
        tri = topo.faces[i]
        for edge in [(tri[1], tri[2]), (tri[2], tri[3]), (tri[3], tri[1])]
            key = minmax(edge[1], edge[2])
            has_neighbor = false
            for (_, n_key) in topo.dual_adj[i]
                if n_key == key; has_neighbor = true; break; end
            end
            if !has_neighbor
                p1, p2 = topo.vertices[edge[1]], topo.vertices[edge[2]]
                global_tan = atan(p2.y - p1.y, p2.x - p1.x)
                re = topo.face_ref_edges[i]
                rp1, rp2 = topo.vertices[re[1]], topo.vertices[re[2]]
                ref_ang = atan(rp2.y - rp1.y, rp2.x - rp1.x)
                constraints[i] = mod2pi(global_tan - ref_ang + π) - π
            end
        end
    end
    return constraints
end

function initialize_field(topo, constraints)
    kappa = compute_transport_angles(topo)
    n = length(topo.faces); theta = zeros(n)
    p_jumps = Dict{Tuple{Int,Int}, Int}()
    for (f, a) in constraints; theta[f] = a; end
    return CrossField(topo, theta, p_jumps, kappa, constraints, Set{Tuple{Int,Int}}())
end

function build_spanning_tree!(field)
    visited = falses(length(field.topology.faces)); queue = Int[]
    start = isempty(field.constrained_faces) ? 1 : first(keys(field.constrained_faces))
    if start > length(visited); start=1; end
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

function solve_greedy!(field)
    println("--- Solving MIQ ---")
    build_spanning_tree!(field)
    
    topo = field.topology; n_faces = length(topo.faces)
    theta_map = Dict{Int,Int}(); p_map = Dict{Tuple{Int,Int}, Int}(); curr_idx = 0
    
    for f in 1:n_faces
        if !haskey(field.constrained_faces, f); curr_idx += 1; theta_map[f] = curr_idx; end
    end
    for i in 1:n_faces
        for (j, _) in topo.dual_adj[i]
            if i < j && !( (i,j) in field.fixed_edges ) && !(haskey(field.constrained_faces, i) && haskey(field.constrained_faces, j))
                curr_idx += 1; p_map[(i,j)] = curr_idx
            end
        end
    end
    
    n_vars = curr_idx; I = Int[]; J = Int[]; V = Float64[]; b = zeros(n_vars)
    function add!(r, c, val); push!(I, r); push!(J, c); push!(V, val); end
    
    for i in 1:n_faces
        for (j, _) in topo.dual_adj[i]
            if i >= j continue end
            edge = (i, j); kappa = field.transport_angles[edge]
            idx_i = get(theta_map, i, 0); idx_j = get(theta_map, j, 0); idx_p = get(p_map, edge, 0)
            
            val_i, val_j = get(field.constrained_faces, i, 0.0), get(field.constrained_faces, j, 0.0)
            val_p = 0.0
            if idx_p == 0
                if edge in field.fixed_edges; val_p = Float64(field.period_jumps[edge])
                elseif haskey(field.constrained_faces, i) && haskey(field.constrained_faces, j)
                    term = (2.0/π) * (val_j - val_i - kappa)
                    val_p = round(term); field.period_jumps[edge] = Int(val_p)
                end
            end
            
            rhs = -kappa - (idx_i==0 ? val_i : 0.0) + (idx_j==0 ? val_j : 0.0) - (idx_p==0 ? (π/2)*val_p : 0.0)
            
            if idx_i>0; add!(idx_i, idx_i, 2.0); if idx_j>0 add!(idx_i, idx_j, -2.0) end; if idx_p>0 add!(idx_i, idx_p, π) end; b[idx_i]+=2.0*rhs; end
            if idx_j>0; if idx_i>0 add!(idx_j, idx_i, -2.0) end; add!(idx_j, idx_j, 2.0); if idx_p>0 add!(idx_j, idx_p, -π) end; b[idx_j]+=-2.0*rhs; end
            if idx_p>0; if idx_i>0 add!(idx_p, idx_i, π) end; if idx_j>0 add!(idx_p, idx_j, -π) end; add!(idx_p, idx_p, π^2/2.0); b[idx_p]+=π*rhs; end
        end
    end
    
    A = sparse(I, J, V, n_vars, n_vars); x = A \ b
    
    # Greedy Rounding (Simplified for brevity)
    int_indices = collect(values(p_map)); fixed_mask = falses(n_vars); diag_A = [A[k,k] for k in 1:n_vars]
    while true
        best_idx = -1; min_err = Inf; best_val = 0.0
        for idx in int_indices
            if !fixed_mask[idx]
                v = x[idx]; r = round(v); err = abs(v-r)
                if err < min_err; min_err = err; best_idx = idx; best_val = r; end
            end
        end
        if best_idx == -1; break; end
        x[best_idx] = best_val; fixed_mask[best_idx] = true
        
        # Fast local update (Gauss-Seidel)
        rows = rowvals(A); vals = nonzeros(A); q = Queue{Int}(); enqueue!(q, best_idx); in_q = Set{Int}([best_idx])
        iter = 0
        while !isempty(q) && iter < 1000
            k = dequeue!(q); delete!(in_q, k); iter += 1
            Ax_k = 0.0; for nz in nzrange(A, k); Ax_k += vals[nz] * x[rows[nz]]; end
            r = b[k] - Ax_k
            if abs(r) > 1e-5
                x[k] += r / diag_A[k]
                for nz in nzrange(A, k); n = rows[nz]; if !fixed_mask[n] && !(n in in_q); enqueue!(q, n); push!(in_q, n); end; end
            end
        end
    end
    
    for (f, idx) in theta_map; field.theta[f] = x[idx]; end
    for (edge, idx) in p_map; field.period_jumps[edge] = Int(round(x[idx])); end
end

# --- Visualization ---
function plot_results(field, filename)
    fig = Figure(size=(1000, 800)); ax = Axis(fig[1, 1], aspect=DataAspect())
    for tri in field.topology.faces
        pts = [field.topology.vertices[i] for i in tri]; push!(pts, pts[1])
        lines!(ax, [p.x for p in pts], [p.y for p in pts], color=:gray90, linewidth=0.5)
    end
    
    xs, ys, us, vs = Float64[], Float64[], Float64[], Float64[]
    for i in 1:length(field.topology.faces)
        re = field.topology.face_ref_edges[i]
        p1, p2 = field.topology.vertices[re[1]], field.topology.vertices[re[2]]
        ref_ang = atan(p2.y - p1.y, p2.x - p1.x); ang = ref_ang + field.theta[i]
        tri = field.topology.faces[i]
        c = Point3D((field.topology.vertices[tri[1]].x + field.topology.vertices[tri[2]].x + field.topology.vertices[tri[3]].x)/3,
                    (field.topology.vertices[tri[1]].y + field.topology.vertices[tri[2]].y + field.topology.vertices[tri[3]].y)/3, 0.0)
        len = 0.02
        push!(xs, c.x); push!(ys, c.y); push!(us, cos(ang)*len); push!(vs, sin(ang)*len)
        push!(xs, c.x); push!(ys, c.y); push!(us, cos(ang+π/2)*len); push!(vs, sin(ang+π/2)*len)
    end
    arrows!(ax, xs, ys, us, vs, color=:blue, arrowsize=0, lengthscale=1.0)
    save(filename, fig); println("Saved to $filename")
end

# --- MAIN ---
function main()
    filename = "triangulations/mesh_airfoil_dae11.su2"
    println("Reading mesh...")
    raw_verts, raw_faces = read_su2(filename)
    verts, faces = safe_repair_mesh(raw_verts, raw_faces) # USE SAFE REPAIR
    topo = build_topology(verts, faces)
    constraints = compute_boundary_constraints(topo)
    println("Constrained $(length(constraints)) boundary faces.")
    
    # Safety: anchor one face if boundaries are confusing
    if isempty(constraints); constraints = Dict(1 => 0.0); end
    
    field = initialize_field(topo, constraints)
    solve_greedy!(field)
    plot_results(field, "miq_result.png")
end

main()