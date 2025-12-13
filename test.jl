using LinearAlgebra
using SparseArrays
using Printf
using CairoMakie

# ==============================================================================
# 1. DATA STRUCTURES & TOPOLOGY
# ==============================================================================

struct Point3D
    x::Float64
    y::Float64
    z::Float64
end

"""
Stores topological connectivity of the mesh.
"""
struct MeshTopology
    vertices::Vector{Point3D}
    faces::Vector{Tuple{Int,Int,Int}} # 1-based indices
    dual_adj::Vector{Vector{Tuple{Int, Tuple{Int,Int}}}} # face -> [(neighbor, edge_key)]
    face_ref_edges::Vector{Tuple{Int,Int}} # Reference edge for local frames
end

"""
Stores the solution state: theta (angles) and p (integer jumps).
"""
mutable struct CrossField
    topology::MeshTopology
    theta::Vector{Float64}                     # Angle per face
    period_jumps::Dict{Tuple{Int,Int}, Int}    # Integer p_ij per edge
    transport_angles::Dict{Tuple{Int,Int}, Float64} # kappa_ij
    constrained_faces::Dict{Int, Float64}      # face_idx => angle
    fixed_edges::Set{Tuple{Int,Int}}           # Edges where p is fixed
end

# ==============================================================================
# 2. FILE I/O (GMSH 4.1 Parser)
# ==============================================================================

"""
    read_msh(filename)

Parses a Gmsh 4.1 ASCII file and returns vertices and triangle faces.
"""
function read_msh(filename::String)
    lines = readlines(filename)
    vertices = Point3D[]
    faces = Tuple{Int,Int,Int}[]
    
    # Map NodeTag (from file) -> Array Index (1-based)
    tag_map = Dict{Int, Int}()
    
    i = 1
    while i <= length(lines)
        line = strip(lines[i])
        
        if line == "\$Nodes"
            i += 1
            # Header: numEntityBlocks numNodes minNodeTag maxNodeTag
            dims = parse.(Int, split(lines[i]))
            num_blocks = dims[1]
            i += 1
            
            for _ in 1:num_blocks
                # Block Header: entityDim entityTag parametric numNodesInBlock
                block_head = parse.(Int, split(lines[i]))
                n_nodes_block = block_head[4]
                i += 1
                
                # In Gmsh 4.1, tags are listed first, then coordinates
                tags = Int[]
                for _ in 1:n_nodes_block
                    push!(tags, parse(Int, lines[i]))
                    i += 1
                end
                
                for k in 1:n_nodes_block
                    coords = parse.(Float64, split(lines[i]))
                    # Gmsh provides x y z
                    push!(vertices, Point3D(coords[1], coords[2], coords[3]))
                    tag_map[tags[k]] = length(vertices)
                    i += 1
                end
            end
            
        elseif line == "\$Elements"
            i += 1
            # Header: numEntityBlocks numElements minTag maxTag
            dims = parse.(Int, split(lines[i]))
            num_blocks = dims[1]
            i += 1
            
            for _ in 1:num_blocks
                # Block Header: entityDim entityTag type numElementsInBlock
                block_head = parse.(Int, split(lines[i]))
                elem_type = block_head[3]
                n_elems_block = block_head[4]
                i += 1
                
                if elem_type == 2 # Type 2 is 3-node Triangle
                    for _ in 1:n_elems_block
                        # Format: elementTag nodeTag1 nodeTag2 nodeTag3
                        vals = parse.(Int, split(lines[i]))
                        # Convert file tags to local 1-based indices
                        n1 = tag_map[vals[2]]
                        n2 = tag_map[vals[3]]
                        n3 = tag_map[vals[4]]
                        push!(faces, (n1, n2, n3))
                        i += 1
                    end
                else
                    # Skip non-triangle elements (points, lines, tets, etc.)
                    i += n_elems_block
                end
            end
        else
            i += 1
        end
    end
    
    return vertices, faces
end

# ==============================================================================
# 3. GEOMETRY KERNEL
# ==============================================================================

function build_topology(vertices::Vector{Point3D}, faces::Vector{Tuple{Int,Int,Int}})
    n_faces = length(faces)
    dual_adj = [Tuple{Int, Tuple{Int,Int}}[] for _ in 1:n_faces]
    face_ref_edges = Vector{Tuple{Int,Int}}(undef, n_faces)
    
    # Map (v_min, v_max) -> face_index to find neighbors
    edge_history = Dict{Tuple{Int,Int}, Int}()
    
    for (i, tri) in enumerate(faces)
        # Define reference edge as the first edge of the triangle (v1->v2)
        face_ref_edges[i] = (tri[1], tri[2])
        
        # Process edges
        vs = [tri[1], tri[2], tri[3]]
        for k in 1:3
            v1, v2 = vs[k], vs[mod1(k+1, 3)]
            edge_key = minmax(v1, v2)
            
            if haskey(edge_history, edge_key)
                neighbor = edge_history[edge_key]
                push!(dual_adj[i], (neighbor, edge_key))
                push!(dual_adj[neighbor], (i, edge_key))
            else
                edge_history[edge_key] = i
            end
        end
    end
    
    return MeshTopology(vertices, faces, dual_adj, face_ref_edges)
end

function get_edge_vector(topo::MeshTopology, v1::Int, v2::Int)
    p1 = topo.vertices[v1]
    p2 = topo.vertices[v2]
    return (p2.x - p1.x, p2.y - p1.y, p2.z - p1.z)
end

function compute_transport_angles(topo::MeshTopology)
    kappa = Dict{Tuple{Int,Int}, Float64}()
    
    for face_i in 1:length(topo.faces)
        for (face_j, _) in topo.dual_adj[face_i]
            if face_i >= face_j continue end # Canonical only
            
            # 1. Identify Shared Edge
            verts_i = topo.faces[face_i]
            verts_j = topo.faces[face_j]
            shared = intersect(verts_i, verts_j)
            u, v = shared[1], shared[2] # Vertices of shared edge
            
            # Shared edge vector
            edge_vec = get_edge_vector(topo, u, v)
            
            # 2. Angle in Face I (relative to ref edge I)
            ref_i = topo.face_ref_edges[face_i]
            ref_vec_i = get_edge_vector(topo, ref_i[1], ref_i[2])
            # atan(y, x) -> result in (-pi, pi]
            angle_i = atan(edge_vec[2], edge_vec[1]) - atan(ref_vec_i[2], ref_vec_i[1])
            
            # 3. Angle in Face J (relative to ref edge J)
            ref_j = topo.face_ref_edges[face_j]
            ref_vec_j = get_edge_vector(topo, ref_j[1], ref_j[2])
            angle_j = atan(edge_vec[2], edge_vec[1]) - atan(ref_vec_j[2], ref_vec_j[1])
            
            # 4. Kappa = Angle_J - Angle_I
            k_val = angle_j - angle_i
            k_val = mod2pi(k_val + π) - π
            
            kappa[(face_i, face_j)] = k_val
        end
    end
    return kappa
end

# ==============================================================================
# 4. MIQ SOLVER (Greedy Mixed-Integer)
# ==============================================================================

function initialize_field(topo::MeshTopology, constraints::Dict{Int, Float64})
    kappa = compute_transport_angles(topo)
    n_faces = length(topo.faces)
    theta = zeros(n_faces)
    p_jumps = Dict{Tuple{Int,Int}, Int}()
    fixed = Set{Tuple{Int,Int}}()
    
    # Apply initial constraints
    for (f, ang) in constraints
        theta[f] = ang
    end
    
    return CrossField(topo, theta, p_jumps, kappa, constraints, fixed)
end

function build_spanning_tree!(field::CrossField)
    visited = falses(length(field.topology.faces))
    queue = Int[]
    
    # Start from a constrained face (or face 1)
    start_node = isempty(field.constrained_faces) ? 1 : first(keys(field.constrained_faces))
    push!(queue, start_node)
    visited[start_node] = true
    
    while !isempty(queue)
        u = popfirst!(queue)
        for (v, _) in field.topology.dual_adj[u]
            if !visited[v]
                visited[v] = true
                push!(queue, v)
                # Fix p=0 on tree edges
                key = minmax(u, v)
                push!(field.fixed_edges, key)
                field.period_jumps[key] = 0
            end
        end
    end
end

function solve_relaxation!(field::CrossField)
    topo = field.topology
    n_faces = length(topo.faces)
    
    # Identify free edges (variables)
    free_edges = Tuple{Int,Int}[]
    for i in 1:n_faces
        for (j, _) in topo.dual_adj[i]
            if i < j && !((i, j) in field.fixed_edges)
                push!(free_edges, (i, j))
            end
        end
    end
    
    n_free_p = length(free_edges)
    n_vars = n_faces + n_free_p
    
    # Map edge -> variable index offset
    edge_to_var = Dict(e => k + n_faces for (k, e) in enumerate(free_edges))
    
    I, J, V = Int[], Int[], Float64[]
    b = zeros(n_vars)
    
    # Build System
    for i in 1:n_faces
        for (j, _) in topo.dual_adj[i]
            if i >= j continue end
            
            k_ij = field.transport_angles[(i, j)]
            edge = (i, j)
            
            p_val = 0.0
            p_var_idx = 0
            
            if edge in field.fixed_edges
                p_val = Float64(field.period_jumps[edge])
            else
                p_var_idx = edge_to_var[edge]
            end
            
            # Energy term: (theta_i - theta_j + p*pi/2 + kappa)^2
            indices = [i, j]
            coeffs = [1.0, -1.0]
            if p_var_idx > 0
                push!(indices, p_var_idx); push!(coeffs, π/2)
            end
            
            rhs = -(k_ij + p_val * π/2)
            
            # Normal Equations Accumulation
            for r in 1:length(indices)
                row = indices[r]
                b[row] += coeffs[r] * rhs
                for c in 1:length(indices)
                    push!(I, row); push!(J, indices[c]); push!(V, coeffs[r] * coeffs[c])
                end
            end
        end
    end
    
    # Add Constraints
    penalty = 10000.0
    for (f, ang) in field.constrained_faces
        push!(I, f); push!(J, f); push!(V, penalty)
        b[f] += penalty * ang
    end
    
    # Solve
    A = sparse(I, J, V, n_vars, n_vars)
    x = A \ b
    
    field.theta = x[1:n_faces]
    return free_edges, x[n_faces+1:end]
end

function solve_miq!(field::CrossField)
    println("--- Starting MIQ Solver ---")
    build_spanning_tree!(field)
    println("Spanning tree fixed $(length(field.fixed_edges)) edges.")
    
    iter = 0
    while true
        iter += 1
        free_edges, p_vals = solve_relaxation!(field)
        
        if isempty(free_edges) break end
        
        # Greedy Selection: Find p closest to integer
        best_err = Inf
        best_idx = -1
        
        for k in 1:length(p_vals)
            err = abs(p_vals[k] - round(p_vals[k]))
            if err < best_err
                best_err = err
                best_idx = k
            end
        end
        
        # Round and fix
        target = free_edges[best_idx]
        val = Int(round(p_vals[best_idx]))
        
        field.period_jumps[target] = val
        push!(field.fixed_edges, target)
        
        println("Iter $iter: Fixed edge $target to $val (Err: $(round(best_err, digits=4)))")
        println("  Remaining free edges: $(length(free_edges) - 1)")
    end
    
    # Final smooth
    solve_relaxation!(field)
    println("--- Solver Complete ---")
end

# ==============================================================================
# 5. EXECUTION
# ==============================================================================

function run_msh_example()
    filename = "triangulations/three-triangles.msh"
    
    if !isfile(filename)
        println("Error: File '$filename' not found.")
        println("Please save your .msh content to '$filename' before running.")
        return
    end
    
    println("Reading $filename...")
    verts, faces = read_msh(filename)
    println("Loaded $(length(verts)) vertices and $(length(faces)) triangles.")
    
    topo = build_topology(verts, faces)
    
    # Constrain the first face to 0 radians (alignment with X axis)
    # In a real app, you would select specific faces or use user input
    constraints = Dict(1 => 0.0)
    field = initialize_field(topo, constraints)
    
    solve_miq!(field)
    
    println("\nTop 5 Face Angles:")
    for i in 1:min(5, length(faces))
        deg = rad2deg(field.theta[i])
        println("Face $i: $(round(deg, digits=2))°")
    end
    
    println("\nDetected Singularities (Non-zero Period Jumps):")
    count = 0
    for (edge, p) in field.period_jumps
        if p != 0
            count += 1
            if count <= 10
                println("Edge $edge: jump $p")
            end
        end
    end
    println("(Total $count singular edges)")
end

# Run
# run_msh_example()


# ==============================================================================
# 6. VISUALIZATION UTILS (UPDATED)
# ==============================================================================

"""
Helper: Compute the centroid of a face for plotting markers.
"""
function get_face_centroid(topo::MeshTopology, face_idx::Int)
    tri = topo.faces[face_idx]
    p1, p2, p3 = topo.vertices[tri[1]], topo.vertices[tri[2]], topo.vertices[tri[3]]
    return Point2f((p1.x+p2.x+p3.x)/3, (p1.y+p2.y+p3.y)/3)
end

"""
Helper: Get the global rotation angle for the cross at a specific face.
"""
function get_global_angle(field::CrossField, face_idx::Int)
    ref_edge = field.topology.face_ref_edges[face_idx]
    p1 = field.topology.vertices[ref_edge[1]]
    p2 = field.topology.vertices[ref_edge[2]]
    ref_angle = atan(p2.y - p1.y, p2.x - p1.x)
    return ref_angle + field.theta[face_idx]
end

"""
    plot_mesh!(ax, topo)
Draws the wireframe of the mesh.
"""
function plot_mesh!(ax, topo::MeshTopology)
    for tri in topo.faces
        pts = [topo.vertices[i] for i in tri]
        push!(pts, pts[1]) 
        lines!(ax, [p.x for p in pts], [p.y for p in pts], color=:gray80, linewidth=1)
    end
end

"""
    plot_field!(ax, field; color=:black, scale=0.1)
Draws the cross field vectors.
"""
function plot_field!(ax, field::CrossField; color=:blue, scale=0.03)
    xs, ys, us, vs = Float64[], Float64[], Float64[], Float64[]
    colors = []

    for i in 1:length(field.topology.faces)
        # Skip constrained faces in the main drawing (we will draw them special)
        if haskey(field.constrained_faces, i)
            continue
        end
        
        c = get_face_centroid(field.topology, i)
        angle = get_global_angle(field, i)
        
        # 4 directions
        for k in 0:3
            a = angle + k * π/2
            push!(xs, c[1]); push!(ys, c[2])
            push!(us, cos(a)*scale); push!(vs, sin(a)*scale)
            if k == 0
                push!(colors, :blue) # Blue for main direction
            else
                push!(colors, :green) # Green for others
            end
        end
    end
    
    arrows2d!(ax, xs, ys, us, vs, color=colors)
end

"""
    plot_constraints!(ax, field)
Highlights the constrained faces and draws their fixed frame.
"""
function plot_constraints!(ax, field::CrossField; scale=0.2)
    topo = field.topology
    
    colors = []

    # 1. Fill the constrained faces
    for (f_idx, _) in field.constrained_faces
        tri = topo.faces[f_idx]
        pts = [Point2f(topo.vertices[i].x, topo.vertices[i].y) for i in tri]
        # Draw filled polygon
        poly!(ax, pts, color=(:red, 0.2), strokewidth=2, strokecolor=:red)
    end
    
    # 2. Draw the constrained vectors (Thicker, Red)
    xs, ys, us, vs = Float64[], Float64[], Float64[], Float64[]
    for (f_idx, _) in field.constrained_faces
        c = get_face_centroid(topo, f_idx)
        angle = get_global_angle(field, f_idx)
        
        for k in 0:3
            a = angle + k * π/2
            push!(xs, c[1]); push!(ys, c[2])
            push!(us, cos(a)*scale); push!(vs, sin(a)*scale)
            if k == 0
                push!(colors, :red) # Red for main direction
            else
                push!(colors, :orange) # Orange for others
            end
        end
    end
    
    arrows2d!(ax, xs, ys, us, vs, color=colors)
end

# ==============================================================================
# 7. EXPORT FUNCTIONS (UPDATED)
# ==============================================================================

function save_comparison_png(initial_field::CrossField, smooth_field::CrossField, filename::String)
    f = Figure(size = (1200, 600))
    
    ax1 = Axis(f[1, 1], title = "Initial (Unoptimized)", aspect = DataAspect())
    plot_mesh!(ax1, initial_field.topology)
    plot_constraints!(ax1, initial_field)  # Highlight constraints
    plot_field!(ax1, initial_field, color=(:blue, 0.5)) # Dim the rest
    hidedecorations!(ax1)
    
    ax2 = Axis(f[1, 2], title = "MIQ Solution", aspect = DataAspect())
    plot_mesh!(ax2, smooth_field.topology)
    plot_constraints!(ax2, smooth_field)   # Highlight constraints
    plot_field!(ax2, smooth_field, color=:blue)
    hidedecorations!(ax2)
    
    save(filename, f)
    println("Saved comparison to $filename")
end

function solve_miq_animated!(field::CrossField, filename::String)
    println("--- Starting Animated Solver ---")
    
    fig = Figure(size = (800, 800))
    ax = Axis(fig[1, 1], title = "MIQ Optimization", aspect = DataAspect())
    hidedecorations!(ax)
    
    build_spanning_tree!(field)
    
    record(fig, filename; framerate = 5) do io
        # Function to redraw frame
        function draw_frame(title, color)
            empty!(ax)
            plot_mesh!(ax, field.topology)
            plot_constraints!(ax, field) # Always show constraints on top
            plot_field!(ax, field, color=color)
            ax.title = title
            recordframe!(io)
        end

        draw_frame("Step 0: Initialization", :red)
        
        iter = 0
        while true
            iter += 1
            free_edges, p_vals = solve_relaxation!(field)
            
            draw_frame("Step $iter: Continuous Relaxation", :orange)
            
            if isempty(free_edges) break end
            
            # Greedy Rounding
            best_err = Inf; best_idx = -1
            for k in 1:length(p_vals)
                err = abs(p_vals[k] - round(p_vals[k]))
                if err < best_err; best_err = err; best_idx = k; end
            end
            
            target = free_edges[best_idx]
            val = Int(round(p_vals[best_idx]))
            
            field.period_jumps[target] = val
            push!(field.fixed_edges, target)
            
            draw_frame("Step $iter: Fixed Edge to $val", :blue)
            println("Iter $iter: Fixed edge $target to $val (Err: $(round(best_err, digits=4)))")
            println("  Remaining free edges: $(length(free_edges) - 1)")
        end
        
        solve_relaxation!(field)
        for _ in 1:5
            draw_frame("Final Solution", :green)
        end
    end
    
    println("Saved animation to $filename")
end

function run_visual_example()
    # 1. Load Data
    # filename = "triangulations/three-triangles.msh"
    # filename = "triangulations/simple-square.msh"
    filename = "triangulations/disk-radial.msh"
    # filename = "triangulations/hole.msh"
    verts, faces = read_msh(filename)
    topo = build_topology(verts, faces)
    constraints = Dict(1 => 2.3 - π/6)#, 50 => -2.3 + π + π/6 - π/3)
    
    # 2. Create Initial State (Copy for comparison)
    field_initial = initialize_field(topo, constraints)
    # Run just the relaxation once on the initial field to propagate constraints roughly
    # (Optional, but makes the "Before" look less empty)
    solve_relaxation!(field_initial) 
    
    # 3. Create Working Field for Solver
    field_smooth = initialize_field(topo, constraints)
    
    # 4. Run Animated Solver (Generates GIF)
    solve_miq_animated!(field_smooth, "solver_process.gif")
    
    # 5. Save Comparison (Generates PNG)
    save_comparison_png(field_initial, field_smooth, "comparison.png")
end

# run_msh_example()
run_visual_example()