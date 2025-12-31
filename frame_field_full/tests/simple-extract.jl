using FrameFieldFull
using FrameFieldFull.Types
using FrameFieldFull.Topology
using FrameFieldFull.Solver
using FrameFieldFull.Analysis
using FrameFieldFull.MeshIO
using FrameFieldFull.Cutting
using FrameFieldFull.Plotting
using FrameFieldFull.Parameterization
using FrameFieldFull.Extraction
using LinearAlgebra
using Printf

# --- 1. Generate Simple Rectangular Mesh ---
function generate_rectangle_mesh(filename; nx=20, ny=20)
    # Generates a grid from (-1,-1) to (1,1)
    vertices = Tuple{Float64, Float64, Float64}[]
    faces = Tuple{Int, Int, Int}[]
    
    # Create vertices
    for j in 1:(ny+1)
        y = 2.0 * (j-1)/ny - 1.0
        for i in 1:(nx+1)
            x = 2.0 * (i-1)/nx - 1.0
            push!(vertices, (x, y, 0.0))
        end
    end
    
    # Create Triangles
    # Grid indices: (j-1)*(nx+1) + i
    get_idx(i, j) = (j-1)*(nx+1) + i
    
    for j in 1:ny
        for i in 1:nx
            # Quad nodes
            p1 = get_idx(i, j)
            p2 = get_idx(i+1, j)
            p3 = get_idx(i+1, j+1)
            p4 = get_idx(i, j+1)
            
            # Split into 2 triangles
            push!(faces, (p1, p2, p3))
            push!(faces, (p1, p3, p4))
        end
    end
    
    # Write SU2
    open(filename, "w") do io
        write(io, "NDIME= 2\n")
        write(io, "NPOIN= $(length(vertices))\n")
        for v in vertices
            write(io, "$(v[1]) $(v[2]) $(v[3])\n")
        end
        write(io, "NELEM= $(length(faces))\n")
        for f in faces
            write(io, "5 $(f[1]-1) $(f[2]-1) $(f[3]-1)\n")
        end
    end
    println("Generated rectangle mesh: $filename ($nx x $ny cells)")
end

# --- 1. Generate L-Shape Mesh ---
function generate_l_shape_mesh(filename; n_side=20)
    # L-Shape defined on [-1, 1] x [-1, 1]
    # Removing the Top-Right Quadrant (x > 0, y > 0)
    
    vertices = Tuple{Float64, Float64, Float64}[]
    faces = Tuple{Int, Int, Int}[]
    
    # We'll use a coordinate map to handle vertex welding easily
    # key: (row_idx, col_idx) -> vertex_array_index
    # ranges: i in 1..2n+1, j in 1..2n+1
    
    # We map indices (i,j) to coordinates (-1..1)
    # n_side is number of cells per unit length (so 2*n_side total width)
    nx = 2 * n_side
    ny = 2 * n_side
    
    node_map = Dict{Tuple{Int,Int}, Int}()
    
    # Create Vertices
    for j in 1:(ny+1)
        y = 2.0 * (j-1)/ny - 1.0
        for i in 1:(nx+1)
            x = 2.0 * (i-1)/nx - 1.0
            
            # Reentrant Corner Logic:
            # The hole is where x > 0 and y > 0
            # epsilon for floating point safety not needed if we use integer logic
            # Center of the domain is i = n_side + 1, j = n_side + 1 => (0.0, 0.0)
            
            # Vertices existing:
            # 1. x <= 0 (Left half)
            # 2. y <= 0 (Bottom half)
            # 3. The corner (0,0) itself is valid
            
            is_in_hole = (x > 1e-10) && (y > 1e-10)
            
            if !is_in_hole
                push!(vertices, (x, y, 0.0))
                node_map[(i, j)] = length(vertices)
            end
        end
    end
    
    # Create Elements
    for j in 1:ny
        for i in 1:nx
            # Check if this QUAD is valid
            # A quad is valid if all its 4 corners exist in node_map
            # (Strictly, if its centroid is not in the hole)
            
            p1_key = (i, j)
            p2_key = (i+1, j)
            p3_key = (i+1, j+1)
            p4_key = (i, j+1)
            
            if haskey(node_map, p1_key) && haskey(node_map, p2_key) && 
               haskey(node_map, p3_key) && haskey(node_map, p4_key)
                
                idx1 = node_map[p1_key]
                idx2 = node_map[p2_key]
                idx3 = node_map[p3_key]
                idx4 = node_map[p4_key]
                
                # Split into 2 triangles
                push!(faces, (idx1, idx2, idx3))
                push!(faces, (idx1, idx3, idx4))
            end
        end
    end
    
    open(filename, "w") do io
        write(io, "NDIME= 2\n")
        write(io, "NPOIN= $(length(vertices))\n")
        for v in vertices
            write(io, "$(v[1]) $(v[2]) $(v[3])\n")
        end
        write(io, "NELEM= $(length(faces))\n")
        for f in faces
            write(io, "5 $(f[1]-1) $(f[2]-1) $(f[3]-1)\n")
        end
    end
    println("Generated L-Shape mesh: $filename (~$(length(faces)) tris)")
end


# --- 2. Setup Constraints (Force Axis Alignment) ---
function setup_aligned_constraints(topo)
    constraints = Dict{Int, Float64}()
    
    for i in 1:length(topo.faces)
        # Check if boundary face
        re = topo.face_ref_edges[i]
        
        # Calculate Reference Angle in Global Space
        p1 = topo.vertices[re[1]]
        p2 = topo.vertices[re[2]]
        ref_ang = atan(p2.y - p1.y, p2.x - p1.x)
        
        # We want the cross field to be aligned with Global X/Y axes.
        # This corresponds to global angle 0.0 (or pi/2, pi, etc.)
        # theta_local = theta_global - ref_angle
        # We enforce theta_global = 0.0
        
        val = mod2pi(0.0 - ref_ang + π) - π
        constraints[i] = val
    end
    
    # Note: We constrain EVERY face to be axis-aligned. 
    # This is a "hard" debug test. For a real solver, we would only constrain boundaries.
    return constraints
end


# --- 2. Setup Constraints (Boundaries Only) ---
function setup_boundary_constraints(topo)
    constraints = Dict{Int, Float64}()
    
    for i in 1:length(topo.faces)
        # Check if face 'i' is on the boundary
        # In a manifold triangle mesh, internal faces have 3 neighbors.
        # Boundary faces have < 3 neighbors in dual_adj.
        if length(topo.dual_adj[i]) < 3
            
            # We found a boundary face. 
            # We want to force it to align with Global X (Angle 0.0).
            
            re = topo.face_ref_edges[i]
            p1 = topo.vertices[re[1]]
            p2 = topo.vertices[re[2]]
            
            # Reference edge angle in global space
            ref_ang = atan(p2.y - p1.y, p2.x - p1.x)
            
            # theta_local = theta_global - ref_angle
            # Target: theta_global = 0.0
            val = mod2pi(0.0 - ref_ang + π) - π
            
            constraints[i] = val
        end
    end
    
    return constraints
end

# --- 1. Generate Quarter Annulus Mesh ---
function generate_curved_mesh(filename; r_in=1.0, r_out=2.5, n_radial=10, n_angular=30)
    # Generates a grid in Polar Coordinates (r, theta)
    # Theta: 0 to pi/2 (90 degrees)
    
    vertices = Tuple{Float64, Float64, Float64}[]
    faces = Tuple{Int, Int, Int}[]
    
    # Logical grid mapping
    get_idx(i, j) = (j-1)*(n_angular+1) + i
    
    for j in 1:(n_radial+1)
        # Normalized Radial coordinate t_r in [0, 1]
        t_r = (j-1)/n_radial
        r = r_in + t_r * (r_out - r_in)
        
        for i in 1:(n_angular+1)
            # Normalized Angular coordinate t_theta in [0, 1]
            t_theta = (i-1)/n_angular
            theta = t_theta * (π/2)
            
            x = r * cos(theta)
            y = r * sin(theta)
            push!(vertices, (x, y, 0.0))
        end
    end
    
    # Generate Triangles
    for j in 1:n_radial
        for i in 1:n_angular
            p1 = get_idx(i, j)
            p2 = get_idx(i+1, j)
            p3 = get_idx(i+1, j+1)
            p4 = get_idx(i, j+1)
            
            push!(faces, (p1, p2, p3))
            push!(faces, (p1, p3, p4))
        end
    end
    
    open(filename, "w") do io
        write(io, "NDIME= 2\n")
        write(io, "NPOIN= $(length(vertices))\n")
        for v in vertices
            write(io, "$(v[1]) $(v[2]) $(v[3])\n")
        end
        write(io, "NELEM= $(length(faces))\n")
        for f in faces
            write(io, "5 $(f[1]-1) $(f[2]-1) $(f[3]-1)\n")
        end
    end
    println("Generated curved mesh: $filename")
end

function setup_smooth_boundary_constraints(topo)
    constraints = Dict{Int, Float64}()
    
    # 1. Identify Boundary Edges
    # Map: Vertex_Start -> (Vertex_End, Face_Index)
    next_edge_map = Dict{Int, Tuple{Int, Int}}()
    
    for i in 1:length(topo.faces)
        tri = topo.faces[i]
        edges = [(tri[1], tri[2]), (tri[2], tri[3]), (tri[3], tri[1])]
        
        for edge in edges
            # Check if boundary (no neighbor sharing this edge key)
            key = minmax(edge[1], edge[2])
            has_neighbor = false
            for (_, n_key) in topo.dual_adj[i]
                if n_key == key; has_neighbor = true; break; end
            end
            
            if !has_neighbor
                # Found boundary edge!
                # Orient it CCW for the face (already is, by definition of 'edges' list)
                next_edge_map[edge[1]] = (edge[2], i)
            end
        end
    end
    
    # 2. Trace Loops & Propagate
    visited_starts = Set{Int}()
    
    for (start_node, _) in next_edge_map
        if start_node in visited_starts; continue; end
        
        # Start a new loop
        curr = start_node
        
        # State: The "Global" angle we want to maintain
        # Initialize with the first edge's tangent
        # (We could align this to a global guide if requested, but Tangent is a good default start)
        first_next, first_face = next_edge_map[curr]
        p1 = topo.vertices[curr]
        p2 = topo.vertices[first_next]
        current_target_angle = atan(p2.y - p1.y, p2.x - p1.x)
        
        loop_edges = []
        
        while true
            push!(visited_starts, curr)
            if !haskey(next_edge_map, curr); break; end
            
            next_node, face_idx = next_edge_map[curr]
            
            # --- Constraint Logic ---
            # 1. Compute Tangent of current edge
            p_start = topo.vertices[curr]
            p_end   = topo.vertices[next_node]
            tangent_ang = atan(p_end.y - p_start.y, p_end.x - p_start.x)
            
            # 2. We want to pick a constraint 'alpha' such that:
            #    alpha == tangent_ang (mod pi/2)
            #    AND alpha is closest to 'current_target_angle'
            
            # Candidates: Tangent, Normal, -Tangent, -Normal
            # We just iterate k in 0..3
            best_diff = Inf
            best_ang = 0.0
            
            # Check 4 cardinal directions relative to tangent
            for k in 0:3
                candidate = tangent_ang + k * (π/2)
                diff = abs(mod2pi(candidate - current_target_angle + π) - π)
                if diff < best_diff
                    best_diff = diff
                    best_ang = candidate
                end
            end
            
            # Update target for the *next* edge to match this one
            current_target_angle = best_ang
            
            # 3. Convert 'best_ang' (Global) to Local Frame for the solver
            re = topo.face_ref_edges[face_idx]
            rp1, rp2 = topo.vertices[re[1]], topo.vertices[re[2]]
            ref_ang = atan(rp2.y - rp1.y, rp2.x - rp1.x)
            
            val = mod2pi(best_ang - ref_ang + π) - π
            constraints[face_idx] = val
            
            # Move to next
            curr = next_node
            if curr == start_node; break; end # Loop closed
        end
    end
    
    return constraints
end

function setup_analytic_constraints(topo)
    constraints = Dict{Int, Float64}()
    
    for i in 1:length(topo.faces)
        # Check if boundary face
        if length(topo.dual_adj[i]) < 3
            
            # 1. Get Face Centroid
            tri = topo.faces[i]
            p1 = topo.vertices[tri[1]]
            p2 = topo.vertices[tri[2]]
            p3 = topo.vertices[tri[3]]
            
            cx = (p1.x + p2.x + p3.x) / 3.0
            cy = (p1.y + p2.y + p3.y) / 3.0
            
            # 2. Determine Target Global Angle
            
            # Check if we are on the bottom edge (y approx 0, x > 0)
            if (abs(p1.y) < 1e-5 || abs(p2.y) < 1e-5 || abs(p3.y) < 1e-5) && cx > 0
                # Force Negative X Axis Direction
                target_global = π
                constraints[i] = target_global
            else
                # Standard Polar Alignment for other boundaries
                target_global = atan(cy, cx)
                # 3. Convert to Local Frame
                re = topo.face_ref_edges[i]
                rp1, rp2 = topo.vertices[re[1]], topo.vertices[re[2]]
                ref_ang = atan(rp2.y - rp1.y, rp2.x - rp1.x)
                
                # Constraint: Field Angle = Target - Reference
                val = mod2pi(target_global - ref_ang + π) - π
                constraints[i] = val
            end
            
            
        end
    end
    return constraints
end


# --- 2. Setup Tangent Constraints ---
function setup_tangent_constraints(topo)
    constraints = Dict{Int, Float64}()
    
    for i in 1:length(topo.faces)
        # Check if boundary face (fewer than 3 neighbors)
        if length(topo.dual_adj[i]) < 3
            
            # Identify the boundary edge
            tri = topo.faces[i]
            edges = [(tri[1], tri[2]), (tri[2], tri[3]), (tri[3], tri[1])]
            
            boundary_edge = nothing
            for edge in edges
                key = minmax(edge[1], edge[2])
                # Check if this edge has a neighbor
                has_neighbor = false
                for (_, n_key) in topo.dual_adj[i]
                    if n_key == key; has_neighbor = true; break; end
                end
                
                if !has_neighbor
                    boundary_edge = edge
                    break
                end
            end
            
            if boundary_edge !== nothing
                # Calculate Tangent Angle
                p1 = topo.vertices[boundary_edge[1]]
                p2 = topo.vertices[boundary_edge[2]]
                
                dx = p2.x - p1.x
                dy = p2.y - p1.y
                tangent_ang = atan(dy, dx)
                
                # Convert to Local Frame Angle
                re = topo.face_ref_edges[i]
                rp1, rp2 = topo.vertices[re[1]], topo.vertices[re[2]]
                ref_ang = atan(rp2.y - rp1.y, rp2.x - rp1.x)
                
                # Constraint: Theta_local = Tangent_global - Reference_global
                # (plus pi for periodicity, handled by mod2pi)
                val = mod2pi(tangent_ang - ref_ang + π) - π
                constraints[i] = val
            end
        end
    end
    return constraints
end



function generate_curved_mesh_with_constraints(filename; r_in=1.0, r_out=2.5, n_radial=10, n_angular=30)
    # Generates a grid in Polar Coordinates (r, theta)
    # Returns: constraints (Dict{Int, Float64})
    
    vertices = Tuple{Float64, Float64, Float64}[]
    faces = Tuple{Int, Int, Int}[]
    constraints = Dict{Int, Float64}()
    
    # Logical grid mapping
    get_idx(i, j) = (j-1)*(n_angular+1) + i
    
    # 1. Generate Vertices
    for j in 1:(n_radial+1)
        t_r = (j-1)/n_radial
        r = r_in + t_r * (r_out - r_in)
        for i in 1:(n_angular+1)
            t_theta = (i-1)/n_angular
            theta = t_theta * (π/2)
            push!(vertices, (r*cos(theta), r*sin(theta), 0.0))
        end
    end
    
    # 2. Generate Elements & Constraints
    face_idx = 0
    for j in 1:n_radial
        for i in 1:n_angular
            # Grid indices for the current quad
            p1 = get_idx(i, j)      # Bottom-Left
            p2 = get_idx(i+1, j)    # Bottom-Right
            p3 = get_idx(i+1, j+1)  # Top-Right
            p4 = get_idx(i, j+1)    # Top-Left
            
            # --- Triangle 1: (p1, p2, p3) ---
            push!(faces, (p1, p2, p3))
            face_idx += 1
            
            # CONSTRAINT A: Inner Radius (Bottom Arc)
            # This triangle touches the inner radius if j=1
            if j == 1
                t_theta_mid = (i - 0.5) / n_angular
                theta_mid = t_theta_mid * (π/2)
                
                # Normal Inward (Theta + Pi)
                target_global = theta_mid + π
                
                v1 = vertices[p1]
                v2 = vertices[p2]
                ref_ang = atan(v2[2] - v1[2], v2[1] - v1[1])
                
                constraints[face_idx] = mod2pi(target_global - ref_ang + π) - π
            end
            
            # --- Triangle 2: (p1, p3, p4) ---
            push!(faces, (p1, p3, p4))
            face_idx += 1
            
            if j == n_radial
                # CONSTRAINT B: Outer Radius (Top Arc)
                t_theta_mid = (i - 0.5) / n_angular
                theta_mid = t_theta_mid * (π/2)
                
                # Normal Outward (Theta)
                target_global = theta_mid
                
                v1 = vertices[p3]
                v2 = vertices[p4]
                ref_ang = atan(v2[2] - v1[2], v2[1] - v1[1])
                
                constraints[face_idx] = mod2pi(target_global - ref_ang + π) - π
            end

            push!(faces, (p1, p3, p4))
            face_idx += 1
            
            # CONSTRAINT C: Left Edge / Inlet (i=1 implies theta=0)
            if i == 1
                target_global = 0.0 # Pointing Right (+X)
                
                v1 = vertices[p1]
                v2 = vertices[p3]
                
                # Calculate angle of the Reference Edge (p1 -> p3)
                ref_ang = atan(v2[2] - v1[2], v2[1] - v1[1])
                
                # Calculate relative constraint
                constraints[face_idx] = mod2pi(target_global - ref_ang + π) - π
            end
        end
    end
    
    # 3. Write File
    open(filename, "w") do io
        write(io, "NDIME= 2\n")
        write(io, "NPOIN= $(length(vertices))\n")
        for v in vertices
            write(io, "$(v[1]) $(v[2]) $(v[3])\n")
        end
        write(io, "NELEM= $(length(faces))\n")
        for f in faces
            write(io, "5 $(f[1]-1) $(f[2]-1) $(f[3]-1)\n")
        end
    end
    
    println("Generated mesh $filename with $(length(constraints)) constraints.")
    return constraints
end




# --- Main Debug Pipeline ---
function main()
    # mesh_file = "debug_l_shape.su2"
    mesh_file = "debug_curved_constrained.su2"
    # generate_rectangle_mesh(mesh_file; nx=10, ny=10) # Low res for easy debugging
    # generate_l_shape_mesh(mesh_file; n_side=10) # Low res for easy debugging
    # generate_curved_mesh(mesh_file; r_in=1.0, r_out=2.5, n_radial=10, n_angular=10)
    constraints = generate_curved_mesh_with_constraints(mesh_file; 
                                                        r_in=1.0, r_out=2.5, 
                                                        n_radial=10, n_angular=10)
    
    # 1. Load
    verts, faces = read_mesh(mesh_file)
    topo = build_topology(verts, faces)
    
    # 2. Field Solve
    # Constrain ALL faces to ensure perfect alignment
    # constraints = setup_aligned_constraints(topo)
    # constraints = setup_boundary_constraints(topo)
    # constraints = setup_smooth_boundary_constraints(topo)
    # constraints = setup_analytic_constraints(topo)
    field = initialize_field(topo, constraints)
    
    # Just one solve iteration needed since everything is constrained
    solve_greedy!(field; verbose=true)
    optimize_singularities!(field; verbose=true)

    # 3. Parametrization 
    # Since the field is smooth (0 index), cuts should be empty
    sings = compute_singularities(field)
    cuts = compute_cut_graph(topo, sings)

    plot_results(field, "output/debug_frame_field.png", sings; cut_edges=cuts, show_period_jumps=true)
    # println("\n--- Solving Parameterization ---")
    # u_param, v_param = solve_mixed_integer_parameterization(field, cuts)
    
    # # # DEBUG: Print UV range
    # u_min, u_max = minimum(u_param), maximum(u_param)
    # v_min, v_max = minimum(v_param), maximum(v_param)
    # @printf("UV Range: U[%.2f, %.2f], V[%.2f, %.2f]\n", u_min, u_max, v_min, v_max)
    # @printf("UV Dimensions: %.2f x %.2f\n", u_max-u_min, v_max-v_min)
    
    # # 4. Visualization (UV Space)
    # plot_quad_mesh(topo, u_param, v_param, Tuple{Int,Int,Int,Int}[], "output/debug_uv_space.png"; verbose=true)
    
    # # 5. Extraction
    # quad_mesh = extract_quad_mesh(topo, u_param, v_param)
    
    # # 6. Visualization (Final)
    # # plot_extracted_quads(quad_mesh.vertices, quad_mesh.quads, "output/debug_rect_quads.png"; topo=topo, verbose=true)
    # plot_pipeline_overview(field, quad_mesh, "output/pipeline_debug.png"; verbose=true)
end

main()