using FrameFieldFull
using FrameFieldFull.Types
using FrameFieldFull.Topology
using FrameFieldFull.Solver
using FrameFieldFull.Analysis
using FrameFieldFull.MeshIO
using FrameFieldFull.Cutting
using LinearAlgebra

include("visualize_cross_field.jl")





# --- 1. Procedural Mesh Generator (NACA 0012 in a Box) ---
function naca0012_y(x)
    t = 0.12
    return 5*t * (0.2969*sqrt(x) - 0.1260*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^4)
end

function generate_airfoil_mesh(filename; n_angular=60, n_radial=5)
    # n_angular: Points around the airfoil (circumferential)
    # n_radial: Layers from airfoil to farfield (normal direction)
    
    vertices = Tuple{Float64, Float64, Float64}[]
    faces = Tuple{Int, Int, Int}[]
    
    # We will store vertex indices in a 2D grid: grid[layer, angular_idx]
    # This makes triangulation much easier.
    grid_indices = zeros(Int, n_radial, n_angular)
    
    # --- 1. Generate Vertices Layer by Layer ---
    for i in 1:n_radial
        # t goes from 0.0 (Airfoil Surface) to 1.0 (Outer Boundary)
        # Using a power law (t^1.2) clusters points slightly closer to the airfoil
        raw_t = (i - 1) / (n_radial - 1)
        t = raw_t^1.2 
        
        for j in 1:n_angular
            # angular parameter 0 to 2pi (CCW from Trailing Edge)
            # j=1 -> theta=0, j=n_angular -> theta ~= 2pi
            theta = ((j - 1) / n_angular) * 2π
            
            # -- A. Inner Point (Airfoil) --
            x_in = 0.5 * (1.0 + cos(theta))
            yt = naca0012_y(x_in)
            y_in = (theta <= π) ? yt : -yt
            
            # Fix LE/TE precision
            if abs(theta - π) < 1e-6; y_in = 0.0; x_in = 0.0; end
            if abs(theta) < 1e-6 || abs(theta - 2π) < 1e-6; y_in = 0.0; x_in = 1.0; end
            
            # -- B. Outer Point (Circle) --
            R_outer = 3.0
            center_x = 0.5
            x_out = center_x + R_outer * cos(theta)
            y_out = R_outer * sin(theta)
            
            # -- C. Interpolate --
            x = (1 - t) * x_in + t * x_out
            y = (1 - t) * y_in + t * y_out
            
            push!(vertices, (x, y, 0.0))
            grid_indices[i, j] = length(vertices)
        end
    end
    
    # --- 2. Triangulate Layers ---
    # Connect Layer i to Layer i+1
    for i in 1:(n_radial - 1)
        for j in 1:n_angular
            j_next = mod1(j + 1, n_angular)
            
            # Indices in the grid
            v_curr_inner = grid_indices[i, j]
            v_next_inner = grid_indices[i, j_next]
            v_curr_outer = grid_indices[i+1, j]
            v_next_outer = grid_indices[i+1, j_next]
            
            # Quad: (v_curr_inner, v_curr_outer, v_next_outer, v_next_inner)
            # Split into 2 Triangles (CCW)
            
            # Tri 1: (Inner1, Outer1, Inner2)
            push!(faces, (v_curr_inner, v_curr_outer, v_next_inner))
            
            # Tri 2: (Outer1, Outer2, Inner2)
            push!(faces, (v_curr_outer, v_next_outer, v_next_inner))
        end
    end
    
    # --- 3. Write SU2 ---
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
end
# --- 2. Custom Constraints ---

function setup_airfoil_constraints(topo; constrain_outer=true)
    constraints = Dict{Int, Float64}()
    
    for i in 1:length(topo.faces)
        tri = topo.faces[i]
        edges = [(tri[1], tri[2]), (tri[2], tri[3]), (tri[3], tri[1])]
        
        for (k, edge) in enumerate(edges)
            key = minmax(edge[1], edge[2])
            
            # Check if boundary
            has_neighbor = false
            for (_, n_key) in topo.dual_adj[i]
                if n_key == key; has_neighbor = true; break; end
            end
            
            if !has_neighbor
                # It is a boundary. Identify which one.
                p1 = topo.vertices[edge[1]]
                p2 = topo.vertices[edge[2]]
                mid_x = (p1.x + p2.x)/2
                mid_y = (p1.y + p2.y)/2
                dist_from_airfoil = sqrt((mid_x - 0.5)^2 + mid_y^2)
                
                # Setup Local Frame
                re = topo.face_ref_edges[i]
                rp1, rp2 = topo.vertices[re[1]], topo.vertices[re[2]]
                ref_ang = atan(rp2.y - rp1.y, rp2.x - rp1.x)
                
                dx = p2.x - p1.x
                dy = p2.y - p1.y
                tangent_ang = atan(dy, dx)
                
                if dist_from_airfoil < 1.0 # Threshold to distinguish inner/outer
                    # === Inner Boundary (Airfoil) ===
                    # ALWAYS Constrain to Tangent
                    val = mod2pi(tangent_ang - ref_ang + π) - π
                    constraints[i] = val
                    
                else
                    # === Outer Boundary (Farfield) ===
                    if constrain_outer
                        # Force Horizontal (Free Stream)
                        val = mod2pi(0.0 - ref_ang + π) - π
                        constraints[i] = val
                    else
                        # Leave Unconstrained (Natural BC)
                        # Do nothing -> Solver will find optimal angle
                    end
                end
            end
        end
    end
    return constraints
end


# --- 3. Run Pipeline ---
function main()
    upper_constraint = true
    mesh_file = "airfoil_test.su2"
    save_path = "output/airfoil_upper$(upper_constraint)_cross_field.png"
    @info("Generating mesh...")
    generate_airfoil_mesh(mesh_file; n_angular=50, n_radial=10)

    
    @info("Reading mesh...")
    verts, faces = read_mesh(mesh_file)
    topo = build_topology(verts, faces)
    
    @info("Setting constraints...")
    constraints = setup_airfoil_constraints(topo, constrain_outer=upper_constraint)
    
    @info("Initializing field...")
    field = initialize_field(topo, constraints)
    
    @info("Solving...")
    solve_greedy!(field; verbose=false)
    
    @info("Optimizing Singularities...")
    optimize_singularities!(field; verbose=false)
    
    # Check Topology
    sings = compute_singularities(field)
    if isempty(sings)
        @info("Total Singularity Index: 0.0 (Perfect!)")
    else
        total_idx = sum(s[2] for s in sings)
        @info("Total Singularity Index: $total_idx (Expected: 0.0)")
    end

    # cut the mesh M into a disk topology
    cuts = compute_cut_graph(topo, sings)

    @info("Visualizing...")
    visualize_field(field, title="Airfoil Flow (Chi=0)", final_cuts=cuts, save_path=save_path, two_d=true)
end

main()