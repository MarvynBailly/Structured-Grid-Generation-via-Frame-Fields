"""
    visualize_frame_field.jl

Visualizes the frame field solution on a mesh, showing:
- Main frame direction (from θ) as solid lines
- Orthogonal direction as dashed lines
- θ values per face
- κ (transport angles) between faces
- Period jumps p on edges
"""

using Pkg
Pkg.activate(".")

# Add parent directory to load path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using GeometryBasics
using CairoMakie
using LinearAlgebra
using Statistics

include("../src/meshio.jl")
include("../src/mesh_types.jl")

folder_name = "visualize_frame_field"
case_name = "three-triangles"


function main()
    println("Loading three-triangles mesh...")
    # make the directory path if it doesn't exist
    mesh_file = joinpath(@__DIR__, "..", "..", "triangulations", "three-triangles.msh")
    mesh = load_triangulation(mesh_file)
    
    println("Building mesh topology...")
    topology = build_mesh_topology(mesh)
    
    println("Creating initial frame field...")
    # Initialize with zero angles and zero period jumps
    field = CrossField(topology)
    
    # Set some interesting angles for visualization
    println("Setting initial θ values...")
    n_faces = topology.n_faces
    
    println("\nFrame Field Statistics:")
    println("  Faces: ", n_faces)
    println("  Non-zero θ values: ", sum(x -> abs(x) > 1e-10, field.theta))
    println("  Non-zero κ values: ", sum(x -> abs(x) > 1e-10, values(field.transport_angles)))
    println("  Non-zero period jumps: ", sum(x -> x != 0, values(field.period_jumps)))
    
    # Create visualization
    println("\nCreating frame field visualization...")
    fig = Figure(size=(2000, 1000))
    
    # Get mesh data
    vs = coordinates(mesh)
    fs = faces(mesh)
    
    # Convert to 2D for plotting
    points_2d = [Point2f(v[1], v[2]) for v in vs]
    
    # Compute face centroids
    face_centroids = Point2f[]
    for face in fs
        v1, v2, v3 = face[1], face[2], face[3]
        p1, p2, p3 = points_2d[v1], points_2d[v2], points_2d[v3]
        centroid = (p1 + p2 + p3) / 3
        push!(face_centroids, centroid)
    end
    
    # --- Left Panel: Frame Field Directions ---
    ax1 = Axis(fig[1, 1], aspect=DataAspect(), title="Frame Field Directions (θ)")
    
    # Plot triangles
    for face in fs
        v1, v2, v3 = face[1], face[2], face[3]
        p1, p2, p3 = points_2d[v1], points_2d[v2], points_2d[v3]
        lines!(ax1, [p1, p2, p3, p1], color=:gray, linewidth=1)
    end
    
    # Draw reference edges with normals
    for (face_idx, ref_edge) in topology.face_reference_edges
        v1, v2 = ref_edge
        p1, p2 = points_2d[v1], points_2d[v2]
        
        # Draw reference edge as arrow
        edge_vec = p2 - p1
        arrows2d!(ax1, [p1[1]], [p1[2]], [edge_vec[1]], [edge_vec[2]], 
            color=:orange)
        
        # Get the third vertex of the triangle to determine inward normal
        face = fs[face_idx]
        verts = [face[1], face[2], face[3]]
        v3 = setdiff(verts, [v1, v2])[1]
        p3 = points_2d[v3]
        
        # Compute perpendicular direction pointing into triangle
        edge_vec = p2 - p1
        perp_vec = Point2f(-edge_vec[2], edge_vec[1])
        perp_vec = normalize(perp_vec)
        
        # Check if perpendicular points toward p3, if not flip it
        to_p3 = p3 - (p1 + p2) / 2
        if dot(perp_vec, to_p3) < 0
            perp_vec = -perp_vec
        end
        
        # Draw perpendicular arrow pointing into triangle
        perp_arrow_start = (p1 + p2) / 2
        arrows2d!(ax1, [perp_arrow_start[1]], [perp_arrow_start[2]], 
                [perp_vec[1] * 0.05], [perp_vec[2] * 0.05],
                color=:darkorange)
    end
    
    # Draw frame directions at each face centroid
    frame_length = 0.05
    for (i, centroid) in enumerate(face_centroids)
        # Get the 4 frame directions using the new function
        directions = get_frame_directions(field, i)
        
        # Draw main direction (first direction, aligns with reference edge when theta=0)
        dir1 = directions[1]
        arrows2d!(ax1, [centroid[1]], [centroid[2]], 
                  [dir1[1] * frame_length], [dir1[2] * frame_length],
                  color=:blue)
        
        # Draw opposite direction
        dir3 = directions[3]
        arrows2d!(ax1, [centroid[1]], [centroid[2]], 
                  [dir3[1] * frame_length], [dir3[2] * frame_length],
                  color=:red)
        
        # Draw two orthogonal directions
        dir2 = directions[2]
        arrows2d!(ax1, [centroid[1]], [centroid[2]], 
                  [dir2[1] * frame_length], [dir2[2] * frame_length],
                  color=:red)

        dir4 = directions[4]
        arrows2d!(ax1, [centroid[1]], [centroid[2]], 
                  [dir4[1] * frame_length], [dir4[2] * frame_length],
                  color=:red)
        
        # Label with θ value
        theta_deg = round(rad2deg(field.theta[i]), digits=1)
        text!(ax1, centroid, 
              text="θ=$(theta_deg)°", color=:black, fontsize=9)
    end
    
    scatter!(ax1, points_2d, color=:black, markersize=6)
    
    # --- Right Panel: Transport Angles (κ) and Period Jumps (p) ---
    ax2 = Axis(fig[1, 2], aspect=DataAspect(), 
               title="Transport Angles κ and Period Jumps p")
    
    # Plot triangles
    for face in fs
        v1, v2, v3 = face[1], face[2], face[3]
        p1, p2, p3 = points_2d[v1], points_2d[v2], points_2d[v3]
        lines!(ax2, [p1, p2, p3, p1], color=:lightgray, linewidth=1)
    end
    
    # Draw dual edges with κ and p labels
    labeled_edges = Set{Tuple{Int,Int}}()
    for face_i in 1:n_faces
        c_i = face_centroids[face_i]
        
        for (face_j, edge) in topology.dual_adj[face_i]
            if face_i < face_j  # Draw each dual edge once
                c_j = face_centroids[face_j]
                
                # Draw dual edge
                lines!(ax2, [c_i, c_j], color=:purple, linewidth=2, linestyle=:dash)
                
                # Get κ and p values
                kappa_ij = get_transport_angle(field, face_i, face_j)
                p_ij = get_period_jump(field, face_i, face_j)
                
                # Label position (midpoint of dual edge)
                midpoint = (c_i + c_j) / 2
                
                # Format κ value
                kappa_deg = round(rad2deg(kappa_ij), digits=1)
                
                # Offset labels slightly for readability
                text!(ax2, midpoint, 
                      text="κ=$(kappa_deg)°", color=:darkblue, fontsize=9,
                      align=(:center, :center), strokecolor=:black, strokewidth=1)
                
                text!(ax2, midpoint .+ Point2f(0, -0.02),
                      text="p=$p_ij", color=:darkred, fontsize=9,
                      align=(:center, :center), strokecolor=:black, strokewidth=1)
                
                push!(labeled_edges, (face_i, face_j))
            end
        end
    end
    
    # Draw face centroids
    scatter!(ax2, face_centroids, color=:blue, markersize=8)
    for (i, c) in enumerate(face_centroids)
        text!(ax2, c, text="$i", color=:black, fontsize=8, align=(:center, :center))
    end
    
    scatter!(ax2, points_2d, color=:black, markersize=6)
    
    # Save figure
    output_dir = joinpath(@__DIR__, "..", "output", folder_name)
    mkpath(output_dir)
    output_file = joinpath(output_dir, "$(case_name).png")
    save(output_file, fig)
    println("\nVisualization saved to: $output_file")
    
    # Print detailed information
    println("\n" * "="^60)
    println("FRAME FIELD DETAILS")
    println("="^60)
    
    println("\nθ values per face:")
    for i in 1:min(10, n_faces)  # Show first 10
        theta_deg = round(rad2deg(field.theta[i]), digits=2)
        println("  Face $i: θ = $(field.theta[i]) rad = $(theta_deg)°")
    end
    if n_faces > 10
        println("  ... (showing first 10 of $n_faces faces)")
    end
    
    println("\nTransport angles κ (sample):")
    cnt = 0
    for ((i, j), kappa) in field.transport_angles
        if i < j && cnt < 10
            kappa_deg = round(rad2deg(kappa), digits=2)
            println("  Face $i → Face $j: κ = $(round(kappa, digits=4)) rad = $(kappa_deg)°")
            cnt += 1
        end
    end
    
    println("\nPeriod jumps p (non-zero only):")
    non_zero_p = filter(x -> x[2] != 0, field.period_jumps)
    if isempty(non_zero_p)
        println("  (All period jumps are zero)")
    else
        for (edge, p) in non_zero_p
            println("  Edge $edge: p = $p")
        end
    end
    
    println("\n" * "="^60)
    
    return fig
end

# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
