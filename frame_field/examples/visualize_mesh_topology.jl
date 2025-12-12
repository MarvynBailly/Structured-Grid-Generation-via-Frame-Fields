"""
    visualize_mesh_topology.jl

Visualizes mesh topology and cross field structures for the simplest-square mesh.
Creates a PNG showing:
- Triangle faces with indices
- Edges with indices
- Dual adjacency graph
- Reference edges per face
"""

using Pkg
Pkg.activate(".")

# Add parent directory to load path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using GeometryBasics
using CairoMakie
using LinearAlgebra

include("../src/meshio.jl")
include("../src/mesh_types.jl")


folder_name = "visualize_mesh_topology"
case_name = "three-triangles"



function main()
    println("Loading $(case_name) mesh...")
    mesh_file = joinpath(@__DIR__, "..", "..", "triangulations", "$(case_name).msh")
    mesh = load_triangulation(mesh_file)
    
    println("Building mesh topology...")
    topology = build_mesh_topology(mesh)
    
    println("\nMesh Statistics:")
    println("  Vertices: ", length(coordinates(mesh)))
    println("  Faces: ", topology.n_faces)
    println("  Edges: ", topology.n_edges)
    
    # Create visualization
    println("\nCreating visualization...")
    fig = Figure(resolution=(1600, 1200))
    
    # Get mesh data
    vs = coordinates(mesh)
    fs = faces(mesh)
    
    # Convert to 2D for plotting
    points_2d = [Point2f(v[1], v[2]) for v in vs]
    
    # --- Panel 1: Mesh with face indices ---
    ax1 = Axis(fig[1, 1], aspect=DataAspect(), title="Triangle Faces")
    
    # Plot triangles
    for (i, face) in enumerate(fs)
        v1, v2, v3 = face[1], face[2], face[3]
        p1, p2, p3 = points_2d[v1], points_2d[v2], points_2d[v3]
        
        # Draw triangle edges
        lines!(ax1, [p1, p2, p3, p1], color=:black, linewidth=2)
        
        # Compute centroid and label face
        centroid = (p1 + p2 + p3) / 3
        text!(ax1, centroid, text="F$i", color=:blue, fontsize=16, align=(:center, :center))
    end
    
    # Plot vertices
    scatter!(ax1, points_2d, color=:red, markersize=12)
    for (i, p) in enumerate(points_2d)
        text!(ax1, p .+ Point2f(0.05, 0.05), text="v$i", color=:red, fontsize=12)
    end
    
    # --- Panel 2: Edges with indices ---
    ax2 = Axis(fig[1, 2], aspect=DataAspect(), title="Edges")
    
    # Plot triangles lightly
    for face in fs
        v1, v2, v3 = face[1], face[2], face[3]
        p1, p2, p3 = points_2d[v1], points_2d[v2], points_2d[v3]
        lines!(ax2, [p1, p2, p3, p1], color=:gray, linewidth=1)
    end
    
    # Plot and label edges
    for (edge, edge_idx) in topology.edge_map
        v1, v2 = edge
        p1, p2 = points_2d[v1], points_2d[v2]
        midpoint = (p1 + p2) / 2
        
        # Draw edge
        lines!(ax2, [p1, p2], color=:green, linewidth=3)
        
        # Label edge
        text!(ax2, midpoint, text="e$edge_idx", color=:darkgreen, fontsize=14, 
              align=(:center, :center), strokecolor=:white, strokewidth=2)
    end
    
    scatter!(ax2, points_2d, color=:red, markersize=8)
    
    # --- Panel 3: Dual adjacency graph ---
    ax3 = Axis(fig[2, 1], aspect=DataAspect(), title="Dual Adjacency Graph")
    
    # Plot triangles lightly
    for face in fs
        v1, v2, v3 = face[1], face[2], face[3]
        p1, p2, p3 = points_2d[v1], points_2d[v2], points_2d[v3]
        lines!(ax3, [p1, p2, p3, p1], color=:lightgray, linewidth=1)
    end
    
    # Compute face centroids
    face_centroids = Point2f[]
    for face in fs
        v1, v2, v3 = face[1], face[2], face[3]
        p1, p2, p3 = points_2d[v1], points_2d[v2], points_2d[v3]
        centroid = (p1 + p2 + p3) / 3
        push!(face_centroids, centroid)
    end
    
    # Draw dual edges (between face centroids)
    for (face_i, neighbors) in topology.dual_adj
        c_i = face_centroids[face_i]
        for (face_j, edge) in neighbors
            if face_i < face_j  # Draw each dual edge once
                c_j = face_centroids[face_j]
                lines!(ax3, [c_i, c_j], color=:purple, linewidth=2, linestyle=:dash)
                
                # Label with edge info
                midpoint = (c_i + c_j) / 2
                edge_idx = topology.edge_map[edge]
                text!(ax3, midpoint, text="e$edge_idx", color=:purple, fontsize=10,
                      align=(:center, :center), strokecolor=:white, strokewidth=1)
            end
        end
    end
    
    # Draw face centroids
    scatter!(ax3, face_centroids, color=:blue, markersize=15)
    for (i, c) in enumerate(face_centroids)
        text!(ax3, c, text="F$i", color=:black, fontsize=12, align=(:center, :center))
    end
    
    # --- Panel 4: Reference edges ---
    ax4 = Axis(fig[2, 2], aspect=DataAspect(), title="Reference Edges per Face")
    
    # Plot triangles
    for face in fs
        v1, v2, v3 = face[1], face[2], face[3]
        p1, p2, p3 = points_2d[v1], points_2d[v2], points_2d[v3]
        lines!(ax4, [p1, p2, p3, p1], color=:gray, linewidth=1)
    end
    
    # Highlight reference edges with directional arrows into triangles
    for (face_idx, ref_edge) in topology.face_reference_edges
        v1, v2 = ref_edge
        p1, p2 = points_2d[v1], points_2d[v2]
        
        # Draw reference edge thicker
        lines!(ax4, [p1, p2], color=:orange, linewidth=4)
        
        # Get the third vertex of the triangle to determine inward normal
        face = fs[face_idx]
        verts = [face[1], face[2], face[3]]
        v3 = setdiff(verts, [v1, v2])[1]
        p3 = points_2d[v3]
        
        # Compute perpendicular direction pointing into triangle
        edge_vec = p2 - p1
        # Perpendicular vector (rotate 90 degrees)
        perp_vec = Point2f(-edge_vec[2], edge_vec[1])
        perp_vec = normalize(perp_vec)
        
        # Check if perpendicular points toward p3, if not flip it
        to_p3 = p3 - (p1 + p2) / 2
        if dot(perp_vec, to_p3) < 0
            perp_vec = -perp_vec
        end
        
        # Draw arrow along edge showing direction
        edge_direction = normalize(p2 - p1)
        arrow_start = (p1 + p2) / 2
        arrows!(ax4, [arrow_start[1]], [arrow_start[2]], 
                [edge_direction[1] * 0.15], [edge_direction[2] * 0.15],
                color=:orange, arrowsize=15)
        
        # Draw smaller perpendicular arrow pointing into triangle
        perp_arrow_start = (p1 + p2) / 2
        arrows!(ax4, [perp_arrow_start[1]], [perp_arrow_start[2]], 
                [perp_vec[1] * 0.1], [perp_vec[2] * 0.1],
                color=:red, arrowsize=12)
    end
    
    # Label faces
    for (i, face) in enumerate(fs)
        v1, v2, v3 = face[1], face[2], face[3]
        p1, p2, p3 = points_2d[v1], points_2d[v2], points_2d[v3]
        centroid = (p1 + p2 + p3) / 3
        text!(ax4, centroid, text="F$i", color=:blue, fontsize=14, align=(:center, :center))
    end
    
    scatter!(ax4, points_2d, color=:red, markersize=8)
    
    # Save figure
    output_dir = joinpath(@__DIR__, "..", "output", folder_name)
    mkpath(output_dir)
    output_file = joinpath(output_dir, "$(case_name).png")
    save(output_file, fig)
    println("\nVisualization saved to: $output_file")
    
    # Print detailed topology information
    println("\n" * "="^60)
    println("DETAILED TOPOLOGY INFORMATION")
    println("="^60)
    
    println("\nEdge Map:")
    for (edge, idx) in sort(collect(topology.edge_map), by=x->x[2])
        v1, v2 = edge
        adjacent_faces = topology.edge_to_faces[idx]
        println("  Edge $idx: ($v1, $v2) → Faces $adjacent_faces")
    end
    
    println("\nDual Adjacency:")
    for face_i in 1:topology.n_faces
        neighbors = topology.dual_adj[face_i]
        println("  Face $face_i:")
        for (face_j, edge) in neighbors
            edge_idx = topology.edge_map[edge]
            println("    → Face $face_j via edge $edge_idx $edge")
        end
    end
    
    println("\nReference Edges:")
    for face_idx in 1:topology.n_faces
        ref_edge = topology.face_reference_edges[face_idx]
        println("  Face $face_idx: reference edge $ref_edge")
    end
    
    println("\n" * "="^60)
    
    return fig
end

# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
