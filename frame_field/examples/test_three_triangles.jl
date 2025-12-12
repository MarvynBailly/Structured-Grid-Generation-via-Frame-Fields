"""
    test_three_triangles.jl

Test the three-triangle mesh to verify topology and visualize.
"""

using Pkg
Pkg.activate(".")

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using GeometryBasics
using CairoMakie

include("../src/meshio.jl")
include("../src/mesh_types.jl")

function main()
    println("Loading three-triangles mesh...")
    mesh_file = joinpath(@__DIR__, "..", "..", "triangulations", "three-triangles.msh")
    mesh = load_triangulation(mesh_file)
    
    println("Building mesh topology...")
    topology = build_mesh_topology(mesh)
    
    println("\n" * "="^60)
    println("MESH INFORMATION")
    println("="^60)
    
    vs = coordinates(mesh)
    fs = faces(mesh)
    
    println("\nVertices ($(length(vs)) total):")
    for (i, v) in enumerate(vs)
        println("  v$i: ($(v[1]), $(v[2]), $(v[3]))")
    end
    
    println("\nFaces ($(length(fs)) total):")
    for (i, f) in enumerate(fs)
        println("  f$i: vertices $(f[1]), $(f[2]), $(f[3])")
        ref_edge = topology.face_reference_edges[i]
        println("      reference edge: $ref_edge")
    end
    
    println("\nDual Adjacency:")
    for face_i in 1:topology.n_faces
        neighbors = topology.dual_adj[face_i]
        if !isempty(neighbors)
            println("  Face $face_i → ", [n[1] for n in neighbors])
        end
    end
    
    println("\nEdges ($(topology.n_edges) total):")
    for (edge, idx) in sort(collect(topology.edge_map), by=x->x[2])
        faces = topology.edge_to_faces[idx]
        println("  Edge $idx: vertices $edge → faces $faces")
    end
    
    # Create frame field
    println("\nCreating frame field...")
    field = CrossField(topology)
    
    println("\nTransport angles:")
    for face_i in 1:topology.n_faces
        for (face_j, edge) in topology.dual_adj[face_i]
            if face_i < face_j
                kappa_ij = get_transport_angle(field, face_i, face_j)
                println("  κ_($face_i,$face_j) = $(round(rad2deg(kappa_ij), digits=2))°")
            end
        end
    end
    
    # Visualize
    println("\nCreating visualization...")
    fig = Figure(size=(1200, 400))
    ax = Axis(fig[1, 1], aspect=DataAspect(), title="Three Triangles")
    
    # Draw triangles
    for (i, f) in enumerate(fs)
        v1, v2, v3 = vs[f[1]], vs[f[2]], vs[f[3]]
        poly!(ax, [Point2f(v1[1], v1[2]), Point2f(v2[1], v2[2]), Point2f(v3[1], v3[2])],
              color=(:blue, 0.1), strokecolor=:blue, strokewidth=2)
        
        # Face label
        centroid = (Point2f(v1[1], v1[2]) + Point2f(v2[1], v2[2]) + Point2f(v3[1], v3[2])) / 3
        text!(ax, centroid, text="$i", color=:red, fontsize=20, align=(:center, :center))
        
        # Draw reference edge
        ref_edge = topology.face_reference_edges[i]
        p1 = vs[ref_edge[1]]
        p2 = vs[ref_edge[2]]
        lines!(ax, [Point2f(p1[1], p1[2]), Point2f(p2[1], p2[2])], 
               color=:orange, linewidth=4)
    end
    
    # Draw vertices
    for (i, v) in enumerate(vs)
        scatter!(ax, [Point2f(v[1], v[2])], color=:black, markersize=10)
        text!(ax, Point2f(v[1], v[2]) .+ Point2f(0, 0.1), 
              text="v$i", color=:black, fontsize=12, align=(:center, :bottom))
    end
    
    output_file = joinpath(@__DIR__, "..", "output", "three_triangles_test.png")
    save(output_file, fig)
    println("\nVisualization saved to: $output_file")
    
    println("\n" * "="^60)
end

main()
