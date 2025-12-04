"""
    visualize_sphere_mesh.jl

Load and visualize frame field on the sphere mesh.
"""

using GeometryBasics
using LinearAlgebra
using Printf

# Add src to load path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using FrameField

# Try to use CairoMakie if available
MAKIE_AVAILABLE = false
try
    using CairoMakie
    global MAKIE_AVAILABLE = true
catch
    @warn "CairoMakie not available. Install with: using Pkg; Pkg.add(\"CairoMakie\")"
end

"""
    load_gmsh_surface_mesh(filename::String)

Load a triangular surface mesh from a Gmsh .msh file (format 4.1).
Returns a GeometryBasics.Mesh with 3D points and triangular faces.
"""
function load_gmsh_surface_mesh(filename::String)
    if !isfile(filename)
        error("File not found: $filename")
    end
    
    println("Loading mesh from: $filename")
    
    vertices = Point3f[]
    faces = TriangleFace{Int}[]
    node_map = Dict{Int, Int}()  # Map from node ID to vertex index
    
    open(filename, "r") do file
        lines = readlines(file)
        i = 1
        
        # Find nodes section
        while i <= length(lines)
            if startswith(lines[i], "\$Nodes")
                i += 1
                # Parse header: numEntityBlocks numNodes minNodeTag maxNodeTag
                header = split(lines[i])
                num_nodes = parse(Int, header[2])
                println("  Reading $num_nodes nodes...")
                
                i += 1
                
                # Read entity blocks
                while i <= length(lines) && !startswith(lines[i], "\$EndNodes")
                    # Entity block header
                    if !startswith(lines[i], "\$")
                        parts = split(lines[i])
                        if length(parts) >= 4
                            num_nodes_in_block = parse(Int, parts[4])
                            i += 1
                            
                            # Read node IDs
                            node_ids = Int[]
                            for _ in 1:num_nodes_in_block
                                push!(node_ids, parse(Int, lines[i]))
                                i += 1
                            end
                            
                            # Read coordinates
                            for node_id in node_ids
                                coords = split(lines[i])
                                x = parse(Float64, coords[1])
                                y = parse(Float64, coords[2])
                                z = parse(Float64, coords[3])
                                push!(vertices, Point3f(x, y, z))
                                node_map[node_id] = length(vertices)
                                i += 1
                            end
                        else
                            i += 1
                        end
                    else
                        i += 1
                    end
                end
                break
            end
            i += 1
        end
        
        # Find elements section
        i = 1
        while i <= length(lines)
            if startswith(lines[i], "\$Elements")
                i += 1
                # Parse header
                header = split(lines[i])
                num_elements = parse(Int, header[2])
                println("  Reading elements...")
                
                i += 1
                
                # Read entity blocks
                while i <= length(lines) && !startswith(lines[i], "\$EndElements")
                    if !startswith(lines[i], "\$")
                        parts = split(lines[i])
                        if length(parts) >= 4
                            entity_dim = parse(Int, parts[1])
                            element_type = parse(Int, parts[3])
                            num_elements_in_block = parse(Int, parts[4])
                            i += 1
                            
                            # Only read 2D triangular elements (type 2)
                            if entity_dim == 2 && element_type == 2
                                for _ in 1:num_elements_in_block
                                    parts = split(lines[i])
                                    # Format: element-id node1 node2 node3
                                    n1 = parse(Int, parts[2])
                                    n2 = parse(Int, parts[3])
                                    n3 = parse(Int, parts[4])
                                    
                                    # Map to vertex indices
                                    if haskey(node_map, n1) && haskey(node_map, n2) && haskey(node_map, n3)
                                        push!(faces, TriangleFace(node_map[n1], node_map[n2], node_map[n3]))
                                    end
                                    i += 1
                                end
                            else
                                # Skip non-triangle elements
                                i += num_elements_in_block
                            end
                        else
                            i += 1
                        end
                    else
                        i += 1
                    end
                end
                break
            end
            i += 1
        end
    end
    
    println("  Loaded: $(length(vertices)) vertices, $(length(faces)) faces")
    
    return GeometryBasics.Mesh(vertices, faces)
end

"""
    plot_sphere_frame_field(mesh, result; arrow_scale=5.0, show_all_dirs=false)

Visualize the frame field on a 3D surface mesh using 3D projection.
"""
function plot_sphere_frame_field(mesh, result; arrow_scale=5.0, show_all_dirs=false)
    if !MAKIE_AVAILABLE
        @error "CairoMakie is required for visualization"
        return
    end
    
    fig = Figure(size=(1400, 1200))
    
    # Get mesh data
    positions = coordinates(mesh)
    fs = faces(mesh)
    
    # Create 3 views: two 2D projections and one 3D
    
    # XY projection (top view)
    ax1 = Axis(fig[1, 1], 
               aspect=DataAspect(),
               xlabel="x", 
               ylabel="y",
               title="XY Projection (Top View)")
    
    # XZ projection (front view)
    ax2 = Axis(fig[1, 2], 
               aspect=DataAspect(),
               xlabel="x", 
               ylabel="z",
               title="XZ Projection (Front View)")
    
    # YZ projection (side view)
    ax3 = Axis(fig[2, 1], 
               aspect=DataAspect(),
               xlabel="y", 
               ylabel="z",
               title="YZ Projection (Side View)")
    
    # Plot mesh edges for all views
    for face in fs
        v1, v2, v3 = face
        p1, p2, p3 = positions[v1], positions[v2], positions[v3]
        
        # XY projection
        lines!(ax1, [Point2f(p1[1], p1[2]), Point2f(p2[1], p2[2])], color=:black, alpha=0.15, linewidth=0.5)
        lines!(ax1, [Point2f(p2[1], p2[2]), Point2f(p3[1], p3[2])], color=:black, alpha=0.15, linewidth=0.5)
        lines!(ax1, [Point2f(p3[1], p3[2]), Point2f(p1[1], p1[2])], color=:black, alpha=0.15, linewidth=0.5)
        
        # XZ projection
        lines!(ax2, [Point2f(p1[1], p1[3]), Point2f(p2[1], p2[3])], color=:black, alpha=0.15, linewidth=0.5)
        lines!(ax2, [Point2f(p2[1], p2[3]), Point2f(p3[1], p3[3])], color=:black, alpha=0.15, linewidth=0.5)
        lines!(ax2, [Point2f(p3[1], p3[3]), Point2f(p1[1], p1[3])], color=:black, alpha=0.15, linewidth=0.5)
        
        # YZ projection
        lines!(ax3, [Point2f(p1[2], p1[3]), Point2f(p2[2], p2[3])], color=:black, alpha=0.15, linewidth=0.5)
        lines!(ax3, [Point2f(p2[2], p2[3]), Point2f(p3[2], p3[3])], color=:black, alpha=0.15, linewidth=0.5)
        lines!(ax3, [Point2f(p3[2], p3[3]), Point2f(p1[2], p1[3])], color=:black, alpha=0.15, linewidth=0.5)
    end
    
    # Get frame directions
    directions = get_frame_directions(mesh, result)
    
    # For each face, compute tangent frame and project field directions
    for (face_idx, face) in enumerate(fs)
        v1, v2, v3 = face
        p1, p2, p3 = positions[v1], positions[v2], positions[v3]
        
        # Compute face center
        center = (p1 .+ p2 .+ p3) ./ 3
        
        # Compute face normal
        edge1 = p2 .- p1
        edge2 = p3 .- p1
        normal = normalize(cross(edge1, edge2))
        
        # Create local tangent frame
        # Use edge1 as first tangent direction (projected to tangent plane)
        tangent1 = normalize(edge1 - dot(edge1, normal) * normal)
        tangent2 = normalize(cross(normal, tangent1))
        
        # Get angle for this face
        θ = result.angles[face_idx]
        
        # Color based on angle
        color_val = mod(θ, π/2) / (π/2)
        
        # Compute field directions in tangent plane
        n_dirs = show_all_dirs ? 4 : 2
        for i in 0:(n_dirs-1)
            angle = θ + i * π/2
            
            # Direction in tangent plane
            dir_3d = cos(angle) * tangent1 + sin(angle) * tangent2
            dir_3d = normalize(dir_3d) * arrow_scale
            
            # Project to each view
            # XY projection
            arrow_starts_xy = [Point2f(center[1], center[2])]
            arrow_vecs_xy = [Vec2f(dir_3d[1], dir_3d[2])]
            arrows!(ax1, arrow_starts_xy, arrow_vecs_xy,
                   color=[color_val], colormap=:hsv, colorrange=(0, 1),
                   arrowsize=arrow_scale*0.15, linewidth=1.0, alpha=0.7)
            
            # XZ projection
            arrow_starts_xz = [Point2f(center[1], center[3])]
            arrow_vecs_xz = [Vec2f(dir_3d[1], dir_3d[3])]
            arrows!(ax2, arrow_starts_xz, arrow_vecs_xz,
                   color=[color_val], colormap=:hsv, colorrange=(0, 1),
                   arrowsize=arrow_scale*0.15, linewidth=1.0, alpha=0.7)
            
            # YZ projection
            arrow_starts_yz = [Point2f(center[2], center[3])]
            arrow_vecs_yz = [Vec2f(dir_3d[2], dir_3d[3])]
            arrows!(ax3, arrow_starts_yz, arrow_vecs_yz,
                   color=[color_val], colormap=:hsv, colorrange=(0, 1),
                   arrowsize=arrow_scale*0.15, linewidth=1.0, alpha=0.7)
        end
    end
    
    # Statistics text
    stats_ax = Axis(fig[2, 2], aspect=DataAspect())
    hidedecorations!(stats_ax)
    hidespines!(stats_ax)
    
    stats_text = """
    Frame Field Statistics
    
    Mesh:
      Vertices: $(length(positions))
      Faces: $(length(fs))
    
    Optimization:
      Converged: $(result.converged)
      Iterations: $(result.iterations)
      Energy: $(round(result.energy, sigdigits=4))
    
    Angles:
      Range: [$(round(minimum(result.angles), digits=3)), 
              $(round(maximum(result.angles), digits=3))]
    
    Singularities:
      Count: $(length(result.singularities))
      Index range: $(length(result.singularities) > 0 ? 
                     "[$(round(minimum(last.(result.singularities)), digits=3)), " *
                     "$(round(maximum(last.(result.singularities)), digits=3))]" : "N/A")
    """
    
    text!(stats_ax, 0.1, 0.9, text=stats_text, align=(:left, :top), 
          fontsize=14, font=:regular, space=:relative)
    
    # Add colorbar
    Colorbar(fig[1:2, 3], 
             colormap=:hsv,
             limits=(0, π/2),
             label="Angle (mod π/2)")
    
    return fig
end

# Main execution
println("="^60)
println("Sphere Mesh Frame Field Visualization")
println("="^60)

# Load mesh
mesh_file = joinpath(@__DIR__, "..", "..", "triangulations", "300_polygon_sphere_100mm.msh")
mesh = load_gmsh_surface_mesh(mesh_file)

# Generate frame field (no constraints)
println("\nGenerating frame field...")
result = generate_frame_field(mesh, verbose=true)

# Print summary
print_frame_field_summary(result)

# Visualize
println("\nVisualizing frame field...")
if MAKIE_AVAILABLE
    fig = plot_sphere_frame_field(mesh, result, arrow_scale=3.0, show_all_dirs=false)
    output_file = "sphere_frame_field_output.png"
    save(output_file, fig, px_per_unit=2)
    println("Saved visualization to: $output_file")
else
    println("Install CairoMakie to see visualization:")
    println("  using Pkg; Pkg.add(\"CairoMakie\")")
end

println("\n" * "="^60)
println("Visualization complete!")
println("="^60)
