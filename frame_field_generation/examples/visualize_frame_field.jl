"""
    visualize_frame_field.jl

Visualize frame field results with arrows showing field directions.
Uses CairoMakie for 2D visualization.
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
    generate_square_mesh(; nx=5, ny=5)

Generate a triangular mesh of a unit square [0,1]×[0,1].
"""
function generate_square_mesh(; nx=5, ny=5)
    # Generate vertices
    vs = Point3f[]
    for j in 0:ny
        for i in 0:nx
            push!(vs, Point3f(i/nx, j/ny, 0))
        end
    end
    
    # Generate faces
    fs = TriangleFace{Int}[]
    row_stride = nx + 1
    
    for j in 1:ny
        for i in 1:nx
            v_bl = (j-1) * row_stride + i
            v_br = v_bl + 1
            v_tl = v_bl + row_stride
            v_tr = v_tl + 1
            
            # Two triangles per quad
            push!(fs, TriangleFace(v_bl, v_br, v_tr))
            push!(fs, TriangleFace(v_bl, v_tr, v_tl))
        end
    end
    
    return GeometryBasics.Mesh(vs, fs)
end

"""
    plot_frame_field_2d(mesh, result; arrow_scale=0.08, show_all_dirs=true)

Visualize the frame field with CairoMakie.
"""
function plot_frame_field_2d(mesh, result; arrow_scale=0.08, show_all_dirs=true)
    if !MAKIE_AVAILABLE
        @error "CairoMakie is required for visualization. Install with: using Pkg; Pkg.add(\"CairoMakie\")"
        return
    end
    
    fig = Figure(size=(1200, 1000))
    ax = Axis(fig[1, 1], 
              aspect=DataAspect(),
              xlabel="x", 
              ylabel="y",
              title="Frame Field Visualization\nEnergy: $(round(result.energy, sigdigits=4)), Singularities: $(length(result.singularities))")
    
    # Get mesh data
    positions = coordinates(mesh)
    fs = faces(mesh)
    
    # Plot mesh edges
    for face in fs
        v1, v2, v3 = face
        p1, p2, p3 = positions[v1], positions[v2], positions[v3]
        
        # Draw triangle edges
        lines!(ax, [Point2f(p1[1], p1[2]), Point2f(p2[1], p2[2])], color=:black, alpha=0.2, linewidth=0.5)
        lines!(ax, [Point2f(p2[1], p2[2]), Point2f(p3[1], p3[2])], color=:black, alpha=0.2, linewidth=0.5)
        lines!(ax, [Point2f(p3[1], p3[2]), Point2f(p1[1], p1[2])], color=:black, alpha=0.2, linewidth=0.5)
    end
    
    # Get frame directions
    directions = get_frame_directions(mesh, result)
    
    # Collect arrow data
    arrow_starts = Point2f[]
    arrow_vecs = Vec2f[]
    arrow_colors = Float64[]
    
    # Plot frame directions at each face center
    for (face_idx, face) in enumerate(fs)
        v1, v2, v3 = face
        p1, p2, p3 = positions[v1], positions[v2], positions[v3]
        
        # Compute face center
        center = (p1 .+ p2 .+ p3) ./ 3
        cx, cy = center[1], center[2]
        
        # Get angle for this face
        θ = result.angles[face_idx]
        
        # Get all 4 directions (90° rotational symmetry)
        dirs = directions[face_idx]
        
        # Color based on angle (mod π/2 for 4-way symmetry)
        color_val = mod(θ, π/2) / (π/2)
        
        # Plot frame directions
        n_dirs = show_all_dirs ? 4 : 2
        for i in 1:n_dirs
            d = dirs[i]
            dx, dy = d[1] * arrow_scale, d[2] * arrow_scale
            
            push!(arrow_starts, Point2f(cx, cy))
            push!(arrow_vecs, Vec2f(dx, dy))
            push!(arrow_colors, color_val)
        end
    end
    
    # Draw all arrows at once
    arrows!(ax, arrow_starts, arrow_vecs, 
            color=arrow_colors, 
            colormap=:hsv,
            arrowsize=arrow_scale*0.2,
            linewidth=1.5,
            alpha=0.7)
    
    # Plot singularities
    if !isempty(result.singularities)
        sing_pos = Point2f[]
        sing_colors = Symbol[]
        sing_sizes = Float64[]
        
        for (v_idx, index) in result.singularities
            pos = positions[v_idx]
            push!(sing_pos, Point2f(pos[1], pos[2]))
            push!(sing_colors, index > 0 ? :red : :blue)
            push!(sing_sizes, abs(index) * 40 + 15)
        end
        
        scatter!(ax, sing_pos, 
                color=sing_colors,
                marker=:star5,
                markersize=sing_sizes,
                strokewidth=2,
                strokecolor=:black)
        
        # Add labels for singularities
        for (v_idx, index) in result.singularities
            pos = positions[v_idx]
            text!(ax, @sprintf("%.2f", index),
                  position=Point2f(pos[1], pos[2] + 0.05),
                  align=(:center, :bottom),
                  fontsize=10)
        end
    end
    
    # Add colorbar
    Colorbar(fig[1, 2], 
             colormap=:hsv,
             limits=(0, π/2),
             label="Angle (mod π/2)")
    
    display(fig)
    
    return fig
end

"""
    save_frame_field_plot(mesh, result, filename; kwargs...)

Save frame field visualization to file.
"""
function save_frame_field_plot(mesh, result, filename; kwargs...)
    if !MAKIE_AVAILABLE
        @error "CairoMakie is required for visualization"
        return
    end
    
    fig = plot_frame_field_2d(mesh, result; kwargs...)
    save(filename, fig, px_per_unit=2)
    println("Saved visualization to: $filename")
    return fig
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    println("="^60)
    println("Frame Field Visualization")
    println("="^60)
    
    # Generate mesh
    println("\nGenerating 8×8 square mesh...")
    mesh = generate_square_mesh(nx=8, ny=8)
    println("  Vertices: $(length(coordinates(mesh)))")
    println("  Faces: $(length(faces(mesh)))")
    
    # Generate frame field (no constraints)
    println("\nGenerating frame field...")
    result = generate_frame_field(mesh, verbose=true)
    
    # Print summary
    print_frame_field_summary(result)
    
    # Visualize
    println("\nVisualizing frame field...")
    if MAKIE_AVAILABLE
        fig = plot_frame_field_2d(mesh, result, arrow_scale=0.06, show_all_dirs=true)
        save("frame_field_output.png", fig, px_per_unit=2)
        println("Saved visualization to: frame_field_output.png")
    else
        println("Install CairoMakie to see visualization:")
        println("  using Pkg; Pkg.add(\"CairoMakie\")")
    end
    
    println("\n" * "="^60)
    println("Visualization complete!")
    println("="^60)
end
