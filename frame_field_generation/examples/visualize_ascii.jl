"""
    visualize_ascii.jl

Simple ASCII visualization of frame field showing angles and singularities.
"""

using GeometryBasics
using LinearAlgebra
using Printf

# Add src to load path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using FrameField

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
    
    return Mesh(vs, fs)
end

"""
    visualize_frame_field_ascii(mesh, result)

Print ASCII visualization of frame field.
"""
function visualize_frame_field_ascii(mesh, result)
    positions = coordinates(mesh)
    fs = faces(mesh)
    
    # Find grid bounds
    xs = [p[1] for p in positions]
    ys = [p[2] for p in positions]
    
    nx = length(unique(xs)) - 1
    ny = length(unique(ys)) - 1
    
    println("\nFrame Field Visualization (ASCII)")
    println("="^60)
    println("Grid: $(nx)×$(ny), Faces: $(length(fs))")
    println()
    
    # Create grid visualization
    # Show angle at each quad center using ASCII characters
    angle_chars = ['─', '/', '│', '\\']  # 0°, 45°, 90°, 135°
    
    # Build a grid of face centers and their angles
    face_grid = Dict{Tuple{Int,Int}, Int}()
    
    for (face_idx, face) in enumerate(fs)
        v1, v2, v3 = face
        p1, p2, p3 = positions[v1], positions[v2], positions[v3]
        
        # Compute face center
        center = (p1 .+ p2 .+ p3) ./ 3
        
        # Map to grid position
        gi = round(Int, center[1] * nx + 0.5)
        gj = round(Int, center[2] * ny + 0.5)
        
        face_grid[(gi, gj)] = face_idx
    end
    
    # Print grid from top to bottom
    println("Field directions (angles):")
    println()
    for j in ny:-1:1
        print("  ")
        for i in 1:nx
            if haskey(face_grid, (i, j))
                face_idx = face_grid[(i, j)]
                θ = result.angles[face_idx]
                
                # Map angle to character (mod π/2 for 4-fold symmetry)
                angle_norm = mod(θ, π/2)
                char_idx = round(Int, angle_norm / (π/2) * 3) + 1
                char_idx = clamp(char_idx, 1, 4)
                
                print(angle_chars[char_idx], " ")
            else
                print("  ")
            end
        end
        println()
    end
    
    println()
    println("Legend: ─ = 0°, / = 45°, │ = 90°, \\ = 135° (mod π/2)")
    println()
    
    # Show singularities
    if !isempty(result.singularities)
        println("\nSingularities:")
        println("-"^60)
        
        # Group by index value
        index_groups = Dict{Float64, Vector{Int}}()
        for (v_idx, index) in result.singularities
            idx_rounded = round(index, digits=3)
            if !haskey(index_groups, idx_rounded)
                index_groups[idx_rounded] = []
            end
            push!(index_groups[idx_rounded], v_idx)
        end
        
        for (index, vertices) in sort(collect(index_groups))
            valence = 4 + round(Int, 4 * index)
            println(@sprintf("  Index %.3f (valence %d): %d vertices", 
                           index, valence, length(vertices)))
        end
    else
        println("\nNo singularities (perfectly smooth field)")
    end
    
    println()
    println("="^60)
end

# Main execution
println("="^60)
println("Frame Field ASCII Visualization")
println("="^60)

# Generate mesh
println("\nGenerating 6×6 square mesh...")
mesh = generate_square_mesh(nx=6, ny=6)
println("  Vertices: $(length(coordinates(mesh)))")
println("  Faces: $(length(faces(mesh)))")

# Generate frame field (no constraints)
println("\nGenerating frame field...")
result = generate_frame_field(mesh, verbose=false)

# Print summary
print_frame_field_summary(result)

# ASCII Visualization
visualize_frame_field_ascii(mesh, result)

println("\nFor graphical visualization, run:")
println("  julia --project=. examples/visualize_frame_field.jl")
println("  (requires PyPlot)")

println("\n" * "="^60)
println("Complete!")
println("="^60)
