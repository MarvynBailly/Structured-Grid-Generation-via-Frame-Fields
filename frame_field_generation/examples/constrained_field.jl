"""
    constrained_field.jl

Example: Generate frame field with boundary constraints.
"""

using GeometryBasics

# Add src to load path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using FrameField

"""
    generate_square_mesh(; nx=5, ny=5)

Generate a triangular mesh of a unit square [0,1]×[0,1].
"""
function generate_square_mesh(; nx=5, ny=5)
    vs = Point3f[]
    for j in 0:ny
        for i in 0:nx
            push!(vs, Point3f(i/nx, j/ny, 0))
        end
    end
    
    fs = TriangleFace{Int}[]
    row_stride = nx + 1
    
    for j in 1:ny
        for i in 1:nx
            v_bl = (j-1) * row_stride + i
            v_br = v_bl + 1
            v_tl = v_bl + row_stride
            v_tr = v_tl + 1
            
            push!(fs, TriangleFace(v_bl, v_br, v_tr))
            push!(fs, TriangleFace(v_bl, v_tr, v_tl))
        end
    end
    
    return Mesh(vs, fs)
end

println("="^60)
println("Example: Constrained Frame Field")
println("="^60)

# Generate mesh
println("\nGenerating 6×6 square mesh...")
mesh = generate_square_mesh(nx=6, ny=6)
n_faces = length(faces(mesh))
println("  Vertices: $(length(coordinates(mesh)))")
println("  Faces: $n_faces")

# Set up constraints on boundary faces
# Constrain corner faces to specific orientations
println("\nSetting up boundary constraints...")
constrained_faces = [1, 12, n_faces-11, n_faces]  # Four corners
constrained_angles = [0.0, π/4, π/2, 3π/4]  # Different orientations

println("  Constrained faces: $constrained_faces")
println("  Angles: $(round.(constrained_angles .* 180/π, digits=1))°")

# Generate frame field with constraints
println("\nGenerating constrained frame field...")
result = generate_frame_field(
    mesh,
    constrained_faces=constrained_faces,
    constrained_angles=constrained_angles,
    verbose=true
)

# Verify constraints
println("\nVerifying constraints:")
for (i, face_idx) in enumerate(constrained_faces)
    actual = result.angles[face_idx]
    expected = constrained_angles[i]
    error = abs(actual - expected)
    println("  Face $face_idx: expected=$(round(expected*180/π, digits=1))°, " *
            "actual=$(round(actual*180/π, digits=1))°, error=$(round(error, digits=8))")
end

# Print summary
print_frame_field_summary(result)

println("\n" * "="^60)
println("Example complete!")
println("="^60)
