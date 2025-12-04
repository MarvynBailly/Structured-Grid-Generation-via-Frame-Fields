"""
    simple_square.jl

Basic example: Generate frame field on a simple square mesh.
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

println("="^60)
println("Example: Simple Square Mesh Frame Field")
println("="^60)

# Generate mesh
println("\nGenerating 5×5 square mesh...")
mesh = generate_square_mesh(nx=5, ny=5)
println("  Vertices: $(length(coordinates(mesh)))")
println("  Faces: $(length(faces(mesh)))")

# Generate frame field (no constraints)
println("\nGenerating frame field...")
result = generate_frame_field(mesh, verbose=true)

# Print summary
print_frame_field_summary(result)

# Analyze singularities
println("Singularity Details:")
if isempty(result.singularities)
    println("  No singularities detected (perfectly smooth field)")
else
    for (v_idx, index) in result.singularities
        valence = 4 + round(Int, 4 * index)
        println("  Vertex $v_idx: index = $(round(index, digits=4)), expected valence = $valence")
    end
end

println("\n" * "="^60)
println("Example complete!")
println("="^60)
