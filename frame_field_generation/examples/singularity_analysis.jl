"""
    singularity_analysis.jl

Example: Analyze singularity placement in frame fields.
"""

using GeometryBasics
using Statistics

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
println("Example: Singularity Analysis")
println("="^60)

# Test different mesh resolutions
mesh_sizes = [3, 5, 7, 10]

println("\nAnalyzing singularity placement across different mesh resolutions:\n")

for nx in mesh_sizes
    println("-"^60)
    println("Mesh size: $(nx)×$(nx)")
    
    mesh = generate_square_mesh(nx=nx, ny=nx)
    n_vertices = length(coordinates(mesh))
    n_faces = length(faces(mesh))
    
    println("  Vertices: $n_vertices, Faces: $n_faces")
    
    # Generate frame field
    result = generate_frame_field(mesh, verbose=false)
    
    n_sing = length(result.singularities)
    println("  Singularities: $n_sing")
    
    if n_sing > 0
        indices = [idx for (_, idx) in result.singularities]
        println("  Index range: [$(minimum(indices)), $(maximum(indices))]")
        println("  Mean index: $(round(mean(indices), digits=4))")
        
        # Count by valence
        valences = [4 + round(Int, 4 * idx) for (_, idx) in result.singularities]
        val_counts = Dict{Int, Int}()
        for v in valences
            val_counts[v] = get(val_counts, v, 0) + 1
        end
        
        println("  Valence distribution:")
        for (val, count) in sort(collect(val_counts))
            println("    Valence $val: $count vertices")
        end
    else
        println("  No singularities (smooth field)")
    end
    
    println("  Energy: $(round(result.energy, digits=6))")
end

println("\n" * "="^60)
println("\nObservations:")
println("- For flat square meshes, singularities may not be necessary")
println("- Energy generally increases with mesh resolution")
println("- Singularity placement is automatic based on geometry")
println("- For general surfaces, singularities compensate for Gaussian curvature")
println("="^60)

println("\nExample complete!")
