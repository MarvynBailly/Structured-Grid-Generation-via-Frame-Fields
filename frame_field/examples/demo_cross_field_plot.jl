using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

include("../src/meshio.jl")
include("../src/plotters.jl")

using LinearAlgebra

println("=== Cross Field Visualization Demo ===\n")

# Load mesh
mesh_path = joinpath(@__DIR__, "..", "..", "triangulations", "simple-square.msh")
println("Loading mesh: $mesh_path")
mesh = load_triangulation(mesh_path)

fs = faces(mesh)
n_faces = length(fs)
println("  Faces: $n_faces")

# Create synthetic theta values for demonstration
# Option 1: All zeros (aligned with x-axis)
println("\n[Demo 1] Uniform field (θ = 0 for all faces)")
theta_uniform = zeros(Float64, n_faces)

output_dir = joinpath(@__DIR__, "..", "output", "cross_field")
mkpath(output_dir)

savepath1 = joinpath(output_dir, "demo_uniform_field.png")
plot_cross_field(mesh, theta_uniform; scale=0.25, savepath=savepath1, field_color=:darkblue, field_width=1.5)

# Option 2: Gradient field (varies linearly)
println("\n[Demo 2] Gradient field (θ varies by face position)")
vs = coordinates(mesh)
V = zeros(Float64, length(vs), 3)
for (i, p) in enumerate(vs)
    V[i, :] = [p[1], p[2], p[3]]
end

theta_gradient = zeros(Float64, n_faces)
for (fi, f) in enumerate(fs)
    idxs = Tuple(f)
    # Centroid x-coordinate determines angle
    cx = mean(V[[idxs...], 1])
    theta_gradient[fi] = cx * π / 4  # Rotate based on x position
end

savepath2 = joinpath(output_dir, "demo_gradient_field.png")
plot_cross_field(mesh, theta_gradient; scale=0.25, savepath=savepath2, field_color=:darkred, field_width=1.5)

# Option 3: Random field
println("\n[Demo 3] Random field")
theta_random = rand(n_faces) .* π / 2

savepath3 = joinpath(output_dir, "demo_random_field.png")
plot_cross_field(mesh, theta_random; scale=0.25, savepath=savepath3, field_color=:darkgreen, field_width=1.5)

println("\n=== Visualization Demo Complete ===")
println("Generated visualizations:")
println("  1. $savepath1")
println("  2. $savepath2")
println("  3. $savepath3")
