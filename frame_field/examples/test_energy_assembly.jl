include("../src/meshio.jl")
include("../src/dijkstra_forest.jl")
include("../src/cross_field_energy.jl")

using LinearAlgebra
using SparseArrays
using CairoMakie

println("=== Energy Matrix Assembly Test ===\n")

# casename = "simple-square"
casename = "300_polygon_sphere_100mm"


# Load mesh (path relative to root project directory)
mesh_path = joinpath(@__DIR__, "..", "..", "triangulations", "$casename.msh")
println("Loading mesh: $mesh_path")
mesh = load_triangulation(mesh_path)

fs = faces(mesh)
n_faces = length(fs)
println("  Faces: $n_faces")

# Set up constraints (optional - empty for now)
constrained_faces = Int[1, 5]  # Example: constrain faces 1 and 5
constrained_angles = Float64[0.0, π/4]  # Example angles in radians

# Compute spanning forest
println("\nComputing spanning forest...")
potential_fixed_edges = compute_spanning_forest(mesh; constrained_faces=constrained_faces)
println("  Potential fixed edges: $(length(potential_fixed_edges))")

# Select suitable fixed edges (one per face)
println("\nSelecting fixed edges...")
fixed_edges_per_face = fix_suitable_edges(mesh, potential_fixed_edges)
println("  Fixed edges selected: $(length(fixed_edges_per_face))")

# Assemble system matrix
println("\nAssembling system matrix...")
A, b, var_to_p, var_to_theta = assemble_system_matrix(
    mesh, 
    potential_fixed_edges, 
    fixed_edges_per_face,
    constrained_faces,
    constrained_angles,
    true  # debug mode
)

n_vars = length(b)
n_free_p = length(var_to_p)
n_free_theta = length(var_to_theta)

println("  System size: $(size(A, 1)) x $(size(A, 2))")
println("  Total variables: $n_vars")
println("  Period jump variables (integers): $n_free_p")
println("  Angle variables (reals): $n_free_theta")
println("  Non-zeros in A: $(nnz(A))")

# Verify variable ordering
println("\nVariable ordering verification:")
println("  First variable type: ", haskey(var_to_p, 1) ? "period jump" : "angle")
println("  Last variable type: ", haskey(var_to_theta, n_vars) ? "angle" : "period jump")

# Check matrix properties
println("\nMatrix properties:")
println("  A is symmetric: ", issymmetric(A))
println("  Rank of A: ", rank(Matrix(A)))
println("  Condition number: ", cond(Matrix(A)))

# Display first few variable mappings
println("\nFirst 5 period jump variables:")
for i in 1:min(5, n_free_p)
    edge = var_to_p[i]
    println("  Variable $i → Edge $edge")
end

println("\nFirst 5 angle variables:")
for i in 1:min(5, n_free_theta)
    var_idx = n_free_p + i
    face_idx = var_to_theta[var_idx]
    println("  Variable $var_idx → Face $face_idx")
end

# Show RHS statistics
println("\nRight-hand side statistics:")
println("  Norm of b: $(norm(b))")
println("  Min value: $(minimum(b))")
println("  Max value: $(maximum(b))")
println("  Mean value: $(sum(b) / length(b))")

println("\n=== Assembly Complete ===")
println("\nThe system Ax = b is ready for the greedy MIP solver.")
println("Variable ordering: x = [p_1, ..., p_$n_free_p, θ_1, ..., θ_$n_free_theta]")

# Visualize the system matrices
println("\n=== Creating Visualizations ===")

# Create output directory (relative to script location)
output_dir = joinpath(@__DIR__, "..", "output", "energy_system")
mkpath(output_dir)
println("  Output directory: $output_dir")

# Convert sparse matrix to dense for visualization
A_dense = Matrix(A)

# Create figure with heatmaps
fig = Figure(size=(1400, 600))

# Heatmap of A matrix
ax1 = Axis(fig[1, 1], 
    title="System Matrix A ($(size(A, 1)) × $(size(A, 2)))",
    xlabel="Variable Index",
    ylabel="Equation Index",
    yreversed=true,  # Reverse y-axis so row 1 is at top
    aspect=DataAspect()
)
A_max = maximum(abs.(A_dense))
hm1 = heatmap!(ax1, 1:size(A_dense, 2), 1:size(A_dense, 1), A_dense, colormap=:RdBu, colorrange=(-A_max, A_max))
Colorbar(fig[1, 2], hm1, label="Matrix Value")

# Heatmap of b vector (as column)
ax2 = Axis(fig[1, 3],
    title="RHS Vector b ($(length(b)) × 1)",
    xlabel="Column",
    ylabel="Equation Index",
    yreversed=true,  # Reverse y-axis so row 1 is at top
    aspect=DataAspect()
)
b_col = reshape(b, :, 1)  # Make it a column for visualization
b_max = maximum(abs.(b))
# Handle case where b is all zeros
if b_max == 0.0
    b_max = 1.0
end
hm2 = heatmap!(ax2, [1], 1:length(b), b_col, colormap=:RdBu, colorrange=(-b_max, b_max))
Colorbar(fig[1, 4], hm2, label="Vector Value")

# Save figure
savepath = joinpath(output_dir, "system_matrices_$(casename).png")
save(savepath, fig)
println("  Saved matrix visualization to: $savepath")

println("\n=== Visualization Complete ===")
