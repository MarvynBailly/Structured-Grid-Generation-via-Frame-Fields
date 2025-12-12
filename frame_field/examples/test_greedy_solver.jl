include("../src/meshio.jl")
include("../src/dijkstra_forest.jl")
include("../src/cross_field_energy.jl")
include("../src/greedy_solver.jl")
include("../src/singularities.jl")

using LinearAlgebra
using SparseArrays

println("=== Greedy Mixed-Integer Solver Test ===\n")

# Load mesh (path relative to root project directory)
# mesh_path = joinpath(@__DIR__, "..", "..", "triangulations", "simple-square.msh")
mesh_path = joinpath(@__DIR__, "..", "..", "triangulations", "300_polygon_sphere_100mm.msh")
# mesh_path = joinpath(@__DIR__, "..", "..", "triangulations", "hole.msh")
println("Loading mesh: $mesh_path")
mesh = load_triangulation(mesh_path)

fs = faces(mesh)
n_faces = length(fs)
println("  Faces: $n_faces")

# Identify boundary faces and compute boundary alignment constraints
println("\nIdentifying boundary constraints...")
boundary_faces, boundary_angles = identify_boundary_faces_and_angles(mesh)
println("  Boundary faces found: $(length(boundary_faces))")

# Set up additional user constraints (optional)
# Start with NO additional constraints - just use boundary
user_constrained_faces = Int[]
user_constrained_angles = Float64[]

# For closed surfaces (like spheres) without boundary, we need at least one constraint
# to break the rotational symmetry. MIQ typically uses 1-2 constraints.
if isempty(boundary_faces)
    println("\nNo boundary detected - adding constraint to break symmetry...")
    user_constrained_faces = Int[1]  # Constrain first face
    user_constrained_angles = Float64[0.0]  # Align with x-axis
end

# Optionally add some interior constraints for testing
# n_fixed = 2
# user_constrained_faces = rand(1:n_faces, n_fixed) 
# user_constrained_angles = rand(n_fixed) .* (π/2)  # Random angles in [0, π/2]

# Combine boundary constraints with user constraints
constrained_faces = vcat(boundary_faces, user_constrained_faces)
constrained_angles = vcat(boundary_angles, user_constrained_angles)

println("\nTotal constrained faces: $(length(constrained_faces))")
println("  - Boundary faces: $(length(boundary_faces))")
println("  - User constrained: $(length(user_constrained_faces))")

# Compute spanning forest
println("\nComputing spanning forest...")
potential_fixed_edges = compute_spanning_forest(mesh; constrained_faces=constrained_faces)
println("  Potential fixed edges: $(length(potential_fixed_edges))")

# Select suitable fixed edges (one per face)
println("\nSelecting fixed edges...")
fixed_edges_per_face = fix_suitable_edges(mesh, potential_fixed_edges)
println("  Fixed edges selected: $(length(fixed_edges_per_face))")

# Use ALL spanning forest edges as fixed edges (like reference implementation)
# The reference fixes the spanning forest to eliminate rotational freedom
# fixed_edges_to_use = potential_fixed_edges
# Or use just one edge per face:
fixed_edges_to_use = Set(values(fixed_edges_per_face))
println("  Actually fixing edges: $(length(fixed_edges_to_use))")

# Assemble FREE variable system (for greedy solver)
println("\nAssembling FREE variable system...")
A, b, var_to_p, var_to_theta = assemble_free_system(
    mesh, 
    fixed_edges_to_use,  # Use selected fixed edges
    fixed_edges_per_face,  # Still needed for face-edge mapping
    constrained_faces,
    constrained_angles,
    true  # debug mode
)

n_vars = length(b)
n_p = length(var_to_p)
n_theta = length(var_to_theta)

println("\n  System size: $(size(A, 1)) × $(size(A, 2))")
println("  Total variables: $n_vars")
println("  Period jump variables (integers): $n_p")
println("  Angle variables (reals): $n_theta")
println("  Non-zeros in A: $(nnz(A))")


# Solve using greedy mixed-integer solver
println("\n" * "="^60)
println("Solving with Greedy Mixed-Integer Solver")
println("="^60)

x = greedy_mip_solver(A, b, n_p; tau=1e-6, max_iterations=1000000, verbose=false)

# Extract results
p_values = extract_integer_variables(x, n_p)
theta_values = extract_continuous_variables(x, n_p)

println("\n" * "="^60)
println("Results Summary")
println("="^60)

# Display some period jump values
println("\nPeriod jumps (first 10):")
for i in 1:min(10, length(p_values))
    edge = var_to_p[i]
    println("  p[$i] (edge $edge): $(p_values[i])")
end

# Display some angle values
println("\nFace angles (first 10):")
for i in 1:min(10, length(theta_values))
    var_idx = n_p + i
    face_idx = var_to_theta[var_idx]
    angle_deg = rad2deg(theta_values[i])
    println("  θ[face $face_idx]: $(theta_values[i]) rad ($(round(angle_deg, digits=2))°)")
end

# Compute final energy
final_energy = 0.5 * dot(x, A * x) - dot(b, x)
println("\nFinal energy: $final_energy")

# Verify constraints
println("\nVerifying constraints:")
println("  All period jumps are integers: ", all(p_values .== round.(p_values)))
if !isempty(constrained_faces)
    println("  Constrained face angles match:")
    for (i, face_idx) in enumerate(constrained_faces)
        # Check if this face is in the free system (it shouldn't be)
        var_idx = findfirst(==(face_idx), [var_to_theta[idx] for idx in keys(var_to_theta)])
        if !isnothing(var_idx)
            # Face is in free system - get from solution
            actual_var_idx = collect(keys(var_to_theta))[var_idx]
            computed_angle = x[actual_var_idx]
            expected_angle = constrained_angles[i]
            error_val = abs(computed_angle - expected_angle)
            println("    Face $face_idx: expected=$(expected_angle), computed=$(computed_angle), error=$(error_val)")
        else
            # Face is constrained (not in free system) - should exactly match
            println("    Face $face_idx: constrained to $(constrained_angles[i])")
        end
    end
else
    println("  No constrained faces specified")
end

# Visualize the computed cross field
println("\n" * "="^60)
println("Visualizing Cross Field")
println("="^60)

# Reconstruct full theta vector (including constrained faces)
theta_full = zeros(Float64, n_faces)
for (var_idx, face_idx) in var_to_theta
    theta_full[face_idx] = x[var_idx]
end
# Add constrained face angles
for (i, face_idx) in enumerate(constrained_faces)
    theta_full[face_idx] = constrained_angles[i]
end

# Create output directory
output_dir = joinpath(@__DIR__, "..", "output", "cross_field")
mkpath(output_dir)

# Get mesh name for filename
mesh_name = splitext(basename(mesh_path))[1]
savepath = joinpath(output_dir, "cross_field_$(mesh_name).png")

# Compute singularities first
println("\n" * "="^60)
println("Computing Singularities")
println("="^60)

# Build period jump dictionary from solution
period_jump_dict = Dict{Tuple{Int, Int}, Int}()
for (var_idx, edge) in var_to_p
    p_val = round(Int, x[var_idx])
    period_jump_dict[edge] = p_val
end

println("  Computing singularity indices...")
singularities = compute_singularities(mesh, period_jump_dict; debug=true)

print_singularity_report(mesh, singularities)

# Plot with singularities highlighted
include("../src/plotters.jl")
plot_cross_field_with_singularities(
    mesh, 
    theta_full,
    singularities;
    scale=0.25,
    savepath=savepath,
    field_color=:darkblue,
    field_width=1.5,
    zero_rotation_color=:red,
    zero_rotation_width=2.5,
    boundary_faces=boundary_faces,
    boundary_color=:orange,
    boundary_width=2.0,
    show_singularity_indices=true,
    singularity_marker_size=12,
    subsample=1
)

println("  Cross field plot with singularities saved to: $savepath")

println("\n=== Test Complete ===")
