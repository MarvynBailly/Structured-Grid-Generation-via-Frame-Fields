"""
    test_cut_graph.jl

Test cut graph computation for global parameterization using actual frame field solver.
"""

include("../src/meshio.jl")
include("../src/dijkstra_forest.jl")
include("../src/cross_field_energy.jl")
include("../src/greedy_solver.jl")
include("../src/singularities.jl")
include("../src/cut_graph.jl")
include("../src/global_parameterization.jl")

using LinearAlgebra
using SparseArrays

println("=== Cut Graph with Frame Field Solver Test ===\n")

# Load mesh
# mesh_path = joinpath(@__DIR__, "..", "..", "triangulations", "simple-square.msh")
# mesh_path = joinpath(@__DIR__, "..", "..", "triangulations", "300_polygon_sphere_100mm.msh")
mesh_path = joinpath(@__DIR__, "..", "..", "triangulations", "hole.msh")

println("Loading mesh: $mesh_path")
mesh = load_triangulation(mesh_path)

fs = faces(mesh)
n_faces = length(fs)
println("  Faces: $n_faces")

# Set up additional user constraints (optional)
user_constrained_faces = Int[1]
user_constrained_angles = Float64[0.0]

# Combine boundary constraints with user constraints
constrained_faces = user_constrained_faces
constrained_angles = user_constrained_angles

println("\nTotal constrained faces: $(length(constrained_faces))")

# Compute spanning forest
println("\nComputing spanning forest...")
potential_fixed_edges = compute_spanning_forest(mesh; constrained_faces=constrained_faces)
println("  Potential fixed edges: $(length(potential_fixed_edges))")

# Select suitable fixed edges
println("\nSelecting fixed edges...")
fixed_edges_per_face = fix_suitable_edges(mesh, potential_fixed_edges)
println("  Fixed edges selected: $(length(fixed_edges_per_face))")

fixed_edges_to_use = Set(values(fixed_edges_per_face))
println("  Actually fixing edges: $(length(fixed_edges_to_use))")

# Assemble FREE variable system
println("\nAssembling FREE variable system...")
A, b, var_to_p, var_to_theta = assemble_free_system(
    mesh, 
    fixed_edges_to_use,
    fixed_edges_per_face,
    constrained_faces,
    constrained_angles,
    false  # debug mode off
)

n_vars = length(b)
n_p = length(var_to_p)
n_theta = length(var_to_theta)

println("  System size: $(size(A, 1)) × $(size(A, 2))")
println("  Total variables: $n_vars")
println("  Period jump variables: $n_p")
println("  Angle variables: $n_theta")

# Solve with greedy mixed-integer solver
println("\nSolving with greedy mixed-integer solver...")
println("  (This may take some time...)")
x = greedy_mip_solver(A, b, n_p; tau=1e-2, max_iterations=1000000, verbose=false)

# Extract results
p_values = extract_integer_variables(x, n_p)
theta_values = extract_continuous_variables(x, n_p)

println("  Solution computed successfully")

# Reconstruct full theta vector (including constrained faces)
theta_full = zeros(Float64, n_faces)
for (var_idx, face_idx) in var_to_theta
    theta_full[face_idx] = x[var_idx]
end
for (i, face_idx) in enumerate(constrained_faces)
    theta_full[face_idx] = constrained_angles[i]
end

# Build period jump dictionary from solution
period_jump_dict = Dict{Tuple{Int, Int}, Int}()
for (var_idx, edge) in var_to_p
    p_val = round(Int, x[var_idx])
    period_jump_dict[edge] = p_val
end

# Compute singularities from the solution
println("\n" * "="^60)
println("Computing Singularities")
println("="^60)

singularities = compute_singularities(mesh, period_jump_dict; debug=true)
print_singularity_report(mesh, singularities)

println("\nSingularities found: $(length(singularities))")
for (v, idx) in sort(collect(singularities), by=x->x[1])
    println("  Vertex $v: index = $idx")
end

# Compute cut graph
println("\n" * "="^60)
println("Testing Cut Graph Computation")
println("="^60)

boundary_verts = identify_boundary_vertices(mesh)
println("\nBoundary vertices: $(length(boundary_verts))")

cut_graph = compute_cut_graph(mesh, singularities; boundary_vertices=boundary_verts)

println("\nCut graph results:")
println("  Seam edges: $(length(cut_graph.seam_edges))")
println("  Singularity tree nodes: $(length(cut_graph.singularity_tree))")

# Show some seam edges
if !isempty(cut_graph.seam_edges)
    println("\nFirst 10 seam edges:")
    for (i, edge) in enumerate(collect(cut_graph.seam_edges)[1:min(10, length(cut_graph.seam_edges))])
        println("  Edge $i: $edge")
    end
end

# Visualize the cut graph
println("\n" * "="^60)
println("Visualizing Cut Graph")
println("="^60)

include("../src/plotters.jl")

output_dir = joinpath(@__DIR__, "..", "output", "cut_graph")
mkpath(output_dir)

mesh_name = splitext(basename(mesh_path))[1]
savepath = joinpath(output_dir, "cut_graph_$(mesh_name).png")

plot_cut_graph(
    mesh,
    cut_graph,
    singularities;
    savepath=savepath,
    cut_color=:red,
    cut_width=4.0,
    singularity_size=20,
    show_vertex_labels=false
)

println("Cut graph visualization saved to: $savepath")

# Test global parameterization (placeholder)
println("\n" * "="^60)
println("Testing Global Parameterization (Placeholder)")
println("="^60)

param = compute_global_parameterization(
    mesh, theta_full, period_jump_dict, singularities;
    method=:poisson
)

println("\nParameterization results:")
println("  u range: [$(minimum(param.u_coords)), $(maximum(param.u_coords))]")
println("  v range: [$(minimum(param.v_coords)), $(maximum(param.v_coords))]")
println("  Seam matching pairs: $(length(param.seam_matching))")

println("\n=== Test Complete ===")
