"""
Verify that the assembled matrix has the correct block structure as shown in the slides.
"""

include("../src/meshio.jl")
include("../src/mesh_types.jl")
include("../src/dijkstra_forest.jl")
include("../src/cross_field_energy.jl")

using LinearAlgebra
using SparseArrays

println("=== Matrix Block Structure Verification ===\n")

# Load simple mesh
mesh_path = joinpath(@__DIR__, "..", "..", "triangulations", "simplest-square.msh")
println("Loading mesh: $mesh_path")
mesh = load_triangulation(mesh_path)

fs = faces(mesh)
n_faces = length(fs)
println("  Faces: $n_faces\n")

# Build topology and create frame field
println("Building mesh topology...")
topology = build_mesh_topology(mesh)

println("Creating frame field...")
field = CrossField(topology)

# Set up simple constraints
constrained_faces = Int[1]
constrained_angles = Float64[0.0]

# Compute spanning forest
println("Computing spanning forest...")
potential_fixed_edges = compute_spanning_forest(mesh; constrained_faces=constrained_faces)

# Select suitable fixed edges
println("Selecting fixed edges...")
fixed_edges_per_face = fix_suitable_edges(mesh, potential_fixed_edges)

# Assemble system matrix
println("\nAssembling system matrix...")
A, b, var_to_p, var_to_theta = assemble_system_matrix(
    field, 
    fixed_edges_per_face,
    constrained_faces,
    constrained_angles,
    false
)

n_p = length(var_to_p)
n_theta = length(var_to_theta)
n_vars = size(A, 1)

println("  System size: $(size(A))")
println("  Number of p variables: $n_p")
println("  Number of θ variables: $n_theta")
println("  Total variables: $n_vars")
println()

# Extract blocks
A_full = Matrix(A)

# Top-left: p × p block
A_pp = A_full[1:n_p, 1:n_p]
println("Top-left block (p × p): $(size(A_pp))")
println("  Should be diagonal with π²/2 on diagonal for free edges, 1.0 for fixed edges")
println("  Is diagonal: $(isdiag(A_pp))")
if isdiag(A_pp)
    diag_vals = unique(diag(A_pp))
    println("  Diagonal values: $diag_vals")
    println("  Expected: [1.0, $(π^2/2)]")
end
println()

# Top-right: p × θ block (should be π times incidence matrix transpose)
A_pθ = A_full[1:n_p, (n_p+1):n_vars]
println("Top-right block (p × θ): $(size(A_pθ))")
println("  Should be ±π for edge-face incidences")
nonzero_vals = unique(filter(x -> abs(x) > 1e-10, A_pθ[:]))
println("  Non-zero values: $nonzero_vals")
println("  Expected: [±π] = [±$(π)]")
println()

# Bottom-left: θ × p block (should be π times incidence matrix)
A_θp = A_full[(n_p+1):n_vars, 1:n_p]
println("Bottom-left block (θ × p): $(size(A_θp))")
println("  Should be ±π for face-edge incidences")
nonzero_vals = unique(filter(x -> abs(x) > 1e-10, A_θp[:]))
println("  Non-zero values: $nonzero_vals")
println("  Expected: [±π] = [±$(π)]")
println()

# Bottom-right: θ × θ block (should be 2 times graph Laplacian for free, 1 for fixed)
A_θθ = A_full[(n_p+1):n_vars, (n_p+1):n_vars]
println("Bottom-right block (θ × θ): $(size(A_θθ))")
println("  Should have 2×degree on diagonal, -2 off-diagonal for free faces, 1 for fixed")
diag_vals = unique(diag(A_θθ))
println("  Diagonal values: $diag_vals")
offdiag_vals = unique(filter(x -> abs(x) > 1e-10, [A_θθ[i,j] for i in 1:size(A_θθ,1) for j in 1:size(A_θθ,2) if i != j]))
println("  Off-diagonal values: $offdiag_vals")
println("  Expected diagonal: [1.0, 2×degree] for degree ≥ 1")
println("  Expected off-diagonal: [-2.0]")
println()

# Check if bottom-left = (π/2) × transpose(top-right) × 2
println("=== Block Structure Relationships ===")
println()

# For the continuous relaxation, we should have:
# A_θp should be related to transpose of A_pθ
# Specifically: A_θp = (2/π) × A_pθ^T for the energy derivative structure

# Check relationship (only for free variables)
println("Checking block relationship:")

# The relationship should be: bottom-left has π, top-right has π
# So they should be related by transpose (with sign conventions)
println("\nComparing A_pθ and A_θp^T for FREE variables only:")

# Identify free and fixed p indices
free_p_indices = [i for i in 1:n_p if A_pp[i,i] ≈ π^2/2]
fixed_p_indices = [i for i in 1:n_p if A_pp[i,i] ≈ 1.0]

# Identify free and fixed θ indices  
free_theta_indices = [i for i in 1:n_theta if A_θθ[i,i] > 1.5]  # Free faces have degree ≥ 2
fixed_theta_indices = [i for i in 1:n_theta if A_θθ[i,i] ≈ 1.0]

println("  Free p variables: $(length(free_p_indices))")
println("  Fixed p variables: $(length(fixed_p_indices))")
println("  Free θ variables: $(length(free_theta_indices))")
println("  Fixed θ variables: $(length(fixed_theta_indices))")

# Extract submatrices for free variables only
A_pθ_free = A_pθ[free_p_indices, free_theta_indices]
A_θp_free = A_θp[free_theta_indices, free_p_indices]

diff_free = sum(abs.(A_pθ_free - A_θp_free'))
println("\n  Transpose check for FREE variables:")
println("    Sum of |A_pθ_free - A_θp_free^T|: $diff_free")

if diff_free < 1e-10
    println("    ✓✓✓ FREE variable blocks are exact transposes!")
else
    println("    Maximum difference: $(maximum(abs.(A_pθ_free - A_θp_free')))")
end

# Check that fixed p rows in A_pθ are zero
fixed_p_nonzero = sum(abs.(A_pθ[fixed_p_indices, :]))
println("\n  Fixed p variable rows in A_pθ:")
println("    Sum of absolute values: $fixed_p_nonzero")
if fixed_p_nonzero < 1e-10
    println("    ✓ Fixed p rows are zero (correct constraint equations)")
else
    println("    ✗ Fixed p rows should be zero!")
end

println("\n=== Verification Complete ===")
