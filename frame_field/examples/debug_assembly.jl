"""
Debug the matrix assembly to understand rank deficiency
"""

using Pkg
Pkg.activate(".")

using GeometryBasics
using LinearAlgebra
using SparseArrays

# Add src to path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

include(joinpath(@__DIR__, "..", "src", "meshio.jl"))
include(joinpath(@__DIR__, "..", "src", "dijkstra_forest.jl"))
include(joinpath(@__DIR__, "..", "src", "cross_field_energy.jl"))

# Load mesh
mesh_path = joinpath(@__DIR__, "..", "..", "triangulations", "simple-square.msh")
mesh = load_triangulation(mesh_path)
println("Loaded mesh with $(length(faces(mesh))) faces")

# Compute spanning forest
constrained_faces = Int[]
potential_fixed_edges = compute_spanning_forest(mesh, constrained_faces=constrained_faces)
println("Potential fixed edges: $(length(potential_fixed_edges))")

# Select fixed edges
fixed_edges_per_face = fix_suitable_edges(mesh, potential_fixed_edges)
println("Fixed edges selected: $(length(fixed_edges_per_face))")

# Assemble system with debugging
println("\n=== Assembling System ===")
A, b, var_to_p, var_to_theta = assemble_free_system(
    mesh,
    potential_fixed_edges,
    fixed_edges_per_face,
    Int[],
    Float64[],
    true
)

println("\n=== Matrix Analysis ===")
println("Matrix size: $(size(A))")
println("Non-zeros: $(nnz(A))")
println("Density: $(nnz(A) / prod(size(A)))")

# Check symmetry
sym_error = norm(A - A')
println("Symmetry error: $sym_error")

# Check diagonal
diag_A = diag(A)
println("Diagonal min: $(minimum(diag_A))")
println("Diagonal max: $(maximum(diag_A))")
println("Zero diagonal entries: $(count(==(0.0), diag_A))")

# Compute eigenvalues
println("\n=== Eigenvalue Analysis ===")
λ = eigvals(Matrix(A))
println("Eigenvalues (sorted):")
for (i, val) in enumerate(sort(real.(λ)))
    if abs(val) < 1e-10
        println("  λ[$i] = $val (≈ 0)")
    elseif i <= 5 || i >= length(λ) - 4
        println("  λ[$i] = $val")
    elseif i == 6
        println("  ...")
    end
end

# Check rank
r = rank(Matrix(A))
println("\nRank: $r / $(size(A, 1))")
println("Null space dimension: $(size(A, 1) - r)")

# Find near-zero rows
println("\n=== Row Analysis ===")
for i in 1:size(A, 1)
    row_norm = norm(A[i, :])
    if row_norm < 1e-10
        println("Row $i is nearly zero (norm = $row_norm)")
        if i <= length(var_to_p)
            println("  -> Period jump variable for edge $(var_to_p[i])")
        else
            face_idx = var_to_theta[i]
            println("  -> Angle variable for face $face_idx")
        end
    end
end

# Find near-zero columns
println("\n=== Column Analysis ===")
for j in 1:size(A, 2)
    col_norm = norm(A[:, j])
    if col_norm < 1e-10
        println("Column $j is nearly zero (norm = $col_norm)")
        if j <= length(var_to_p)
            println("  -> Period jump variable for edge $(var_to_p[j])")
        else
            face_idx = var_to_theta[j]
            println("  -> Angle variable for face $face_idx")
        end
    end
end

println("\n=== Complete ===")
