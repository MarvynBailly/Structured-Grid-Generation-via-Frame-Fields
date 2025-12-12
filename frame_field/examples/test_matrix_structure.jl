"""
    test_matrix_structure.jl

Visualize the matrix structure of the Ax=b minimization problem
on the simple 3-triangle mesh.

The optimization problem (similar to cross_field_energy.jl):
    min E = Σ (θ_i + κ_ij + (π/2)p_ij - θ_j)²

Taking derivatives with respect to p_ij and θ_i:
    ∂E/∂p_ij = π(θ_i + κ_ij + (π/2)p_ij - θ_j) = 0
    ∂E/∂θ_i = Σ_j 2(θ_i + κ_ij + (π/2)p_ij - θ_j) = 0

Variables ordered as: [p_1, p_2, ..., θ_1, θ_2, ...]
- Period jumps p_ij for dual edges where i < j  
- Angles θ_i for all faces
"""

using Pkg
Pkg.activate(".")

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using GeometryBasics
using CairoMakie
using LinearAlgebra
using SparseArrays

include("../src/meshio.jl")
include("../src/mesh_types.jl")

folder_name = "matrix_structure"
# case_name = "three-triangles"
case_name = "simple-square"

function main()
    println("Loading $case_name mesh...")
    mesh_file = joinpath(@__DIR__, "..", "..", "triangulations", "$case_name.msh")
    mesh = load_triangulation(mesh_file)
    
    println("Building mesh topology...")
    topology = build_mesh_topology(mesh)
    
    println("Creating frame field...")
    field = CrossField(topology)
    
    n_faces = topology.n_faces
    
    println("\n" * "="^60)
    println("PROBLEM SETUP")
    println("="^60)
    
    # Collect all dual edges
    dual_edges = Tuple{Int, Int}[]
    for face_i in 1:n_faces
        for (face_j, edge) in topology.dual_adj[face_i]
            if face_i < face_j
                push!(dual_edges, (face_i, face_j))
            end
        end
    end
    
    n_p_vars = length(dual_edges)
    n_theta_vars = n_faces
    n_vars = n_p_vars + n_theta_vars
    n_equations = length(dual_edges) + n_faces  # One equation per p_ij, one per θ_i
    
    println("\nVariables:")
    println("  Period jumps p_ij: $n_p_vars")
    for (i, (face_i, face_j)) in enumerate(dual_edges)
        println("    x_$i = p_($face_i,$face_j)")
    end
    println("  Angles θ_i: $n_theta_vars")
    for i in 1:n_faces
        println("    x_$(n_p_vars + i) = θ_$i")
    end
    println("  Total variables: $n_vars")
    
    println("\nEquations (from ∂E/∂p and ∂E/∂θ = 0):")
    println("  ∂E/∂p_ij equations: $n_p_vars")
    for (face_i, face_j) in dual_edges
        kappa_ij = get_transport_angle(field, face_i, face_j)
        println("    (π/2)p_($face_i,$face_j) + θ_$face_i - θ_$face_j = -κ_($face_i,$face_j)")
        println("      where κ_($face_i,$face_j) = $(round(rad2deg(kappa_ij), digits=2))°")
    end
    println("  ∂E/∂θ_i equations: $n_theta_vars")
    for i in 1:n_faces
        neighbors = [n for (n, _) in topology.dual_adj[i]]
        println("    Σ_j (θ_$i + κ_($i,j) + (π/2)p_($i,j) - θ_j) = 0 for j ∈ $neighbors")
    end
    
    println("  Kappa_ij values:")
    for (face_i, face_j) in dual_edges
        kappa_ij = get_transport_angle(field, face_i, face_j)
        println("    κ_($face_i,$face_j) = $(round(rad2deg(kappa_ij), digits=2))°")
    end

    # Build the system matrix A and vector b
    println("\n" * "="^60)
    println("BUILDING SYSTEM MATRIX Ax = b")
    println("="^60)
    
    # Create sparse matrix A
    A = spzeros(Float64, n_equations, n_vars)
    b = zeros(Float64, n_equations)
    
    # Map (face_i, face_j) to p variable index
    p_var_idx = Dict{Tuple{Int,Int}, Int}()
    for (idx, (face_i, face_j)) in enumerate(dual_edges)
        p_var_idx[(face_i, face_j)] = idx
    end
    
    row = 0
    
    # First set of equations: ∂E/∂p_ij = π(θ_i + κ_ij + (π/2)p_ij - θ_j) = 0
    # Simplifies to: (π/2)p_ij + θ_i - θ_j = -κ_ij
    for (face_i, face_j) in dual_edges
        row += 1
        
        p_idx = p_var_idx[(face_i, face_j)]
        theta_i_idx = n_p_vars + face_i
        theta_j_idx = n_p_vars + face_j
        
        A[row, p_idx] = π/2
        A[row, theta_i_idx] = 1.0
        A[row, theta_j_idx] = -1.0
        
        kappa_ij = get_transport_angle(field, face_i, face_j)
        b[row] = -kappa_ij
    end
    
    # Second set of equations: ∂E/∂θ_i = Σ_j 2(θ_i + κ_ij + (π/2)p_ij - θ_j) = 0
    # For each face i, sum over all neighbors j
    for face_i in 1:n_faces
        row += 1
        theta_i_idx = n_p_vars + face_i
        
        # Degree of face i
        degree = length(topology.dual_adj[face_i])
        A[row, theta_i_idx] = degree
        
        rhs = 0.0
        for (face_j, edge) in topology.dual_adj[face_i]
            # Determine p_ij index and sign
            if face_i < face_j
                p_idx = p_var_idx[(face_i, face_j)]
                p_sign = 1.0
                kappa_ij = get_transport_angle(field, face_i, face_j)
            else
                p_idx = p_var_idx[(face_j, face_i)]
                p_sign = -1.0
                kappa_ij = get_transport_angle(field, face_i, face_j)
            end
            
            A[row, p_idx] += (π/2) * p_sign
            
            theta_j_idx = n_p_vars + face_j
            A[row, theta_j_idx] -= 1.0
            
            rhs -= 1.0 * kappa_ij
        end
        b[row] = rhs
    end
    
    println("\nMatrix A dimensions: $(size(A))")
    println("Vector b dimensions: $(size(b))")
    println("\nMatrix A (full):")
    display(Matrix(A))
    println("\n\nVector b ($(length(b))):")
    display(rad2deg.(b)) 
    
    # Create visualization
    println("\n" * "="^60)
    println("CREATING VISUALIZATIONS")
    println("="^60)
    
    fig = Figure(size=(1600, 1200))
    
    # Matrix A structure
    ax1 = Axis(fig[1, 1], aspect=DataAspect(), 
               title="Matrix A Structure ($(n_equations)×$(n_vars))",
               xlabel="Variables", ylabel="Equations")
    
    # Plot non-zero entries
    rows, cols, vals = findnz(A)
    colors = [v > 0 ? :blue : :red for v in vals]
    scatter!(ax1, cols, rows, color=colors, markersize=20)
    
    # Add grid
    hlines!(ax1, 1:n_equations, color=:gray, linewidth=0.5, alpha=0.3)
    vlines!(ax1, 1:n_vars, color=:gray, linewidth=0.5, alpha=0.3)
    
    # Add variable labels
    var_labels = ["p_$i" for i in 1:n_p_vars] ∪ ["θ_$i" for i in 1:n_theta_vars]
    ax1.xticks = (1:n_vars, var_labels)
    ax1.xticklabelrotation = π/4
    
    # Invert y-axis so equation 1 is at top
    ylims!(ax1, n_equations + 0.5, 0.5)
    
    # Matrix A with values
    ax2 = Axis(fig[2, 1], aspect=DataAspect(),
               title="Matrix A with Values",
               xlabel="Variables", ylabel="Equations")
    
    heatmap!(ax2, 1:n_vars, 1:n_equations, Matrix(A), colormap=:RdBu)
    
    # Add text annotations for values
    for i in 1:n_equations
        for j in 1:n_vars
            val = A[i, j]
            if abs(val) > 1e-5
                val_str = abs(val - π/2) < 1e-5 ? "π/2" : 
                         abs(val - 1.0) < 1e-5 ? "1" :
                         abs(val + 1.0) < 1e-5 ? "-1" : 
                         abs(val - 2.0) < 1e-5 ? "2" :
                         abs(val + 2.0) < 1e-5 ? "-2" :
                         abs(val - π) < 1e-5 ? "π" :
                         abs(val + π) < 1e-5 ? "-π" :
                         "$(round(val, digits=2))"
                text!(ax2, j, i, text=val_str, color=:white, 
                      fontsize=10, align=(:center, :center))
            end
        end
    end
    
    ax2.xticks = (1:n_vars, var_labels)
    ax2.xticklabelrotation = π/4
    ylims!(ax2, n_equations + 0.5, 0.5)
    
    # Save figure
    output_dir = joinpath(@__DIR__, "..", "output", folder_name)
    mkpath(output_dir)
    output_file = joinpath(output_dir, "$(case_name).png")
    save(output_file, fig)
    println("\nVisualization saved to: $output_file")
    
    println("\n" * "="^60)
end

main()
