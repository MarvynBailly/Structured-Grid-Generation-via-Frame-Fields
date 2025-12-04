"""
    energy.jl

Energy function computation and derivatives for frame field optimization.

Implements the smoothness energy:
E_smooth = Σ (θ_i + κ_ij + (π/2)p_ij - θ_j)²
"""

using LinearAlgebra
using SparseArrays
using GeometryBasics

"""
    compute_smoothness_energy(
        mesh::Mesh,
        topology::MeshTopology,
        theta::Vector{Float64},
        p::Dict{Int, Int}
    ) -> Float64

Compute the frame field smoothness energy.

# Arguments
- `mesh::Mesh` - Triangular mesh
- `topology::MeshTopology` - Mesh topology
- `theta::Vector{Float64}` - Face angles (length = n_faces)
- `p::Dict{Int, Int}` - Period jumps indexed by edge

# Returns
- `Float64` - Smoothness energy value
"""
function compute_smoothness_energy(
    mesh::Mesh,
    topology::MeshTopology,
    theta::Vector{Float64},
    p::Dict{Int, Int}
)
    energy = 0.0
    fs = faces(mesh)
    n_faces = length(fs)
    
    # Sum over all edges
    for face_i in 1:n_faces
        for (face_j, edge_idx) in topology.dual_adj[face_i]
            if face_i < face_j  # Count each edge once
                kappa_ij = topology.kappa[(face_i, face_j)]
                p_ij = get(p, edge_idx, 0)
                
                diff = theta[face_i] + kappa_ij + (π/2) * p_ij - theta[face_j]
                energy += diff^2
            end
        end
    end
    
    return energy
end

"""
    assemble_system_matrix(
        mesh::Mesh,
        topology::MeshTopology,
        forest::DijkstraForest,
        constrained_faces::Vector{Int},
        constrained_angles::Vector{Float64}
    ) -> (SparseMatrixCSC, Vector{Float64}, Dict, Dict)

Assemble the linear system for frame field optimization.

# Arguments
- `mesh::Mesh` - Triangular mesh
- `topology::MeshTopology` - Mesh topology
- `forest::DijkstraForest` - Dijkstra forest (identifies fixed period jumps)
- `constrained_faces::Vector{Int}` - Faces with fixed angles
- `constrained_angles::Vector{Float64}` - Corresponding angle values

# Returns
- `A::SparseMatrixCSC` - System matrix
- `b::Vector{Float64}` - Right-hand side vector
- `var_to_theta::Dict{Int, Int}` - Maps variable index to face index for θ variables
- `var_to_p::Dict{Int, Int}` - Maps variable index to edge index for p variables

# Algorithm
Takes derivatives of E_smooth with respect to each free variable:
- ∂E/∂θ_k = Σ 2(θ_k + κ_kj + (π/2)p_kj - θ_j) = 0
- ∂E/∂p_ij = π(θ_i + κ_ij + (π/2)p_ij - θ_j) = 0

Assembles these into matrix form Ax = b where:
- First k variables are period jumps (integers)
- Remaining variables are face angles (reals)
"""
function assemble_system_matrix(
    mesh::Mesh,
    topology::MeshTopology,
    forest::DijkstraForest,
    constrained_faces::Vector{Int}=Int[],
    constrained_angles::Vector{Float64}=Float64[]
)
    fs = faces(mesh)
    n_faces = length(fs)
    
    # Determine free variables
    free_faces = setdiff(1:n_faces, constrained_faces)
    free_edges = forest.free_edges
    
    n_free_theta = length(free_faces)
    n_free_p = length(free_edges)
    n_vars = n_free_p + n_free_theta
    
    # Create variable mappings
    var_to_p = Dict{Int, Int}()
    p_to_var = Dict{Int, Int}()
    for (i, edge_idx) in enumerate(free_edges)
        var_to_p[i] = edge_idx
        p_to_var[edge_idx] = i
    end
    
    var_to_theta = Dict{Int, Int}()
    theta_to_var = Dict{Int, Int}()
    for (i, face_idx) in enumerate(free_faces)
        var_idx = n_free_p + i
        var_to_theta[var_idx] = face_idx
        theta_to_var[face_idx] = var_idx
    end
    
    # Build sparse system
    I_idx = Int[]
    J_idx = Int[]
    vals = Float64[]
    b = zeros(n_vars)
    
    row = 0
    
    # Equations from ∂E/∂p_ij = 0 for each free edge
    for (var_idx, edge_idx) in sort(collect(var_to_p))
        # Find faces incident to this edge
        incident_faces = topology.edge_to_faces[edge_idx]
        if length(incident_faces) != 2
            # Boundary edge - add trivial equation p_ij = 0
            row += 1
            push!(I_idx, row)
            push!(J_idx, var_idx)
            push!(vals, 1.0)
            b[row] = 0.0
            continue
        end
        
        row += 1
        
        face_i, face_j = incident_faces
        kappa_ij = topology.kappa[(face_i, face_j)]
        
        # Coefficient for p_ij itself: π² / 2
        push!(I_idx, row)
        push!(J_idx, var_idx)
        push!(vals, π^2 / 2)
        
        # Coefficient for θ_i: π (if free)
        if face_i in free_faces
            theta_var = theta_to_var[face_i]
            push!(I_idx, row)
            push!(J_idx, theta_var)
            push!(vals, π)
        else
            # Constrained, move to RHS
            idx = findfirst(==(face_i), constrained_faces)
            if !isnothing(idx)
                b[row] -= π * constrained_angles[idx]
            end
        end
        
        # Coefficient for θ_j: -π (if free)
        if face_j in free_faces
            theta_var = theta_to_var[face_j]
            push!(I_idx, row)
            push!(J_idx, theta_var)
            push!(vals, -π)
        else
            # Constrained, move to RHS
            idx = findfirst(==(face_j), constrained_faces)
            if !isnothing(idx)
                b[row] += π * constrained_angles[idx]
            end
        end
        
        # Constant term from κ_ij
        b[row] -= π * kappa_ij
    end
    
    # Equations from ∂E/∂θ_k = 0 for each free face
    for face_k in free_faces
        row += 1
        theta_var = theta_to_var[face_k]
        
        # Sum over all neighbors of face_k
        degree = length(topology.dual_adj[face_k])
        
        # Diagonal term: 2 * degree
        push!(I_idx, row)
        push!(J_idx, theta_var)
        push!(vals, 2.0 * degree)
        
        for (face_j, edge_idx) in topology.dual_adj[face_k]
            kappa_kj = topology.kappa[(face_k, face_j)]
            
            # Coefficient for p_kj: π (if free)
            if edge_idx in free_edges
                p_var = p_to_var[edge_idx]
                push!(I_idx, row)
                push!(J_idx, p_var)
                push!(vals, π)
            else
                # Fixed to zero, no contribution
            end
            
            # Coefficient for θ_j: -2 (if free)
            if face_j in free_faces
                theta_var_j = theta_to_var[face_j]
                push!(I_idx, row)
                push!(J_idx, theta_var_j)
                push!(vals, -2.0)
            else
                # Constrained, move to RHS
                idx = findfirst(==(face_j), constrained_faces)
                if !isnothing(idx)
                    b[row] += 2.0 * constrained_angles[idx]
                end
            end
            
            # Constant term from κ_kj
            b[row] -= 2.0 * kappa_kj
        end
    end
    
    # Build sparse matrix
    A = sparse(I_idx, J_idx, vals, n_vars, n_vars)
    
    return A, b, var_to_p, var_to_theta
end

"""
    compute_singularity_index(
        mesh::Mesh,
        topology::MeshTopology,
        p::Dict{Int, Int},
        vertex_idx::Int
    ) -> Float64

Compute the cross field index at a vertex.

I(v) = I₀(v) + Σ(p_ij / 4)

where the sum is over edges in star(v).
"""
function compute_singularity_index(
    mesh::Mesh,
    topology::MeshTopology,
    p::Dict{Int, Int},
    vertex_idx::Int
)
    vs = coordinates(mesh)
    fs = faces(mesh)
    
    # Compute geometric index I₀(v)
    angle_defect = compute_angle_defect(mesh, vertex_idx)
    
    # Sum κ_ij over edges in star(v)
    kappa_sum = 0.0
    p_sum = 0.0
    
    for face_idx in 1:length(fs)
        face = fs[face_idx]
        if vertex_idx in face
            # Check edges incident to this vertex
            for (neighbor_idx, edge_idx) in topology.dual_adj[face_idx]
                kappa_ij = topology.kappa[(face_idx, neighbor_idx)]
                kappa_sum += kappa_ij
                
                p_ij = get(p, edge_idx, 0)
                p_sum += p_ij / 4.0
            end
        end
    end
    
    I_0 = (angle_defect + kappa_sum) / (2π)
    I = I_0 + p_sum
    
    return I
end
