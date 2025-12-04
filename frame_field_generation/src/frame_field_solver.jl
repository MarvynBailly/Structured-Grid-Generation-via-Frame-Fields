"""
    frame_field_solver.jl

Main solver for frame field generation using greedy mixed-integer optimization.
"""

using LinearAlgebra
using SparseArrays
using GeometryBasics

# Import GreedyMIPSolver from parent directory
include("../../greedy_mip_solver/src/GreedyMIPSolver.jl")
using .GreedyMIPSolver

"""
    FrameFieldResult

Container for frame field generation results.

# Fields
- `angles::Vector{Float64}` - Face angles θ_f for each face
- `period_jumps::Dict{Int, Int}` - Period jumps p_e indexed by edge
- `singularities::Vector{Tuple{Int, Float64}}` - (vertex_idx, index) for singular vertices
- `energy::Float64` - Final smoothness energy value
- `converged::Bool` - Whether solver converged
- `iterations::Int` - Number of iterations
"""
struct FrameFieldResult
    angles::Vector{Float64}
    period_jumps::Dict{Int, Int}
    singularities::Vector{Tuple{Int, Float64}}
    energy::Float64
    converged::Bool
    iterations::Int
end

"""
    generate_frame_field(
        mesh::Mesh;
        constrained_faces::Vector{Int}=Int[],
        constrained_angles::Vector{Float64}=Float64[],
        τ::Float64=1e-6,
        max_iterations::Int=10000,
        verbose::Bool=false
    ) -> FrameFieldResult

Generate a smooth frame field over a triangular mesh.

# Arguments
- `mesh::Mesh` - Input triangular mesh
- `constrained_faces::Vector{Int}` - Face indices with fixed orientation
- `constrained_angles::Vector{Float64}` - Corresponding angle values
- `τ::Float64` - Convergence tolerance for MIP solver
- `max_iterations::Int` - Maximum iterations for MIP solver
- `verbose::Bool` - Print progress information

# Returns
- `FrameFieldResult` - Complete frame field solution

# Algorithm
1. Build mesh topology (edges, dual graph, reference edges)
2. Construct Dijkstra forest to identify fixed period jumps
3. Assemble linear system from energy derivatives
4. Solve mixed-integer problem using greedy solver
5. Extract θ and p values from solution
6. Compute singularities
7. Calculate final energy

# Example
```julia
mesh = generate_square_mesh(nx=10, ny=10)
result = generate_frame_field(mesh, verbose=true)
println("Energy: ", result.energy)
println("Singularities: ", length(result.singularities))
```
"""
function generate_frame_field(
    mesh::Mesh;
    constrained_faces::Vector{Int}=Int[],
    constrained_angles::Vector{Float64}=Float64[],
    τ::Float64=1e-6,
    max_iterations::Int=10000,
    verbose::Bool=false
)
    verbose && println("=== Frame Field Generation ===")
    
    # Step 1: Build mesh topology
    verbose && println("\n1. Building mesh topology...")
    topology = build_mesh_topology(mesh)
    verbose && println("   - Faces: $(length(faces(mesh)))")
    verbose && println("   - Edges: $(topology.n_edges)")
    
    # Step 2: Build Dijkstra forest
    verbose && println("\n2. Building Dijkstra forest...")
    forest = build_dijkstra_forest(mesh, topology, constrained_faces)
    n_free_p = count_free_edges(forest)
    n_fixed = length(forest.fixed_edges_per_face)
    verbose && println("   - Fixable edges (in tree): $(length(forest.fixable_edges))")
    verbose && println("   - Fixed edges (one per triangle): $n_fixed")
    verbose && println("   - Free edges (p variables): $n_free_p")
    
    # Step 3: Assemble system matrix
    verbose && println("\n3. Assembling linear system...")
    A, b, var_to_p, var_to_theta = assemble_system_matrix(
        mesh, topology, forest, constrained_faces, constrained_angles
    )
    n_vars = length(b)
    n_free_theta = length(var_to_theta)
    verbose && println("   - Total variables: $n_vars")
    verbose && println("   - Integer variables (p): $n_free_p")
    verbose && println("   - Real variables (θ): $n_free_theta")
    
    # Step 4: Solve mixed-integer problem
    verbose && println("\n4. Solving mixed-integer problem...")
    mip_result = solve_greedy_mip(A, b, n_free_p, τ=τ, max_iterations=max_iterations, verbose=verbose)
    
    if !mip_result.converged
        @warn "MIP solver did not fully converge"
    end
    
    verbose && println("   - Converged: $(mip_result.converged)")
    verbose && println("   - Iterations: $(mip_result.iterations)")
    
    # Step 5: Extract solution
    verbose && println("\n5. Extracting solution...")
    x_sol = mip_result.x
    
    # Initialize arrays
    n_faces = length(faces(mesh))
    angles = zeros(n_faces)
    period_jumps = Dict{Int, Int}()
    
    # Extract period jumps (first n_free_p variables)
    for (var_idx, edge_idx) in var_to_p
        period_jumps[edge_idx] = round(Int, x_sol[var_idx])
    end
    
    # Fixed edges (one per face) have p = 0
    for (face_idx, edge_idx) in forest.fixed_edges_per_face
        period_jumps[edge_idx] = 0
    end
    
    # Extract angles (remaining variables)
    for (var_idx, face_idx) in var_to_theta
        angles[face_idx] = x_sol[var_idx]
    end
    
    # Constrained faces keep their specified angles
    for (i, face_idx) in enumerate(constrained_faces)
        angles[face_idx] = constrained_angles[i]
    end
    
    # Step 6: Compute singularities
    verbose && println("\n6. Computing singularities...")
    singularities = compute_singularities(mesh, topology, period_jumps)
    verbose && println("   - Singular vertices: $(length(singularities))")
    
    # Step 7: Compute final energy
    energy = compute_smoothness_energy(mesh, topology, angles, period_jumps)
    verbose && println("\n7. Final smoothness energy: $energy")
    
    verbose && println("\n=== Frame Field Generation Complete ===\n")
    
    return FrameFieldResult(
        angles,
        period_jumps,
        singularities,
        energy,
        mip_result.converged,
        mip_result.iterations
    )
end

"""
    compute_singularities(
        mesh::Mesh,
        topology::MeshTopology,
        p::Dict{Int, Int};
        tol::Float64=1e-6
    ) -> Vector{Tuple{Int, Float64}}

Find all singular vertices in the frame field.

# Arguments
- `mesh::Mesh` - Triangular mesh
- `topology::MeshTopology` - Mesh topology
- `p::Dict{Int, Int}` - Period jumps
- `tol::Float64` - Tolerance for considering index as non-zero

# Returns
- `Vector{Tuple{Int, Float64}}` - List of (vertex_index, singularity_index) pairs

A vertex is singular if its cross field index is non-zero (within tolerance).
"""
function compute_singularities(
    mesh::Mesh,
    topology::MeshTopology,
    p::Dict{Int, Int};
    tol::Float64=1e-6
)
    vs = coordinates(mesh)
    n_vertices = length(vs)
    
    singularities = Tuple{Int, Float64}[]
    
    for v_idx in 1:n_vertices
        index = compute_singularity_index(mesh, topology, p, v_idx)
        
        if abs(index) > tol
            push!(singularities, (v_idx, index))
        end
    end
    
    return singularities
end

"""
    get_frame_directions(result::FrameFieldResult, face_idx::Int) -> Vector{Float64}

Get the four frame directions at a given face.

Returns angles [θ_f, θ_f + π/2, θ_f + π, θ_f + 3π/2]
"""
function get_frame_directions(result::FrameFieldResult, face_idx::Int)
    θ = result.angles[face_idx]
    return [θ, θ + π/2, θ + π, θ + 3π/2]
end

"""
    get_frame_directions(mesh::Mesh, result::FrameFieldResult) -> Vector{Vector{Vector{Float64}}}

Get frame direction vectors for all faces in the mesh.

Returns a vector where each element is a vector of 4 direction vectors (2D unit vectors)
corresponding to the 4 frame directions at that face.
"""
function get_frame_directions(mesh::Mesh, result::FrameFieldResult)
    fs = faces(mesh)
    n_faces = length(fs)
    
    directions = Vector{Vector{Vector{Float64}}}(undef, n_faces)
    
    for face_idx in 1:n_faces
        θ = result.angles[face_idx]
        
        # Compute 4 unit direction vectors (90° apart)
        dirs = Vector{Vector{Float64}}(undef, 4)
        for i in 0:3
            angle = θ + i * π/2
            dirs[i+1] = [cos(angle), sin(angle)]
        end
        
        directions[face_idx] = dirs
    end
    
    return directions
end
