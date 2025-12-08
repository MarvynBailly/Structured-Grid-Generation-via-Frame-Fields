"""
    global_parameterization.jl

Global seamless parameterization for quad mesh generation.

After computing the cross field, we need to integrate it to get a global
parameterization (u,v) coordinates. The key challenges:

1. **Cutting**: Closed surfaces must be cut to create a topological disk
2. **Seamless constraints**: Grid lines must match across cut boundaries
3. **Integer constraints**: Jumps across seams must be grid symmetries

The parameterization satisfies: ∇u × ∇v aligns with the cross field direction.
"""

using LinearAlgebra
using SparseArrays
using GeometryBasics

include("cut_graph.jl")


"""
    GlobalParameterization

Stores the global (u,v) parameterization of the mesh.

# Fields
- `u_coords::Vector{Float64}` - u-coordinate at each vertex
- `v_coords::Vector{Float64}` - v-coordinate at each vertex  
- `seam_matching::Dict{Int,Int}` - Maps vertices on one side of seam to other side
- `cut_graph::CutGraph` - The cut graph used
"""
struct GlobalParameterization
    u_coords::Vector{Float64}
    v_coords::Vector{Float64}
    seam_matching::Dict{Int,Int}
    cut_graph::CutGraph
end


"""
    compute_global_parameterization(
        mesh::Mesh,
        theta::Vector{Float64},
        period_jumps::Dict{Tuple{Int,Int},Int},
        singularities::Dict{Int,Float64};
        method=:poisson
    ) -> GlobalParameterization

Compute global seamless parameterization by integrating the cross field.

# Algorithm
1. Compute cut graph connecting singularities
2. Duplicate vertices along seam (topologically cut the mesh)
3. Set up Poisson system: Δu = ∇×f_u, Δv = ∇×f_v
4. Add seam compatibility constraints (integer translations + rotations)
5. Solve for (u,v) coordinates at all vertices

# Arguments
- `mesh`: Triangular mesh
- `theta`: Face angles from frame field optimization
- `period_jumps`: Period jump integers on edges
- `singularities`: Detected singularities
- `method`: Integration method (:poisson or :greedy)

# Returns
- `GlobalParameterization` with (u,v) coordinates and seam information
"""
function compute_global_parameterization(
    mesh::GeometryBasics.Mesh,
    theta::Vector{Float64},
    period_jumps::Dict{Tuple{Int,Int},Int},
    singularities::Dict{Int,Float64};
    method::Symbol=:poisson
)
    println("\n" * "="^60)
    println("Global Parameterization")
    println("="^60)
    
    # Step 1: Compute cut graph
    println("\nStep 1: Computing cut graph...")
    boundary_verts = identify_boundary_vertices(mesh)
    has_boundary = !isempty(boundary_verts)
    
    if has_boundary
        println("  Mesh has boundary with $(length(boundary_verts)) vertices")
    else
        println("  Closed surface - will compute cut graph")
    end
    
    cut_graph = compute_cut_graph(mesh, singularities; 
                                   boundary_vertices=boundary_verts)
    
    # Step 2: Set up parameterization system
    println("\nStep 2: Setting up parameterization system...")
    
    if method == :poisson
        u_coords, v_coords, seam_matching = solve_poisson_parameterization(
            mesh, theta, period_jumps, cut_graph
        )
    else
        error("Method $method not implemented yet")
    end
    
    println("\nGlobal parameterization computed successfully")
    println("  u range: [$(minimum(u_coords)), $(maximum(u_coords))]")
    println("  v range: [$(minimum(v_coords)), $(maximum(v_coords))]")
    
    return GlobalParameterization(u_coords, v_coords, seam_matching, cut_graph)
end


"""
    solve_poisson_parameterization(mesh, theta, period_jumps, cut_graph)
    -> (u_coords, v_coords, seam_matching)

Solve for global parameterization using Poisson integration.

The cross field defines a direction field. We integrate this to get (u,v):
- ∂u/∂x, ∂u/∂y aligned with cross direction
- ∂v/∂x, ∂v/∂y perpendicular to cross direction

This becomes a Poisson equation: Δu = curl(field), with seam constraints.
"""
function solve_poisson_parameterization(
    mesh::GeometryBasics.Mesh,
    theta::Vector{Float64},
    period_jumps::Dict{Tuple{Int,Int},Int},
    cut_graph::CutGraph
)
    vs = coordinates(mesh)
    fs = faces(mesh)
    n_verts = length(vs)
    n_faces = length(fs)
    
    println("  Building parameterization system...")
    println("    Vertices: $n_verts")
    println("    Faces: $n_faces")
    println("    Seam edges: $(length(cut_graph.seam_edges))")
    
    # For now, return trivial parameterization
    # TODO: Implement full Poisson solve with seam constraints
    
    # Placeholder: Use face centroids projected to 2D
    u_coords = zeros(Float64, n_verts)
    v_coords = zeros(Float64, n_verts)
    
    for i in 1:n_verts
        p = vs[i]
        u_coords[i] = p[1]
        v_coords[i] = p[2]
    end
    
    seam_matching = Dict{Int,Int}()
    
    println("  [TODO] Full Poisson integration not yet implemented")
    println("  Returning placeholder (x,y) projection")
    
    return u_coords, v_coords, seam_matching
end


"""
    compute_seam_constraints(
        mesh, theta, period_jumps, cut_graph
    ) -> (constraint_matrix, constraint_rhs)

Compute constraints that ensure seamless parameterization across cuts.

For each seam edge:
- The jump in (u,v) must be an integer grid translation
- The rotation must be a multiple of 90° (determined by period jumps)

This gives linear constraints: (u_right - u_left, v_right - v_left) = (n, m) + R(period)
where R is rotation by period*90° and (n,m) are integers to be determined.
"""
function compute_seam_constraints(
    mesh::GeometryBasics.Mesh,
    theta::Vector{Float64},
    period_jumps::Dict{Tuple{Int,Int},Int},
    cut_graph::CutGraph
)
    # TODO: Implement seam constraints
    println("  Computing seam constraints...")
    
    # For each seam edge:
    # 1. Find the corresponding period jump
    # 2. Determine required rotation (0°, 90°, 180°, 270°)
    # 3. Add constraint: jump = integer + rotation
    
    # Return empty for now
    return sparse(Float64[], Int[], Int[], 0, 0), Float64[]
end


"""
    extract_quad_mesh(mesh, param::GlobalParameterization)
    -> (quad_vertices, quad_faces)

Extract a pure quad mesh from the parameterization by:
1. Rounding (u,v) to integer grid
2. Creating quads at each grid cell
3. Merging vertices that should be identified

This is the final output: a quad mesh suitable for simulation.
"""
function extract_quad_mesh(
    mesh::GeometryBasics.Mesh,
    param::GlobalParameterization
)
    println("\n" * "="^60)
    println("Quad Mesh Extraction")
    println("="^60)
    
    # TODO: Implement quad extraction
    println("  [TODO] Quad mesh extraction not yet implemented")
    
    # Placeholder
    quad_vertices = Point3f[]
    quad_faces = Vector{Int}[]
    
    return quad_vertices, quad_faces
end
