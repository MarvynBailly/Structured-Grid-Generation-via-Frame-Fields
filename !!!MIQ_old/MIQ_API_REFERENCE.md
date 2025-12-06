# MIQ (Mixed-Integer Quadrangulation) API Reference

This document provides a comprehensive API reference for all grid functions, geometry operations, and mesh topology utilities required to perform the MIQ procedure for cross-field generation and quadrangulation.

---

## Table of Contents

1. [Data Structures](#data-structures)
2. [Mesh I/O](#mesh-io)
3. [Mesh Generation & Preprocessing](#mesh-generation--preprocessing)
4. [Topology Construction](#topology-construction)
5. [Geometry & Local Frames](#geometry--local-frames)
6. [Cross Field Computation](#cross-field-computation)
7. [Spanning Forest](#spanning-forest)
8. [System Assembly & Solving](#system-assembly--solving)
9. [Singularity Detection](#singularity-detection)
10. [Visualization](#visualization)
11. [Utility Functions](#utility-functions)

---

## Data Structures

### `MeshTopology`

A struct that encapsulates the dual graph structure and edge topology of a triangle mesh.

```julia
struct MeshTopology
    edge_map::Dict{Tuple{Int, Int}, Int}        # Maps sorted vertex pairs (v1, v2) to unique edge index
    dual_adj::Vector{Vector{Tuple{Int, Int}}}   # Maps face index to list of (neighbor_face, shared_edge)
    edge_to_faces::Dict{Int, Vector{Int}}       # Maps edge index to list of adjacent face indices
    n_edges::Int                                 # Total number of edges in the mesh
end
```

**Fields:**
- `edge_map`: Dictionary mapping sorted vertex tuples to unique edge IDs (1-based)
- `dual_adj`: Adjacency list for the dual graph, indexed by face ID
- `edge_to_faces`: Maps each edge to the faces it connects (1 or 2 faces)
- `n_edges`: Total count of unique edges in the mesh

---

## Mesh I/O

### `load_triangulation`

Loads a mesh from file and standardizes it into a `GeometryBasics.Mesh` with 1-based integer `TriangleFace` indices.

```julia
load_triangulation(filepath::String) -> Mesh
```

**Arguments:**
- `filepath::String`: Path to the mesh file (.msh, .obj, .stl, etc.)

**Returns:**
- `Mesh`: A `GeometryBasics.Mesh` containing `Point3f` vertices and `TriangleFace{Int}` faces

**Throws:**
- Error if file not found

**Example:**
```julia
mesh = load_triangulation("output/meshes/hole.msh")
```

---

### `save_cross_field_to_vtu`

Exports cross field results to VTU format for visualization in ParaView.

```julia
save_cross_field_to_vtu(
    filename::String, 
    mesh::Mesh, 
    thetas::Vector{Float64}, 
    constraint_vecs::Dict{Int, Vec3f}, 
    topo::MeshTopology, 
    p_map::Dict{Int, Int}, 
    x_sol::Vector{Float64}
)
```

**Arguments:**
- `filename`: Output path (without .vtu extension)
- `mesh`: The triangle mesh
- `thetas`: Per-face angle values (in radians)
- `constraint_vecs`: User-specified direction constraints (face_id → direction vector)
- `topo`: Mesh topology structure
- `p_map`: Mapping from edge indices to matrix variable indices for period jumps
- `x_sol`: Solution vector from the optimization

**Exports:**
- Cross field directions (4 vectors per face)
- Constraint markers and desired directions
- Period jumps (singularity indicators)
- Face angles (theta values)

**Example:**
```julia
save_cross_field_to_vtu("output/cross_field", mesh, thetas, constraints, topo, p_map, x_sol)
```

---

## Mesh Generation & Preprocessing

### `generate_square_mesh`

Generates a clean, connected triangulated unit square mesh [0,1]×[0,1].

```julia
generate_square_mesh(; nx=6, ny=6) -> Mesh
```

**Keyword Arguments:**
- `nx::Int`: Number of subdivisions along X-axis (default: 6)
- `ny::Int`: Number of subdivisions along Y-axis (default: 6)

**Returns:**
- `Mesh`: Triangle mesh with counter-clockwise face winding

**Notes:**
- Each quad is split into 2 triangles
- Total vertices: `(nx+1) × (ny+1)`
- Total faces: `2 × nx × ny`

**Example:**
```julia
mesh = generate_square_mesh(nx=10, ny=10)  # 10×10 grid → 200 triangles
```

---

### `weld_mesh`

Merges duplicate vertices within a tolerance and updates face connectivity.

```julia
weld_mesh(mesh::Mesh; tol=1e-5) -> Mesh
```

**Arguments:**
- `mesh::Mesh`: Input mesh potentially with duplicate vertices
- `tol::Float64`: Distance threshold for considering vertices identical (default: 1e-5)

**Returns:**
- `Mesh`: New mesh with welded vertices

**Algorithm:**
- Uses KDTree for efficient nearest-neighbor search
- Remaps face indices to canonical vertex indices
- Removes degenerate faces (where 2+ vertices collapse)

**Example:**
```julia
welded_mesh = weld_mesh(raw_mesh; tol=1e-6)
```

---

## Topology Construction

### `build_topology`

Constructs the dual graph adjacency structure and assigns unique edge IDs.

```julia
build_topology(mesh::Mesh) -> MeshTopology
```

**Arguments:**
- `mesh::Mesh`: Triangle mesh with 1-based face indices

**Returns:**
- `MeshTopology`: Complete topology structure

**Algorithm:**
1. Iterates over all triangles to identify unique edges
2. Assigns sequential edge IDs (1-based)
3. Builds dual adjacency by connecting faces sharing edges
4. Maps edges to their incident faces

**Example:**
```julia
topo = build_topology(mesh)
println("Mesh has $(topo.n_edges) edges")
```

---

## Geometry & Local Frames

### `get_local_frame`

Computes a local 2D coordinate system (X, Y) in 3D space for a triangle.

```julia
get_local_frame(p1::Point3f, p2::Point3f, p3::Point3f) -> (center, x_axis, y_axis)
```

**Arguments:**
- `p1, p2, p3`: Vertices of the triangle in 3D

**Returns:**
- `center::Point3f`: Centroid of the triangle
- `x_axis::Vec3f`: Normalized direction from p1 to p2
- `y_axis::Vec3f`: In-plane direction perpendicular to x_axis

**Notes:**
- X-axis aligned with edge p1→p2
- Y-axis computed via cross product to ensure orthogonality
- Forms a right-handed coordinate system

**Example:**
```julia
center, x, y = get_local_frame(p1, p2, p3)
angle = atan(dot(direction, y), dot(direction, x))  # Project 3D direction to 2D angle
```

---

## Cross Field Computation

### `compute_kappas`

Computes the transport angle (κ_ij) for every interior edge, measuring the rotation between adjacent face frames.

```julia
compute_kappas(mesh::Mesh, topo::MeshTopology) -> Dict{Int, Float64}
```

**Arguments:**
- `mesh::Mesh`: Triangle mesh
- `topo::MeshTopology`: Topology structure

**Returns:**
- `Dict{Int, Float64}`: Map from edge_id to transport angle in radians (wrapped to [-π, π])

**Algorithm:**
1. For each interior edge connecting faces i and j
2. Compute local frames for both faces
3. Find shared edge vector
4. Measure angle of shared edge in each local frame
5. κ_ij = angle_j - angle_i (coordinate system mismatch)

**Notes:**
- Only computes for interior edges (connecting 2 faces)
- Angles wrapped to [-π, π] for consistency

**Example:**
```julia
kappas = compute_kappas(mesh, topo)
for (edge_id, kappa) in kappas
    println("Edge $edge_id: transport angle = $(rad2deg(kappa))°")
end
```

---

### `constraints_to_angles`

Converts 3D direction constraint vectors to 2D angles in each face's local frame.

```julia
constraints_to_angles(mesh::Mesh, constraints::Dict{Int, Vec3f}) -> Dict{Int, Float64}
```

**Arguments:**
- `mesh::Mesh`: Triangle mesh
- `constraints::Dict{Int, Vec3f}`: Map from face_id to desired 3D direction vector

**Returns:**
- `Dict{Int, Float64}`: Map from face_id to local angle in radians

**Algorithm:**
1. For each constrained face, compute its local frame
2. Project the 3D constraint vector onto the local (X, Y) plane
3. Compute angle: θ = atan2(local_y, local_x)

**Example:**
```julia
constraint_vecs = Dict(1 => Vec3f(1, 0, 0), 4 => Vec3f(0, 1, 0))
constrained_angles = constraints_to_angles(mesh, constraint_vecs)
```

---

## Spanning Forest

### `compute_spanning_forest`

Computes a spanning forest of the dual graph, identifying edges where period jumps should be fixed to zero.

```julia
compute_spanning_forest(
    topology::MeshTopology, 
    n_faces::Int, 
    constrained_faces::Vector{Int}
) -> BitVector
```

**Arguments:**
- `topology::MeshTopology`: Mesh topology structure
- `n_faces::Int`: Total number of faces in the mesh
- `constrained_faces::Vector{Int}`: List of face IDs with user constraints (roots of trees)

**Returns:**
- `BitVector`: Length `n_edges`, true indicates edge is in the spanning forest (p_ij = 0)

**Algorithm:**
- Uses BFS (Breadth-First Search) starting from constrained faces
- Grows trees by visiting unvisited neighbors via dual graph
- Marks traversed edges as part of the forest
- Handles disconnected components if no constraints in a region

**Notes:**
- Ensures the system has a unique solution by fixing a maximal set of period jumps
- Standard meshes are usually one connected component

**Example:**
```julia
constrained_faces = [1, 10, 25]
fixed_edges = compute_spanning_forest(topo, length(faces(mesh)), constrained_faces)
println("Fixed $(sum(fixed_edges)) edges in spanning forest")
```

---

## System Assembly & Solving

### `assemble_system`

Assembles the sparse linear system for the mixed-integer quadrangulation optimization.

```julia
assemble_system(
    mesh::Mesh, 
    topo::MeshTopology, 
    fixed_edges::BitVector, 
    kappas::Dict{Int, Float64}, 
    constrained_angles::Dict{Int, Float64}
) -> (A, b, theta_map, p_map)
```

**Arguments:**
- `mesh`: Triangle mesh
- `topo`: Topology structure
- `fixed_edges`: Boolean vector indicating which edges have p_ij = 0
- `kappas`: Transport angles for each edge
- `constrained_angles`: Fixed angle values for constrained faces

**Returns:**
- `A::SparseMatrixCSC`: Symmetric positive-definite system matrix
- `b::Vector{Float64}`: Right-hand side vector
- `theta_map::Dict{Int, Int}`: Map from free face_id to variable index in solution vector
- `p_map::Dict{Int, Int}`: Map from free edge_id to variable index in solution vector

**Energy Function:**
```
E = Σ (θ_j - θ_i - κ_ij - (π/2)p_ij)²
```

**Variables:**
- θ_i: Angle at face i (continuous, free or constrained)
- p_ij: Period jump across edge ij (integer, free or fixed)

**Constraints:**
- Constrained faces: θ_i fixed to user-specified value
- Spanning forest edges: p_ij = 0
- Edges between two constrained faces: p_ij = round((2/π)(θ_j - θ_i - κ_ij))

**Notes:**
- System is derived from energy gradient (∂E/∂θ_i, ∂E/∂p_ij = 0)
- Only free variables included in the system
- Fixed variables moved to the right-hand side

**Example:**
```julia
A, b, theta_map, p_map = assemble_system(mesh, topo, fixed_edges, kappas, constrained_angles)
println("System size: $(size(A, 1)) variables")
```

---

### `solve_greedy`

Solves the mixed-integer optimization using greedy rounding with local relaxation.

```julia
solve_greedy(A::SparseMatrixCSC, b::Vector{Float64}, p_var_map::Dict{Int, Int}) -> Vector{Float64}
```

**Arguments:**
- `A`: Sparse system matrix (SPD)
- `b`: Right-hand side vector
- `p_var_map`: Mapping from edge IDs to integer variable indices

**Returns:**
- `Vector{Float64}`: Solution vector (continuous θ values + integer p values)

**Algorithm:**
1. Factorize system: `x = A \ b` (initial continuous solution)
2. Iteratively:
   - Find p variable closest to an integer
   - Round it and fix it
   - Locally relax neighboring variables using Gauss-Seidel
3. Repeat until all p variables are integers

**Notes:**
- Uses Cholesky factorization (requires SPD matrix)
- Local relaxation prevents error propagation

**Example:**
```julia
x_sol = solve_greedy(A, b, p_map)
```

---

### `solve_greedy_global`

Enhanced greedy solver with global relaxation for better convergence.

```julia
solve_greedy_global(A::SparseMatrixCSC, b::Vector{Float64}, p_var_map::Dict{Int, Int}) -> Vector{Float64}
```

**Arguments:**
- Same as `solve_greedy`

**Returns:**
- `Vector{Float64}`: Solution vector with integer period jumps

**Algorithm:**
1. Solve initial continuous problem
2. Iteratively:
   - Round closest-to-integer p variable
   - Relax **all** free variables (not just local neighborhood)
   - Uses queue-based Gauss-Seidel for efficiency
3. Repeat until all p variables are integers

**Advantages over `solve_greedy`:**
- More robust for complex meshes
- Better handles long-range coupling
- Logs progress for debugging

**Example:**
```julia
x_sol = solve_greedy_global(A, b, p_map)
```

---

### `local_gauss_seidel!`

Performs local Gauss-Seidel relaxation starting from a single node.

```julia
local_gauss_seidel!(
    A::SparseMatrixCSC, 
    b::Vector{Float64}, 
    x::Vector{Float64}, 
    fixed_mask::BitVector, 
    start_node::Int; 
    max_iter=500000, 
    tol=1e-5
)
```

**Arguments:**
- `A`: Sparse system matrix
- `b`: Right-hand side
- `x`: Solution vector (modified in-place)
- `fixed_mask`: Boolean mask indicating which variables are fixed
- `start_node`: Starting variable index for relaxation
- `max_iter`: Maximum iterations (default: 500,000)
- `tol`: Convergence tolerance (default: 1e-5)

**Notes:**
- Uses queue-based propagation (only updates variables with large residuals)
- Stops when all residuals below tolerance
- Modifies `x` in-place

---

### `run_gauss_seidel_queue!`

Queue-based Gauss-Seidel relaxation for global updates.

```julia
run_gauss_seidel_queue!(
    A::SparseMatrixCSC, 
    b::Vector{Float64}, 
    x::Vector{Float64}, 
    fixed_mask::BitVector, 
    q::Queue{Int}, 
    in_queue::BitVector; 
    max_iter=1_000_000, 
    tol=1e-5
) -> Int
```

**Arguments:**
- `A, b, x, fixed_mask`: Same as above
- `q`: Queue of variable indices to process
- `in_queue`: Boolean array tracking queue membership
- `max_iter`: Maximum iterations
- `tol`: Residual tolerance

**Returns:**
- `Int`: Number of iterations performed

**Algorithm:**
1. Pop variable from queue
2. Compute residual: r = b[k] - Σ A[k,j] x[j]
3. If |r| > tol: update x[k], add neighbors to queue
4. Repeat until queue empty or max_iter reached

---

## Singularity Detection

### `compute_singularities`

Identifies singularities (irregular vertices) in the cross field by summing period jumps around 1-rings.

```julia
compute_singularities(
    mesh::Mesh, 
    topo::MeshTopology, 
    p_values::Dict{Int, Int}
) -> Dict{Int, Float64}
```

**Arguments:**
- `mesh`: Triangle mesh
- `topo`: Topology structure
- `p_values`: Period jump values (edge_id → integer)

**Returns:**
- `Dict{Int, Float64}`: Map from vertex_id to singularity index

**Singularity Indices:**
- `+0.25`: Valence-3 singularity (triangle vertex in quad mesh)
- `-0.25`: Valence-5 singularity (pentagon vertex in quad mesh)
- `±0.5, ±0.75, ...`: Higher-order singularities (rare)

**Algorithm:**
1. Build vertex-to-faces adjacency
2. For each interior vertex:
   - Sort incident faces in CCW order (geometric 1-ring)
   - Sum oriented period jumps around the loop
3. Singularity index = Σp / 4

**Notes:**
- Only reports vertices with non-zero index
- Boundary vertices are skipped (open 1-ring)
- Uses robust geometric sorting to handle non-manifold cases

**Example:**
```julia
singularities = compute_singularities(mesh, topo, final_p_values)
for (v_id, index) in singularities
    println("Vertex $v_id: index = $index")
end
```

---

### `get_ordered_1ring`

Sorts faces incident to a vertex in counter-clockwise order by following half-edge connectivity.

```julia
get_ordered_1ring(
    pivot_v::Int, 
    adj_faces::Vector{Int}, 
    all_faces::Vector{TriangleFace{Int}}, 
    topo::MeshTopology
) -> Vector{Int}
```

**Arguments:**
- `pivot_v`: Vertex ID
- `adj_faces`: List of face IDs incident to the vertex
- `all_faces`: Complete list of mesh faces
- `topo`: Topology structure

**Returns:**
- `Vector{Int}`: Sorted face IDs in CCW order, or empty vector if boundary vertex

**Notes:**
- Uses geometric traversal (not combinatorial)
- Detects boundary by checking for open loops
- Returns empty vector for non-manifold or boundary vertices

---

### `find_shared_edge`

Finds the edge ID shared by two adjacent faces.

```julia
find_shared_edge(f1::Int, f2::Int, topo::MeshTopology) -> Int
```

**Arguments:**
- `f1, f2`: Face IDs
- `topo`: Topology structure

**Returns:**
- `Int`: Edge ID if faces are adjacent, -1 otherwise

**Example:**
```julia
edge_id = find_shared_edge(5, 12, topo)
if edge_id != -1
    println("Faces 5 and 12 share edge $edge_id")
end
```

---

## Visualization

### `compute_face_centroids`

Computes the centroid of every triangle in the mesh.

```julia
compute_face_centroids(mesh::Mesh) -> Vector{Point3f}
```

**Arguments:**
- `mesh::Mesh`: Triangle mesh

**Returns:**
- `Vector{Point3f}`: Centroid of each face (indexed by face ID)

**Example:**
```julia
centroids = compute_face_centroids(mesh)
```

---

### `get_dual_graph_segments`

Generates line segments representing the dual graph edges.

```julia
get_dual_graph_segments(topo::MeshTopology, centroids::Vector{Point3f}) -> Vector{Point3f}
```

**Arguments:**
- `topo`: Topology structure
- `centroids`: Face centroids

**Returns:**
- `Vector{Point3f}`: Flat array of point pairs [p1, p2, p1, p2, ...] for line rendering

**Notes:**
- Each segment connects centroids of adjacent faces
- Avoids duplicates (only includes edge once)

**Example:**
```julia
segments = get_dual_graph_segments(topo, centroids)
linesegments!(ax, segments, color=:blue)
```

---

### `visualize_forest`

Visualizes the spanning forest with color-coded trees.

```julia
visualize_forest(
    mesh::Mesh, 
    topo::MeshTopology, 
    fixed_edges::BitVector, 
    constrained_faces::Vector{Int}
) -> Figure
```

**Arguments:**
- `mesh`: Triangle mesh
- `topo`: Topology structure
- `fixed_edges`: Boolean vector indicating spanning forest edges
- `constrained_faces`: Root faces (constraints)

**Returns:**
- `Figure`: Makie figure with visualization

**Visualization:**
- Each spanning tree colored uniquely
- Non-tree faces in gray
- Constrained roots highlighted in green
- Dual graph edges shown in blue

---

### `plot_cross_field`

Visualizes the generated cross field with optional singularities and constraints.

```julia
plot_cross_field(
    mesh::Mesh, 
    thetas::Vector{Float64}, 
    constrained_faces::Vector{Int}; 
    constraint_vecs::Union{Nothing, Dict{Int, Vec3f}}=nothing, 
    sings::Union{Nothing, Dict{Int, Float64}}=nothing,
    topo::Union{Nothing, MeshTopology}=nothing,
    p_values::Union{Nothing, Dict{Int, Int}}=nothing
) -> Figure
```

**Arguments:**
- `mesh`: Triangle mesh
- `thetas`: Per-face angle values
- `constrained_faces`: List of constrained face IDs
- `constraint_vecs`: User direction constraints (optional)
- `sings`: Singularities from `compute_singularities` (optional)
- `topo`: Topology structure (optional, for dual edge visualization)
- `p_values`: Period jumps (optional, for dual edge coloring)

**Returns:**
- `Figure`: Interactive 3D visualization

**Visualization Elements:**
- Mesh wireframe (gray)
- Cross field (4 black lines per face at ±π/2)
- User constraints (thick red arrows)
- Singularities (colored markers: green for +, blue for -)
- Dual edges (colored by period jump: red for +, blue for -, gray for 0)

**Example:**
```julia
fig = plot_cross_field(mesh, thetas, [1, 5]; 
                       constraint_vecs=constraints, 
                       sings=singularities,
                       topo=topo,
                       p_values=final_p_values)
display(fig)
```

---

### `plot_field_smoothness`

Visualizes the smoothness of the cross field by showing alignment errors.

```julia
plot_field_smoothness(
    mesh::Mesh, 
    thetas::Vector{Float64}, 
    topo::MeshTopology, 
    kappas::Dict{Int, Float64}, 
    p_values::Dict{Int, Int}
) -> Figure
```

**Arguments:**
- `mesh`: Triangle mesh
- `thetas`: Per-face angle values
- `topo`: Topology structure
- `kappas`: Transport angles
- `p_values`: Period jumps

**Returns:**
- `Figure`: Heat map visualization (white = smooth, red = error)

**Algorithm:**
- For each face, computes maximum mismatch with neighbors
- Mismatch = |θ_j - θ_i - κ_ij - (π/2)p_ij|
- Colors faces by maximum error

**Example:**
```julia
fig_smooth = plot_field_smoothness(mesh, thetas, topo, kappas, final_p_values)
```

---

## Utility Functions

### `log_euler_char`

Computes and logs the Euler characteristic (χ = V - E + F).

```julia
log_euler_char(mesh::Mesh, topology::MeshTopology)
```

**Arguments:**
- `mesh`: Triangle mesh
- `topology`: Topology structure

**Prints:**
- Number of vertices (V)
- Number of edges (E)
- Number of faces (F)
- Euler characteristic (χ)

**Notes:**
- For a closed orientable surface: χ = 2(1 - g), where g is genus
- Sphere: χ = 2
- Torus: χ = 0

**Example:**
```julia
log_euler_char(mesh, topo)
# Output: |V| = 642, |E| = 1920, |F| = 1280, χ = 2
```

---

## Complete Workflow Example

```julia
using GeometryBasics
using LinearAlgebra

# 1. Load or generate mesh
mesh = load_triangulation("output/meshes/hole.msh")
# OR
# mesh = generate_square_mesh(nx=10, ny=10)

# 2. Build topology
topo = build_topology(mesh)
log_euler_char(mesh, topo)

# 3. Define constraints
n_constraints = 3
total_faces = length(faces(mesh))
constrained_faces = rand(1:total_faces, n_constraints)
constraint_vecs = Dict(f => normalize(Vec3f(1, 0, 0)) for f in constrained_faces)

# 4. Compute geometry
kappas = compute_kappas(mesh, topo)
constrained_angles = constraints_to_angles(mesh, constraint_vecs)

# 5. Compute spanning forest
fixed_edges = compute_spanning_forest(topo, length(faces(mesh)), constrained_faces)

# 6. Assemble system
A, b, theta_map, p_map = assemble_system(mesh, topo, fixed_edges, kappas, constrained_angles)

# 7. Solve
x_sol = solve_greedy_global(A, b, p_map)

# 8. Extract results
final_thetas = zeros(length(faces(mesh)))
final_p_values = Dict{Int, Int}()

for (e_idx, mat_idx) in p_map
    final_p_values[e_idx] = round(Int, x_sol[mat_idx])
end

for (f_idx, mat_idx) in theta_map
    final_thetas[f_idx] = x_sol[mat_idx]
end

for (f_idx, val) in constrained_angles
    final_thetas[f_idx] = val
end

# 9. Compute singularities
singularities = compute_singularities(mesh, topo, final_p_values)
println("Found $(length(singularities)) singularities")

# 10. Visualize
fig = plot_cross_field(mesh, final_thetas, constrained_faces; 
                       constraint_vecs=constraint_vecs, 
                       sings=singularities,
                       topo=topo,
                       p_values=final_p_values)
display(fig)

# 11. Export
save_cross_field_to_vtu("output/cross_field", mesh, final_thetas, 
                        constraint_vecs, topo, p_map, x_sol)
```

---

## Dependencies

- **GeometryBasics**: Mesh data structures (`Mesh`, `Point3f`, `TriangleFace`, `Vec3f`)
- **LinearAlgebra**: Vector/matrix operations, `normalize`, `dot`, `cross`
- **SparseArrays**: Sparse matrix assembly (`sparse`, `rowvals`, `nonzeros`)
- **DataStructures**: Queue for BFS (`Queue`, `enqueue!`, `dequeue!`)
- **NearestNeighbors**: KDTree for vertex welding (`KDTree`, `nn`)
- **WriteVTK**: VTU file export (`vtk_grid`, `MeshCell`)
- **FileIO / MeshIO**: Mesh loading (`load`, `decompose`)
- **WGLMakie**: Visualization (`Figure`, `Axis3`, `mesh!`, `linesegments!`, `scatter!`)

---

## References

- Bommes et al. (2009): "Mixed-Integer Quadrangulation"
- Palacios & Zhang (2007): "Rotational Symmetry Field Design on Surfaces"
- Knöppel et al. (2013): "Globally Optimal Direction Fields"

---

## Notes

- All indices are **1-based** (Julia convention)
- Angles are in **radians** unless specified
- Face winding should be **counter-clockwise** for correct normal computation
- The system matrix is **symmetric positive-definite** (SPD)
- Period jumps (p_ij) must be **integers** for a valid quad mesh
- Singularities are unavoidable on surfaces with non-zero genus or boundary

