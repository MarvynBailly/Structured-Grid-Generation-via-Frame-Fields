# Frame Field Generation API Reference

Complete function documentation for the FrameField module.

## Main Functions

### `generate_frame_field`

```julia
generate_frame_field(
    mesh::Mesh;
    constrained_faces::Vector{Int}=Int[],
    constrained_angles::Vector{Float64}=Float64[],
    τ::Float64=1e-6,
    max_iterations::Int=10000,
    verbose::Bool=false
) -> FrameFieldResult
```

Generate a smooth frame field over a triangular mesh.

**Arguments:**
- `mesh::Mesh` - Input triangular mesh (GeometryBasics.Mesh)
- `constrained_faces::Vector{Int}` - Face indices with fixed orientations
- `constrained_angles::Vector{Float64}` - Corresponding angle values (radians)
- `τ::Float64` - Convergence tolerance for MIP solver (default: 1e-6)
- `max_iterations::Int` - Maximum iterations for MIP solver (default: 10000)
- `verbose::Bool` - Print progress information (default: false)

**Returns:**
- `FrameFieldResult` - Complete solution with angles, period jumps, and singularities

**Example:**
```julia
result = generate_frame_field(mesh, verbose=true)
println("Energy: ", result.energy)
```

---

## Data Structures

### `FrameFieldResult`

Container for frame field solution.

**Fields:**
- `angles::Vector{Float64}` - Face angles θ_f (length = number of faces)
- `period_jumps::Dict{Int, Int}` - Period jumps p_e indexed by edge
- `singularities::Vector{Tuple{Int, Float64}}` - (vertex_idx, index) pairs
- `energy::Float64` - Final smoothness energy value
- `converged::Bool` - Whether MIP solver converged
- `iterations::Int` - Number of MIP solver iterations

**Example:**
```julia
for (v_idx, index) in result.singularities
    println("Vertex $v_idx has index $index")
end
```

### `MeshTopology`

Complete mesh topology information.

**Fields:**
- `edge_map::Dict{Tuple{Int,Int}, Int}` - Maps vertex pairs to edge indices
- `dual_adj::Vector{Vector{Tuple{Int,Int}}}` - Face adjacency in dual graph
- `edge_to_faces::Dict{Int, Vector{Int}}` - Maps edges to incident faces
- `n_edges::Int` - Total number of edges
- `reference_edges::Vector{Int}` - Reference edge for each face
- `kappa::Dict{Tuple{Int,Int}, Float64}` - Angles between reference edges

### `DijkstraForest`

Spanning forest for dimensionality reduction.

**Fields:**
- `tree_edges::Set{Int}` - Edges in spanning forest
- `fixed_edges::Set{Int}` - Edges with p_e = 0
- `free_edges::Vector{Int}` - Edges with p_e as free variables
- `parent::Dict{Int, Int}` - Parent face in tree (-1 for roots)

---

## Mesh Topology Functions

### `build_mesh_topology`

```julia
build_mesh_topology(mesh::Mesh) -> MeshTopology
```

Construct complete mesh topology from triangular mesh.

**Process:**
1. Enumerate unique edges
2. Build face adjacency (dual graph)
3. Assign reference edges
4. Compute angles κ_ij between reference edges

**Example:**
```julia
topology = build_mesh_topology(mesh)
println("Edges: ", topology.n_edges)
```

### `get_face_neighbors`

```julia
get_face_neighbors(topology::MeshTopology, face_idx::Int) -> Vector{Int}
```

Get all neighboring face indices for a given face.

### `get_shared_edge`

```julia
get_shared_edge(topology::MeshTopology, face_i::Int, face_j::Int) -> Union{Int, Nothing}
```

Get edge index shared between two faces, or nothing if not adjacent.

### `compute_angle_defect`

```julia
compute_angle_defect(mesh::Mesh, vertex_idx::Int) -> Float64
```

Compute angle defect at a vertex: 2π - (sum of incident face angles).

---

## Dijkstra Forest Functions

### `build_dijkstra_forest`

```julia
build_dijkstra_forest(
    mesh::Mesh,
    topology::MeshTopology,
    constrained_faces::Vector{Int}=Int[]
) -> DijkstraForest
```

Build spanning forest over dual graph starting from constrained faces.

**Purpose:** Identify edges whose period jumps can be fixed to zero.

**Example:**
```julia
forest = build_dijkstra_forest(mesh, topology, [1, 2])
println("Free edges: ", length(forest.free_edges))
```

### `is_tree_edge`

```julia
is_tree_edge(forest::DijkstraForest, edge_idx::Int) -> Bool
```

Check if edge is in spanning forest.

### `get_free_variables`

```julia
get_free_variables(forest::DijkstraForest) -> Vector{Int}
```

Get list of edge indices with free period jump variables.

### `count_free_edges`

```julia
count_free_edges(forest::DijkstraForest) -> Int
```

Count number of free edge variables.

---

## Energy Functions

### `compute_smoothness_energy`

```julia
compute_smoothness_energy(
    mesh::Mesh,
    topology::MeshTopology,
    theta::Vector{Float64},
    p::Dict{Int, Int}
) -> Float64
```

Compute frame field smoothness energy: Σ (θ_i + κ_ij + (π/2)p_ij - θ_j)²

### `assemble_system_matrix`

```julia
assemble_system_matrix(
    mesh::Mesh,
    topology::MeshTopology,
    forest::DijkstraForest,
    constrained_faces::Vector{Int}=Int[],
    constrained_angles::Vector{Float64}=Float64[]
) -> (SparseMatrixCSC, Vector{Float64}, Dict, Dict)
```

Assemble linear system for frame field optimization.

**Returns:**
- `A::SparseMatrixCSC` - System matrix
- `b::Vector{Float64}` - Right-hand side
- `var_to_theta::Dict{Int, Int}` - Variable to face mapping
- `var_to_p::Dict{Int, Int}` - Variable to edge mapping

### `compute_singularity_index`

```julia
compute_singularity_index(
    mesh::Mesh,
    topology::MeshTopology,
    p::Dict{Int, Int},
    vertex_idx::Int
) -> Float64
```

Compute cross field index at a vertex: I(v) = I₀(v) + Σ(p_ij / 4)

---

## Singularity Functions

### `compute_singularities`

```julia
compute_singularities(
    mesh::Mesh,
    topology::MeshTopology,
    p::Dict{Int, Int};
    tol::Float64=1e-6
) -> Vector{Tuple{Int, Float64}}
```

Find all singular vertices (vertices with non-zero index).

**Returns:** List of (vertex_index, singularity_index) pairs

**Example:**
```julia
singularities = compute_singularities(mesh, topology, result.period_jumps)
for (v, idx) in singularities
    valence = 4 + round(Int, 4 * idx)
    println("Vertex $v: valence $valence")
end
```

### `get_frame_directions`

```julia
get_frame_directions(result::FrameFieldResult, face_idx::Int) -> Vector{Float64}
```

Get four frame directions at a given face.

**Returns:** [θ_f, θ_f + π/2, θ_f + π, θ_f + 3π/2]

---

## Visualization Functions

### `plot_frame_field`

```julia
plot_frame_field(
    mesh::Mesh,
    result::FrameFieldResult;
    arrow_scale::Float64=0.1,
    show_all_directions::Bool=false
)
```

Plot frame field on mesh (requires GLMakie or WGLMakie).

### `plot_singularities`

```julia
plot_singularities(
    mesh::Mesh,
    result::FrameFieldResult;
    marker_size::Float64=20.0,
    show_labels::Bool=true
)
```

Plot singular vertices on mesh.

### `plot_period_jumps`

```julia
plot_period_jumps(mesh::Mesh, topology::MeshTopology, result::FrameFieldResult)
```

Plot period jumps on mesh edges.

### `print_frame_field_summary`

```julia
print_frame_field_summary(result::FrameFieldResult)
```

Print text summary of frame field solution.

**Example:**
```julia
result = generate_frame_field(mesh)
print_frame_field_summary(result)
```

### `export_frame_field_vtk`

```julia
export_frame_field_vtk(filename::String, mesh::Mesh, result::FrameFieldResult)
```

Export frame field to VTK format for ParaView (requires WriteVTK).

---

## Usage Patterns

### Basic Usage

```julia
using FrameField
using GeometryBasics

# Load mesh
mesh = load_mesh("model.obj")

# Generate frame field
result = generate_frame_field(mesh, verbose=true)

# Analyze results
print_frame_field_summary(result)
plot_frame_field(mesh, result)
```

### With Constraints

```julia
# Constrain boundary faces
boundary_faces = [1, 2, 3, 4]
boundary_angles = [0.0, 0.0, π/2, π/2]

result = generate_frame_field(
    mesh,
    constrained_faces=boundary_faces,
    constrained_angles=boundary_angles,
    verbose=true
)
```

### Custom Analysis

```julia
# Check specific face
face_idx = 10
dirs = get_frame_directions(result, face_idx)
println("Directions: ", dirs)

# Analyze singularities
for (v_idx, index) in result.singularities
    if abs(index) > 0.3
        println("Strong singularity at vertex $v_idx")
    end
end
```

---

## Performance Tips

1. **Mesh Quality:** Well-shaped triangles improve convergence
2. **Constraints:** Fewer constraints = faster solve
3. **Tolerance:** Looser tolerance (1e-4) speeds up MIP solver
4. **Sparse Systems:** Solver exploits sparsity automatically
5. **Large Meshes:** Consider iterative refinement or domain decomposition

---

## Error Handling

Common issues and solutions:

**"Matrix A must have non-zero diagonal"**
- Mesh has degenerate triangles
- Solution: Clean mesh, remove zero-area faces

**"MIP solver did not converge"**
- System is ill-conditioned
- Solution: Increase max_iterations or relax tolerance

**"No free variables found"**
- All faces are constrained
- Solution: Remove some constraints or check mesh topology

---

## See Also

- [Algorithm Details](algorithm_details.md) - Detailed algorithm explanation
- [README](../README.md) - Module overview
- [Examples](../examples/) - Usage examples
- [GreedyMIPSolver](../../greedy_mip_solver/) - Mixed-integer solver documentation
