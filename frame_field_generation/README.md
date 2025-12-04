# Frame Field Generation

This module implements smooth frame field generation over triangular meshes for quadrilateral mesh generation, based on the methodology described in Bommes et al. (2009).

## Overview

A frame field is a 4-symmetry direction field (cross field) where at each face in the triangular mesh, four orthogonal directions are defined. The frame field is represented by:
- **Angle field** θ: F → ℝ - assigns an angle θ_f to each face f
- **Period-jump field** p: E → ℤ - assigns an integer period-jump p_e to each edge e

## Mathematical Formulation

### Smoothness Energy

The smoothness of a frame field is measured by:

```
E_smooth = Σ (θ_i + κ_ij + (π/2)p_ij - θ_j)²
```

where:
- θ_i, θ_j are face angles relative to reference edges
- κ_ij ∈ (-π, π] is the angle between reference edges
- p_ij is the integer period-jump across edge e_ij
- Sum is over all edges in the mesh

### Singularities

The cross field index at vertex v is:

```
I(v) = I₀(v) + Σ(p_ij/4)
```

where:
- I₀(v) = (1/2π)(A_d(v) + Σκ_ij) is the geometric index
- A_d(v) is the angle defect at vertex v
- Vertices with non-zero index become singularities in the frame field

### Optimization Problem

Given:
- Triangular mesh M = (V, E, F)
- Constrained faces F_c ⊂ F with fixed angles θ̂_i

Find θ and p that minimize E_smooth subject to θ_i = θ̂_i for i ∈ F_c.

This is a **Mixed-Integer Problem**: θ values are real, p values are integers.

## Module Structure

```
frame_field_generation/
├── README.md                    # This file
├── src/
│   ├── FrameField.jl           # Main module
│   ├── mesh_topology.jl        # Mesh topology utilities
│   ├── dijkstra_forest.jl      # Dijkstra tree for fixing period jumps
│   ├── energy.jl               # Energy computation and derivatives
│   ├── frame_field_solver.jl   # Frame field optimization
│   └── visualization.jl        # Plotting and visualization
├── test/
│   ├── test_frame_field.jl     # Unit tests
│   └── test_meshes/            # Test mesh files
├── examples/
│   ├── simple_square.jl        # Basic square mesh example
│   ├── constrained_field.jl    # Example with boundary constraints
│   └── singularity_analysis.jl # Singularity placement example
└── docs/
    ├── algorithm_details.md     # Detailed algorithm explanation
    └── api_reference.md         # Function documentation
```

## Algorithm Overview

1. **Construct Mesh Topology**
   - Build edge map and dual graph adjacency
   - Identify face neighbors across each edge

2. **Build Dijkstra Forest**
   - Starting from constrained faces
   - Identify edges whose period jumps can be fixed to zero
   - Reduces problem dimensionality without changing energy

3. **Set Up Optimization Problem**
   - Compute derivatives: ∂E/∂θ_k = 0, ∂E/∂p_ij = 0
   - Assemble linear system for free variables
   - Integer variables: period jumps on non-tree edges
   - Real variables: face angles on non-constrained faces

4. **Solve Mixed-Integer Problem**
   - Use greedy mixed-integer solver
   - Returns optimal θ and p values

5. **Extract Frame Field**
   - Compute four directions per face: θ_f + k·π/2 for k=0,1,2,3
   - Verify smoothness and singularity placement

## Usage

```julia
using FrameField

# Load or generate triangular mesh
mesh = generate_square_mesh(nx=10, ny=10)

# Optional: Set boundary constraints
constrained_faces = [1, 2, 3]  # Face indices
constrained_angles = [0.0, π/4, π/2]  # Desired angles

# Generate smooth frame field
result = generate_frame_field(
    mesh, 
    constrained_faces=constrained_faces,
    constrained_angles=constrained_angles
)

# Access results
θ = result.angles          # Face angles
p = result.period_jumps    # Edge period jumps
singularities = result.singularities  # Singular vertex indices

# Visualize
plot_frame_field(mesh, result)
plot_singularities(mesh, result)
```

## Dependencies

- `LinearAlgebra` - Matrix operations
- `SparseArrays` - Sparse matrix support
- `DataStructures` - Queue, priority queue for Dijkstra
- `GeometryBasics` - Mesh data structures
- `Statistics` - Basic statistics
- `Graphs` (optional) - Graph algorithms for dual graph
- Plotting: `GLMakie` or `WGLMakie` for visualization

## Integration with Greedy MIP Solver

This module uses the `GreedyMIPSolver` module from `../greedy_mip_solver/` to solve the mixed-integer optimization problem. The frame field solver:

1. Constructs the system matrix A and RHS vector b
2. Identifies which variables are integers (period jumps)
3. Calls `solve_greedy_mip(A, b, k)` where k is the number of period jump variables
4. Unpacks the solution into θ and p arrays

## Key Features

- ✅ Smooth frame field generation with minimal energy
- ✅ Support for orientation constraints on boundary faces
- ✅ Automatic singularity placement based on mesh geometry
- ✅ Efficient mixed-integer solver integration
- ✅ Comprehensive visualization tools
- ✅ Modular design for easy extension

## References

- Bommes, D., Zimmer, H., & Kobbelt, L. (2009). Mixed-integer quadrangulation. ACM Transactions on Graphics (TOG), 28(3), 1-10.
- Ray, N., Vallet, B., Li, W. C., & Lévy, B. (2008). N-symmetry direction field design. ACM Transactions on Graphics (TOG), 27(2), 1-13.

## Testing

Run the test suite:

```julia
using Pkg
Pkg.activate(".")
Pkg.test("FrameField")
```

Or run individual test files:

```julia
include("test/test_frame_field.jl")
```

## Status

- [x] Project structure created
- [ ] Mesh topology implementation
- [ ] Dijkstra forest implementation
- [ ] Energy computation
- [ ] Frame field solver
- [ ] Visualization tools
- [ ] Test suite
- [ ] Documentation
- [ ] Example scripts

## Future Enhancements

- Support for anisotropic frame fields
- GPU acceleration for large meshes
- Interactive constraint specification
- Integration with quad mesh extraction
- Mesh quality metrics and analysis
