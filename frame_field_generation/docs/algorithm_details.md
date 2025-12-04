# Frame Field Algorithm Details

This document provides a detailed explanation of the frame field generation algorithm implemented in this module.

## Background

Frame fields are N-symmetry direction fields where N directions are defined at each point, symmetric under rotation by 2π/N. For quadrilateral mesh generation, we use 4-symmetry frame fields (cross fields) where N=4.

## Mathematical Foundation

### Frame Field Representation

A frame field over a triangular mesh is represented by two functions:

1. **Angle field** θ: F → ℝ
   - Maps each face f to an angle θ_f
   - Represents rotation from a reference direction

2. **Period-jump field** p: E → ℤ
   - Maps each edge e to an integer p_e
   - Determines how frame directions align across edges

### Four Directions per Face

At each face f, the four frame directions are:
```
d₁ = θ_f
d₂ = θ_f + π/2
d₃ = θ_f + π
d₄ = θ_f + 3π/2
```

These are obtained by rotating a reference edge by θ_f, then adding multiples of π/2.

### Cross-Edge Consistency

For adjacent faces f_i and f_j sharing edge e_ij, the frame fields are related by:
```
θ_j = θ_i + κ_ij + (π/2)p_ij
```

where:
- κ_ij ∈ (-π, π] is the angle between reference edges
- p_ij ∈ ℤ is the period jump

This ensures that direction d₁ in face f_i aligns with direction d_{1+p_ij} in face f_j.

## Energy Minimization

### Smoothness Energy

The smoothness of a frame field is measured by:

```
E_smooth = Σ (θ_i + κ_ij + (π/2)p_ij - θ_j)²
```

where the sum is over all interior edges.

**Goal:** Find θ and p that minimize E_smooth subject to boundary constraints.

### Mixed-Integer Problem

This is a mixed-integer optimization problem because:
- θ values are real-valued (continuous)
- p values are integers (discrete)

General mixed-integer problems are NP-hard, so we use an approximation algorithm.

## Algorithm Steps

### 1. Mesh Topology Construction

**Input:** Triangular mesh M = (V, E, F)

**Process:**
- Enumerate all edges
- Build dual graph (face adjacency)
- Assign reference edges (one per face)
- Compute κ_ij for each edge

**Output:** MeshTopology structure

**Complexity:** O(|E| + |F|)

### 2. Dijkstra Forest Construction

**Purpose:** Reduce problem dimensionality by fixing some period jumps to zero.

**Key Insight:** Period jumps on spanning tree edges can be set to zero without changing the energy minimum.

**Process:**
1. Start from constrained faces (or face 1 if none)
2. Run Dijkstra's algorithm on dual graph
3. Mark spanning tree edges
4. Fix p_e = 0 for all tree edges
5. Only non-tree edges have p_e as free variables

**Output:** DijkstraForest structure

**Complexity:** O(|F| log |F| + |E|)

**Why it works:** 
- The energy E_smooth depends only on differences θ_i - θ_j
- Adding a constant to all θ values doesn't change the energy
- This rotational invariance allows us to fix one period jump per tree path

### 3. System Assembly

**Process:** Take derivatives of E_smooth with respect to each free variable.

**For period jumps (∂E/∂p_ij = 0):**
```
π(θ_i + κ_ij + (π/2)p_ij - θ_j) = 0
```

**For face angles (∂E/∂θ_k = 0):**
```
Σ 2(θ_k + κ_kj + (π/2)p_kj - θ_j) = 0
```

where sum is over neighbors of face k.

**Matrix Form:** Ax = b

where:
- First k variables are period jumps (integers)
- Remaining variables are face angles (reals)
- A is sparse, symmetric, positive semi-definite

**Output:** Sparse matrix A, vector b, variable mappings

**Complexity:** O(|E| + |F|)

### 4. Mixed-Integer Solving

**Method:** Greedy rounding with continuous relaxation

**Process:**
1. Solve continuous relaxation: A x₀ = b (all variables real)
2. For i = 1 to k:
   - Find integer variable j with minimum rounding error: |x_j - round(x_j)|
   - Round x_j to nearest integer and fix
   - Update continuous variables (using direct resolution or Gauss-Seidel)
3. Return solution

**Why greedy works:** 
- Rounding the "easiest" variable first minimizes local perturbation
- Updating continuous variables maintains optimality for fixed integers
- Not globally optimal, but provides good approximation

**Integration:** Uses GreedyMIPSolver module from ../greedy_mip_solver/

**Complexity:** O(k × n³) for direct resolution, O(k × nnz × iter) for Gauss-Seidel

### 5. Solution Extraction

**Process:**
- Map solution vector back to θ and p arrays
- Set p_e = 0 for tree edges
- Verify constraints are satisfied

**Output:** Angle field θ, period jump field p

### 6. Singularity Computation

**Definition:** Vertex v is singular if its cross field index I(v) ≠ 0.

**Index Formula:**
```
I(v) = I₀(v) + Σ (p_ij / 4)
```

where:
- I₀(v) = (A_d(v) + Σ κ_ij) / (2π) is the geometric index
- A_d(v) is the angle defect: 2π - (sum of incident face angles)
- Sum is over edges in star(v)

**Physical Meaning:**
- I(v) = +1/4: Singularity with valence 3 in quad mesh
- I(v) = -1/4: Singularity with valence 5 in quad mesh
- I(v) = 0: Regular vertex (valence 4)

**Gauss-Bonnet:** Σ I(v) = χ(M) / 4, where χ is Euler characteristic

**Process:**
- For each vertex v:
  - Compute angle defect A_d(v)
  - Sum κ_ij and p_ij over incident edges
  - Calculate I(v)
  - If |I(v)| > tolerance, mark as singular

**Output:** List of (vertex_index, singularity_index) pairs

## Convergence and Optimality

### Convergence Guarantee

The algorithm always terminates because:
1. Finite number of integer variables to round
2. Each rounding step fixes one variable
3. Continuous updates converge (under mild conditions)

### Optimality

**Not globally optimal** because:
- Greedy choice of rounding order is heuristic
- No backtracking if poor choices are made

**Typically good in practice** because:
- Starting from continuous optimum provides good initialization
- Greedy selection minimizes local perturbation
- Continuous updates maintain optimality for fixed integers

**Alternative:** Could use branch-and-bound for global optimum, but much slower.

## Computational Complexity

| Step | Complexity | Dominant Factor |
|------|------------|-----------------|
| Topology | O(\|E\| + \|F\|) | Edge enumeration |
| Dijkstra | O(\|F\| log \|F\|) | Priority queue |
| Assembly | O(\|E\| + \|F\|) | Matrix construction |
| MIP Solve | O(k × n³) or O(k × nnz × iter) | System resolution |
| Singularities | O(\|V\| × \|E\|) | Index computation |

**Overall:** O(k × n³) for dense, O(k × nnz × iter) for sparse

For typical meshes:
- |F| ≈ 2|V| (triangular mesh)
- |E| ≈ 3|V| (triangular mesh)
- k ≈ |E| - |F| + 1 (free edges)

## References

1. Bommes, D., Zimmer, H., & Kobbelt, L. (2009). Mixed-integer quadrangulation. ACM TOG.
2. Ray, N., Vallet, B., Li, W. C., & Lévy, B. (2008). N-symmetry direction field design. ACM TOG.
3. Palacios, J., & Zhang, E. (2007). Rotational symmetry field design on surfaces. ACM TOG.
