# Implementation Comparison: Current vs Reference (MIQ_old)

## Critical Differences Found

### 1. **Matrix Assembly - Energy Derivatives**

#### Reference Implementation (`!!!MIQ_old/code/CrossField.jl`)
```julia
# Row i (for ∂E/∂θ_i)
if row_i > 0
    add!(row_i, row_i, 2.0)           # Coefficient for θ_i
    if row_j>0
        add!(row_i, row_j, -2.0)      # Coefficient for θ_j
    end
    if row_p>0
        add!(row_i, row_p, pi)        # Coefficient for p_ij
    end
    b[row_i] += 2.0 * const_RHS
end

# Row j (for ∂E/∂θ_j)
if row_j > 0
    if row_i>0; add!(row_j, row_i, -2.0); end
    add!(row_j, row_j, 2.0)
    if row_p>0; add!(row_j, row_p, -pi); end
    b[row_j] += -2.0 * const_RHS
end

# Row p (for ∂E/∂p_ij)
if row_p > 0
    if row_i>0; add!(row_p, row_i, pi); end
    if row_j>0; add!(row_p, row_j, -pi); end
    add!(row_p, row_p, pi^2/2.0)
    b[row_p] += pi * const_RHS
end
```

**Key Pattern**: For each edge (i,j):
- Creates **3 equation rows** (one for each of: θ_i, θ_j, p_ij)
- Each row represents the derivative of E with respect to that variable
- The derivatives are **accumulated** across all edges

#### Current Implementation (`frame_field/src/cross_field_energy.jl`)

**`assemble_system_matrix` version**:
```julia
# One row per variable (mixing energy derivatives with constraints)
# For p_ij variables:
for (var_idx, edge) in var_to_p
    row += 1
    if p_is_fixed[var_idx]
        # Constraint: p_ij = 0
        push!(vals, 1.0)
        b[row] = 0.0
    else
        # ∂E/∂p_ij equation
        push!(vals, π^2 / 2)   # p_ij coefficient
        push!(vals, π)          # θ_i coefficient
        push!(vals, -π)         # θ_j coefficient
    end
end

# For θ_k variables:
for (var_idx, face_k) in var_to_theta
    row += 1
    if theta_is_fixed[var_idx]
        # Constraint: θ_k = value
        push!(vals, 1.0)
        b[row] = constrained_angles[idx]
    else
        # ∂E/∂θ_k equation (summed over all neighbors)
        push!(vals, 2.0 * degree)  # Diagonal
        # ... neighbor contributions ...
    end
end
```

**`assemble_free_system` version**: Similar but only free variables.

### 2. **The Fundamental Problem**

#### Reference Approach:
- **Accumulates** derivatives across edges
- For a face with degree d, its θ variable appears in **d equations** (one per incident edge)
- Each edge contributes to **3 equations** (∂E/∂θ_i, ∂E/∂θ_j, ∂E/∂p_ij)
- Creates a **rectangular-ish system** where number of equations ≈ 3 × n_edges
- Each face variable accumulates contributions from all its neighbors

#### Current Approach:
- **One equation per variable**
- Tries to sum all contributions into a single equation per face
- For ∂E/∂θ_k, manually computes the sum: Σ_neighbors 2(θ_k + κ_kj + (π/2)p_kj - θ_j)
- Creates exactly **n_vars equations for n_vars unknowns**

### 3. **Why Current Implementation is Rank-Deficient**

The current implementation creates a square system, but:

1. **The energy formulation is inherently redundant**: 
   - E = Σ_edges (θ_i + κ_ij + (π/2)p_ij - θ_j)²
   - Taking derivatives: ∂E/∂θ_k involves summing over all edges incident to k
   - These summations create **linear dependencies** between equations

2. **Missing the constraint structure**:
   - The reference uses fixed edges (p_ij = 0) to break the rotational symmetry
   - But it also leverages the **over-determined** nature of the system
   - By having more equations than variables, the system can be properly solved via least squares

3. **Numerical issues**:
   - Current system is nearly singular (rank-deficient by 2-3)
   - Condition number ~1e17-1e18
   - The Gauss-Seidel relaxation diverges because the system has near-zero eigenvalues

### 4. **Solver Differences**

#### Reference:
```julia
function solve_greedy(A, b, p_var_map)
    F = cholesky(Hermitian(A))  # Assumes A is SPD
    x = F \ b
    
    # Greedy rounding...
    while true
        # Find best integer to round
        # Fix it
        # Relax with local_gauss_seidel!
    end
end
```

**Critical assumption**: `A` is **symmetric positive definite (SPD)**
- This requires the matrix to be well-conditioned
- The system must be properly constrained

#### Current:
```julia
function greedy_mip_solver(A, b, n_integer)
    # First n_integer variables are integers
    # Tries to solve with A\b (no SPD assumption)
    # Rounds and relaxes
end
```

## Root Cause Summary

**The current implementation fails because**:

1. **Wrong system structure**: Creates a square n×n system instead of the over-determined system the reference uses

2. **Missing accumulation**: Doesn't properly accumulate derivatives across edges - tries to compute sums in one equation per variable instead of building up contributions edge-by-edge

3. **Insufficient constraints**: Even with fixed edges, the system remains under-constrained because it doesn't have enough equations

4. **Violates SPD requirement**: The reference solver assumes A is SPD (from cholesky), but current system is rank-deficient

## Solution Path

To fix the current implementation:

1. **Adopt the reference's matrix assembly**:
   - Build system edge-by-edge, accumulating contributions
   - Create multiple equation rows per variable (one per incident edge)
   - Result will be over-determined (more equations than variables)

2. **Use proper solver**:
   - For the continuous part: use least squares or Cholesky on A^T A
   - The reference uses `cholesky(Hermitian(A))` - this assumes A is already SPD
   - Likely the reference's A is actually A^T A or already accumulated properly

3. **Alternative**: Keep current structure but add regularization:
   - Add small diagonal terms to break near-singularity
   - This is a hack and may not produce correct results

**Recommendation**: Rewrite `assemble_free_system` to match the reference's accumulation pattern.
