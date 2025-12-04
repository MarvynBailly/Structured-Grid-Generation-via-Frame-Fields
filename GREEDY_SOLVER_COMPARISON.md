# Greedy Mixed-Integer Solver Comparison

## Overview

This document compares the greedy mixed-integer solver implementations in:
1. **MIQ Code** (`MIQ/code/CrossField.jl`)
2. **Our Implementation** (`greedy_mip_solver/src/GreedyMIPSolver.jl`)

## Implementation Comparison

### Common Algorithm Structure

Both implementations follow the same high-level greedy algorithm:

1. **Continuous Relaxation**: Solve the full system Ax = b with all variables as reals
2. **Greedy Selection**: Find the integer variable with minimum rounding error (closest to an integer)
3. **Rounding**: Fix that variable to its nearest integer value
4. **Local Update**: Update continuous variables to accommodate the newly fixed integer
5. **Repeat**: Continue until all k integer variables are rounded

---

## Key Differences

### 1. **System Resolution After Rounding**

**MIQ Implementation (`solve_greedy` and `solve_greedy_global`):**
- Uses **Gauss-Seidel iterative method** for local updates
- After rounding a variable, adds affected neighbors to a queue
- Iteratively updates variables until residuals converge
- Two variants:
  - `solve_greedy`: Only adds direct neighbors of fixed variable
  - `solve_greedy_global`: Adds ALL continuous variables to queue for more robust updates

**Our Implementation:**
- Originally used Gauss-Seidel (similar to MIQ)
- **Changed to direct system resolution** after encountering numerical stability issues
- Builds a **reduced system** with only continuous variables:
  ```julia
  A_reduced = A[continuous_indices, continuous_indices]
  b_reduced = b[continuous_indices] - A[continuous_indices, rounded_indices] * x[rounded_indices]
  x[continuous_indices] = A_reduced \ b_reduced
  ```
- Solves the reduced system directly using backslash operator

### 2. **Queue Management**

**MIQ Implementation:**
```julia
# Add all affected neighbors to queue
for idx in nzrange(A, k)
    neighbor = rows[idx]
    if !fixed_mask[neighbor] && !in_queue[neighbor]
        enqueue!(q, neighbor)
        in_queue[neighbor] = true
    end
end
```

**Our Implementation (original Gauss-Seidel version):**
```julia
# Used PriorityQueue with residual-based priority
queue = PriorityQueue{Int, Float64}()
queue[m] = -abs(r_m)  # Negative for max-heap (process largest residuals first)
```
- Priority-based: processes variables with largest residuals first
- Added oscillation prevention: limit updates per variable to 50

### 3. **Matrix Handling**

**MIQ Implementation:**
- Assumes matrix is Symmetric Positive Definite (SPD)
- Uses Cholesky factorization in `solve_greedy`:
  ```julia
  F = cholesky(Hermitian(A))
  x = F \ b
  ```
- `solve_greedy_global` uses direct solve without factorization

**Our Implementation:**
- More general: works with any square matrix
- Uses direct solve with backslash operator
- No special assumptions about matrix structure

### 4. **Convergence Checking**

**MIQ Implementation:**
```julia
if abs(r_k) > tol
    x[k] += r_k / diag_A[k]
    # ... add neighbors to queue
end
```
- Checks residual against tolerance before updating
- Default tolerance: `tol=1e-5`

**Our Implementation:**
- Direct resolution doesn't require residual checking
- Solves reduced system to machine precision in one step
- More robust but potentially more expensive per iteration

### 5. **Sparse Matrix Optimization**

**MIQ Implementation:**
- Explicitly handles sparse matrices efficiently:
  ```julia
  rows = rowvals(A)
  vals = nonzeros(A)
  for idx in nzrange(A, k)
      Ax_k += vals[idx] * x[rows[idx]]
  end
  ```
- Manually computes matrix-vector products using sparse structure
- Avoids dense operations

**Our Implementation:**
- Relies on Julia's native sparse matrix operations
- Less manually optimized for sparsity
- Builds reduced dense/sparse systems depending on input

---

## Algorithmic Trade-offs

### MIQ's Gauss-Seidel Approach

**Advantages:**
- ✅ Sparse-aware: Only touches non-zero entries
- ✅ Memory efficient: No system reconstruction
- ✅ Incremental: Updates only what changed
- ✅ Scalable: O(nnz × iterations) complexity

**Disadvantages:**
- ⚠️ Convergence issues: May not converge for ill-conditioned systems
- ⚠️ Numerical instability: Can produce Inf/NaN values
- ⚠️ Iteration count unpredictable: May require many iterations
- ⚠️ Quality depends on system structure: Works best for diagonally dominant matrices

### Our Direct Resolution Approach

**Advantages:**
- ✅ Numerically stable: Direct solve is more robust
- ✅ Predictable: One solve per rounding step
- ✅ Guaranteed optimal: Solves reduced system exactly
- ✅ Simpler code: Fewer edge cases to handle

**Disadvantages:**
- ⚠️ More expensive: Full factorization at each step
- ⚠️ Less sparse-aware: May densify during factorization
- ⚠️ Memory intensive: Builds multiple reduced systems
- ⚠️ Complexity: O(n³) per rounding step in worst case

---

## Why We Switched to Direct Resolution

During testing, the Gauss-Seidel approach encountered several issues:

1. **Infinite Values**: Variables became Inf during updates
   ```
   Test 1: simple 2x2 system - Failed
     Expected: x ≈ [2.0, 1.0]
     Got: x = [Inf, Inf]
   ```

2. **Non-convergence**: Residuals remained high even after many iterations
   ```
   Test 2: 4x4 with 2 integers - Failed
     Expected residual < 1e-4
     Got residual: 2.47
   ```

3. **Oscillation**: Some variables oscillated without converging

The direct resolution approach fixed all these issues and provides exact solutions for the reduced continuous system after each integer rounding.

---

## Performance Characteristics

### Time Complexity

| Approach | Per-Round Complexity | Total Complexity |
|----------|---------------------|------------------|
| MIQ Gauss-Seidel | O(nnz × iterations) | O(k × nnz × iterations) |
| Our Direct Solve | O(n³) worst case | O(k × n³) worst case |

For **sparse systems** with many variables, Gauss-Seidel is typically faster.
For **dense systems** or systems requiring high accuracy, direct solve is better.

### Memory Usage

| Approach | Memory |
|----------|--------|
| MIQ Gauss-Seidel | O(nnz) - only stores sparse matrix and queue |
| Our Direct Solve | O(n²) - builds reduced systems |

---

## Accuracy Comparison

### MIQ Implementation:
- **Approximate**: Gauss-Seidel stops when residuals < tolerance
- **Residuals**: Typically 1e-5 to 1e-6 for well-conditioned systems
- **Trade-off**: Faster but less accurate

### Our Implementation:
- **Exact** (within machine precision): Direct solve finds optimal continuous solution
- **Residuals**: Machine epsilon (~1e-15 for Float64)
- **Trade-off**: Slower but more accurate

However, both are **approximate for the original MIP problem** because:
- Greedy selection is heuristic (not globally optimal)
- Order of rounding affects final solution
- No backtracking or branch-and-bound

---

## Code Structure Comparison

### MIQ (`CrossField.jl`)

```julia
function solve_greedy_global(A, b, p_var_map)
    # 1. Initial solve
    x = A \ b
    
    # 2. Rounding loop
    while true
        # Find best candidate
        best_idx = find_min_rounding_error(int_indices, x, fixed_mask)
        if best_idx == -1; break; end
        
        # Fix variable
        x[best_idx] = round(x[best_idx])
        fixed_mask[best_idx] = true
        
        # Relax with Gauss-Seidel
        q = Queue{Int}()
        # Add all continuous variables
        for i in all_vars
            if !fixed_mask[i]
                enqueue!(q, i)
            end
        end
        run_gauss_seidel_queue!(A, b, x, fixed_mask, q, in_queue)
    end
    return x
end
```

### Our Implementation (`GreedyMIPSolver.jl`)

```julia
function solve_greedy_mip(A, b, k; τ=1e-6, max_iterations=10000, verbose=false)
    # 1. Initial solve
    x = A \ b
    
    # 2. Rounding loop
    for round_iter in 1:k
        # Find best candidate
        best_idx = find_min_rounding_error(1:k, x, rounded)
        
        # Fix variable
        x[best_idx] = round(x[best_idx])
        push!(rounded, best_idx)
        
        # Direct resolution for continuous variables
        continuous_indices = [i for i in 1:n if i ∉ rounded]
        if !isempty(continuous_indices)
            A_reduced = A[continuous_indices, continuous_indices]
            b_reduced = b[continuous_indices] - A[continuous_indices, rounded] * x[rounded]
            x[continuous_indices] = A_reduced \ b_reduced
        end
    end
    return GreedyMIPResult(x, converged, total_iterations, rounding_order)
end
```

---

## Recommendations

### Use MIQ's Gauss-Seidel Approach When:
- ✅ Working with large sparse systems (> 10,000 variables)
- ✅ System is well-conditioned (diagonally dominant)
- ✅ Speed is more important than accuracy
- ✅ Moderate tolerance (1e-5 to 1e-6) is acceptable

### Use Our Direct Resolution Approach When:
- ✅ Numerical stability is critical
- ✅ High accuracy required (residuals < 1e-10)
- ✅ System is small to medium size (< 1,000 variables)
- ✅ System may be ill-conditioned
- ✅ Debugging/prototyping (simpler to understand and validate)

---

## Potential Improvements

### For MIQ's Implementation:
1. Add **preconditioning** to improve Gauss-Seidel convergence
2. Implement **successive over-relaxation (SOR)** for faster convergence
3. Add **residual monitoring** to detect non-convergence early
4. Use **Cholesky factorization** consistently for SPD systems

### For Our Implementation:
1. Add **sparse matrix optimization** to avoid densification
2. Implement **iterative refinement** after direct solve
3. Add **hybrid approach**: Use Gauss-Seidel for large sparse systems, direct solve for small systems
4. Cache factorizations when possible (e.g., when only RHS changes)

---

## Conclusion

Both implementations follow the same greedy algorithm but differ in how they update continuous variables after rounding:

- **MIQ uses iterative Gauss-Seidel**: Fast, memory-efficient, but can have convergence issues
- **Our implementation uses direct resolution**: Slower, more memory-intensive, but numerically stable

The choice depends on problem size, sparsity, conditioning, and accuracy requirements. For the frame field application with moderate-sized meshes (~1000 faces), our direct approach provides reliable results. For very large meshes, MIQ's sparse Gauss-Seidel approach would be more scalable.

A **hybrid approach** combining both methods could be optimal:
- Use direct solve for the first few roundings (when system is still well-conditioned)
- Switch to Gauss-Seidel for later roundings (when system becomes sparse and structured)
