# Fix Summary: Frame Field Energy System

## Problem
The frame field energy system was severely ill-conditioned:
- Rank-deficient (33-44 out of 47-49 needed)
- Condition number ~1e17-1e18
- Gauss-Seidel diverged during greedy rounding

## Root Causes Identified

### 1. Insufficient Edge Constraints
**Issue**: Only fixing a small subset (~12) of edges from the spanning forest  
**Reference approach**: Fixes ALL spanning forest edges (~25 for 26 faces)  
**Impact**: Without enough fixed edges, the system has rotational freedom (can add constant to all angles)

### 2. Incorrect Use of fixed_edges Parameter
**Issue**: `assemble_free_system` was using `fixed_edges_per_face` values (subset) instead of `potential_fixed_edges` (full spanning forest)  
**Fix**: Changed to use `potential_fixed_edges` directly as the set of fixed edges

### 3. Gauss-Seidel Implementation Issues
**Issue 1**: Max iterations too low (1000 vs reference's 500000-1000000)  
**Issue 2**: Residual computation had sign error in update step  
**Issue 3**: Tolerance too strict (1e-6 vs reference's 1e-5)

## Solutions Implemented

### Fix 1: Use Full Spanning Forest as Fixed Edges
```julia
# In assemble_free_system:
fixed_edges = potential_fixed_edges  # Use ALL spanning forest edges
```

### Fix 2: Edge-by-Edge Accumulation (Already Correct)
The matrix assembly already follows the reference pattern:
- Iterate over all edges
- Accumulate contributions to equation rows for θ_i, θ_j, and p_ij
- Creates symmetric positive definite system

### Fix 3: Corrected Gauss-Seidel
```julia
# Increased limits and relaxed tolerance:
tau::Float64=1e-5,          # Was 1e-6
max_iterations::Int=1000000, # Was 1000

# Fixed residual computation:
Ax_m = Σ_n A_mn * x_n      # Full matrix-vector product
residual = b_m - Ax_m       # Residual
x_m ← x_m + residual / A_mm  # Update (was x_m - residual)
```

## Results

### Before Fix
- System size: 49 × 49 variables
- Rank: 33/49 (16 dimensional null space)
- Condition number: 2.85e17
- 15 near-zero eigenvalues
- Greedy solver: **FAILED** (Gauss-Seidel diverged)

### After Fix
- System size: 33 × 33 variables (fewer free edges)
- Rank: 33/33 (**full rank!**)
- Condition number: **~380** (excellent!)
- All eigenvalues positive
- Greedy solver: **SUCCESS** (2 iterations per rounding step)

### Verification
- All period jumps are integers ✓
- Constrained face angles match exactly ✓
- Face angles nearly uniform (~20.35°) ✓
- Energy properly minimized ✓
- Cross field visualization generated ✓

## Key Insight from Reference Code

The reference implementation (`!!!MIQ_old/code/CrossField.jl`) has a critical step at lines 743-753:

```julia
# let's pick only one edge per triangle to be fixed
for f_idx in 1:length(faces(mesh_info))
    # find the first edge that is not fixed yet and fix it
    for edge in face_edges
        e_idx = topo.edge_map[edge]
        if !fixed_edges[e_idx]
            fixed_edges[e_idx] = true
            break
        end
    end
end
```

This adds ADDITIONAL fixed edges beyond the spanning forest, fixing at least one edge per face. This heavily constrains the system and breaks any remaining rotational symmetry.

Our approach of fixing all spanning forest edges achieves a similar effect - approximately n-1 fixed edges for n faces provides enough constraint to make the system well-conditioned.

## Remaining Work
The `fix_suitable_edges` function and `fixed_edges_per_face` parameter are now somewhat redundant since we're fixing the entire spanning forest. Consider:
- Simplifying the API to remove the intermediate selection step
- Or: Adding even MORE fixed edges (like the reference does) for additional stability
