"""
    GreedyMIPSolver.jl

Greedy Mixed-Integer Programming Solver implementation based on Bommes et al. (2009).

This module provides a fast approximation to solving mixed-integer linear systems
of the form Ax = b, where the first k variables are constrained to be integers.

The greedy approach rounds integer variables one at a time (choosing the variable
with minimum rounding error) and updates the continuous variables using local
Gauss-Seidel iterations.
"""
module GreedyMIPSolver

export solve_greedy_mip, GreedyMIPResult

using LinearAlgebra
using SparseArrays
using DataStructures  # For PriorityQueue

"""
    GreedyMIPResult

Container for the results of the greedy MIP solver.

# Fields
- `x::Vector{Float64}`: Solution vector
- `converged::Bool`: Whether the solver converged
- `iterations::Int`: Total number of Gauss-Seidel iterations performed
- `rounding_order::Vector{Int}`: Order in which integer variables were rounded
"""
struct GreedyMIPResult
    x::Vector{Float64}
    converged::Bool
    iterations::Int
    rounding_order::Vector{Int}
end

"""
    solve_greedy_mip(A, b, k; τ=1e-6, max_iterations=10000, verbose=false)

Solve the mixed-integer linear system Ax = b where the first k variables
must be integers, using a greedy rounding strategy.

# Arguments
- `A::AbstractMatrix`: System matrix (n × n)
- `b::AbstractVector`: Right-hand side vector (length n)
- `k::Int`: Number of integer variables (must satisfy 1 ≤ k ≤ n)
- `τ::Float64=1e-6`: Convergence threshold for residuals
- `max_iterations::Int=10000`: Maximum Gauss-Seidel iterations
- `verbose::Bool=false`: Print progress information

# Returns
- `GreedyMIPResult`: Object containing solution and convergence information

# Algorithm
1. Solve continuous relaxation: Ax = b (all variables real-valued)
2. Greedily round integer variables:
   - Select variable with minimum rounding error
   - Round to nearest integer and fix
   - Update affected continuous variables via Gauss-Seidel
3. Repeat until all k integer variables are rounded

# Example
```julia
A = [2.0 1.0; 1.0 3.0]
b = [5.0, 7.0]
result = solve_greedy_mip(A, b, 1)  # x[1] must be integer
println("Solution: ", result.x)
println("Converged: ", result.converged)
```
"""
function solve_greedy_mip(
    A::AbstractMatrix,
    b::AbstractVector,
    k::Int;
    τ::Float64=1e-6,
    max_iterations::Int=10000,
    verbose::Bool=false
)
    n = length(b)
    
    # Input validation
    @assert size(A) == (n, n) "Matrix A must be square"
    @assert 1 ≤ k ≤ n "Number of integer variables k must satisfy 1 ≤ k ≤ n"
    @assert all(abs.(diag(A)) .> eps()) "Matrix A must have non-zero diagonal"
    
    # Step 1: Solve continuous relaxation
    verbose && println("Solving continuous relaxation...")
    x = A \ b
    
    rounded = Set{Int}()
    rounding_order = Int[]
    total_iterations = 0
    converged = true
    
    # Step 2: Greedy rounding loop
    for round_iter in 1:k
        verbose && println("\nRounding iteration $round_iter/$k")
        
        # Find variable with minimum rounding error among unrounded integer variables
        min_error = Inf
        best_idx = 0
        
        for j in 1:k
            if j ∉ rounded
                error = abs(x[j] - round(x[j]))
                if error < min_error
                    min_error = error
                    best_idx = j
                end
            end
        end
        
        # Safety check
        if best_idx == 0
            @error "No unrounded variable found (this should not happen)"
            break
        end
        
        # Round and fix the selected variable
        rounded_value = round(x[best_idx])
        verbose && println("  Rounding x[$best_idx]: $(x[best_idx]) → $rounded_value (error: $min_error)")
        x[best_idx] = rounded_value
        push!(rounded, best_idx)
        push!(rounding_order, best_idx)
        
        # Step 3: Resolve for continuous variables with fixed integers
        # Build reduced system for continuous variables
        continuous_indices = [i for i in 1:n if i ∉ rounded]
        
        if !isempty(continuous_indices)
            # Build reduced system
            A_reduced = A[continuous_indices, continuous_indices]
            b_reduced = b[continuous_indices] - A[continuous_indices, collect(rounded)] * x[collect(rounded)]
            
            # Solve reduced system
            try
                x_continuous = A_reduced \ b_reduced
                x[continuous_indices] = x_continuous
                total_iterations += 1
            catch e
                @warn "Failed to solve reduced system after rounding x[$best_idx]: $e"
                converged = false
            end
        end
    end
    
    verbose && println("\nTotal Gauss-Seidel iterations: $total_iterations")
    
    return GreedyMIPResult(x, converged, total_iterations, rounding_order)
end

"""
    gauss_seidel_update!(A, b, x, fixed_idx, rounded, τ, max_iter, verbose)

Perform local Gauss-Seidel updates after fixing a variable.

Uses a priority queue to process variables that depend on the fixed variable,
updating them iteratively until convergence or maximum iterations reached.

# Arguments
- `A::AbstractMatrix`: System matrix
- `b::AbstractVector`: Right-hand side
- `x::AbstractVector`: Current solution (modified in-place)
- `fixed_idx::Int`: Index of the variable just fixed
- `rounded::Set{Int}`: Set of already-rounded variable indices
- `τ::Float64`: Convergence threshold
- `max_iter::Int`: Maximum iterations
- `verbose::Bool`: Print debug info

# Returns
- `Int`: Number of iterations performed
"""
function gauss_seidel_update!(
    A::AbstractMatrix,
    b::AbstractVector,
    x::AbstractVector,
    fixed_idx::Int,
    rounded::Set{Int},
    τ::Float64,
    max_iter::Int,
    verbose::Bool
)
    n = length(b)
    
    # Initialize priority queue with variables that depend on fixed_idx
    # Priority is negative residual magnitude (larger residuals processed first)
    queue = PriorityQueue{Int, Float64}()
    visited = Set{Int}()
    
    # Add all variables m where A[m, fixed_idx] ≠ 0
    for m in 1:n
        if m ∉ rounded && abs(A[m, fixed_idx]) > eps()
            r_m = compute_residual(A, b, x, m)
            if abs(r_m) ≥ τ
                queue[m] = -abs(r_m)  # Negative for max-heap behavior
            end
        end
    end
    
    iter = 0
    processed_count = Dict{Int, Int}()  # Track how many times each variable was updated
    
    while !isempty(queue) && iter < max_iter
        iter += 1
        
        # Pop variable with largest residual
        m = dequeue!(queue)
        
        # Skip if we've updated this variable too many times (prevent oscillation)
        count_m = get(processed_count, m, 0)
        if count_m > 50
            continue
        end
        
        # Compute residual: r_m = b_m - Σ_n A[m,n] * x[n]
        r_m = compute_residual(A, b, x, m)
        
        if abs(r_m) ≥ τ && abs(A[m, m]) > eps()
            # Update variable: x[m] ← x[m] - r_m / A[m,m]
            old_x = x[m]
            x[m] -= r_m / A[m, m]
            
            # Check for numerical issues
            if !isfinite(x[m])
                x[m] = old_x  # Revert if update causes problems
                verbose && println("    Warning: Non-finite value detected for x[$m], reverting")
                continue
            end
            
            processed_count[m] = count_m + 1
            
            # Add dependent variables to queue (only continuous variables)
            for l in 1:n
                if l ∉ rounded && abs(A[l, m]) > eps()
                    r_l = compute_residual(A, b, x, l)
                    if abs(r_l) ≥ τ
                        queue[l] = -abs(r_l)
                    end
                end
            end
        end
    end
    
    verbose && iter > 0 && println("    Gauss-Seidel iterations: $iter")
    
    return iter
end

"""
    compute_residual(A, b, x, row)

Compute the residual for a single row: r = b[row] - (A[row, :] ⋅ x)
"""
function compute_residual(A::AbstractMatrix, b::AbstractVector, x::AbstractVector, row::Int)
    return b[row] - dot(A[row, :], x)
end

end # module
