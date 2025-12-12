"""
    greedy_solver.jl

Greedy mixed-integer solver for frame field optimization.

Implements the algorithm described by Bommes et al. for solving mixed-integer
quadratic programs. The solver greedily rounds integer variables one at a time,
using local Gauss-Seidel updates to maintain consistency with continuous variables.
"""

using LinearAlgebra
using SparseArrays
using DataStructures

"""
    greedy_mip_solver(
        A::SparseMatrixCSC,
        b::Vector{Float64},
        n_integer::Int;
        tau::Float64=1e-6,
        max_iterations::Int=1000,
        verbose::Bool=false
    ) -> Vector{Float64}

Solve a mixed-integer quadratic program using greedy rounding.

The problem is: minimize (1/2)x'Ax - b'x
subject to x[1:n_integer] ∈ ℤ (integers), x[n_integer+1:end] ∈ ℝ (reals).

# Algorithm
1. Solve continuous relaxation: Ax = b
2. For each integer variable (in order of smallest rounding error):
   - Round to nearest integer
   - Update affected continuous variables using Gauss-Seidel iterations
   - Continue until local residuals are below threshold

# Arguments
- `A` - Sparse system matrix (n × n)
- `b` - Right-hand side vector (length n)
- `n_integer` - Number of integer variables (first n_integer variables)
- `tau` - Threshold for Gauss-Seidel convergence (default: 1e-6)
- `max_iterations` - Maximum Gauss-Seidel iterations per rounding step (default: 1000)
- `verbose` - Print progress information (default: false)

# Returns
- Solution vector x with x[1:n_integer] rounded to integers
"""
function greedy_mip_solver(
    A::SparseMatrixCSC,
    b::Vector{Float64},
    n_integer::Int;
    tau::Float64=1e-5,
    max_iterations::Int=1000000,
    verbose::Bool=false
)
    n = length(b)
    
    if n_integer > n
        error("n_integer ($n_integer) cannot exceed total variables ($n)")
    end
    
    verbose && println("=== Greedy Mixed-Integer Solver ===")
    verbose && println("Total variables: $n")
    verbose && println("Integer variables: $n_integer")
    verbose && println("Continuous variables: $(n - n_integer)")
    
    # Step 1: Solve continuous relaxation
    verbose && println("\n[Step 1] Solving continuous relaxation...")
    x = A \ b
    
    if verbose
        initial_energy = 0.5 * dot(x, A * x) - dot(b, x)
        println("  Initial energy: $initial_energy")
    end
    
    # Track which variables have been rounded
    rounded = Set{Int}()
    
    # Step 2: Greedily round integer variables
    verbose && println("\n[Step 2] Greedy rounding of integer variables...")
    
    for round_step in 1:n_integer
        # Find integer variable with minimum rounding error
        min_error = Inf
        best_idx = -1
        
        for j in 1:n_integer
            if j ∉ rounded
                error_val = abs(x[j] - round(x[j]))
                if error_val < min_error
                    min_error = error_val
                    best_idx = j
                end
            end
        end
        
        if best_idx == -1
            error("Failed to find unrounded integer variable")
        end
        
        # Round and fix this variable
        old_val = x[best_idx]
        new_val = round(x[best_idx])
        x[best_idx] = new_val
        push!(rounded, best_idx)
        
        if verbose && (round_step <= 5 || round_step % 10 == 0 || round_step == n_integer)
            println("  Round $round_step/$n_integer: Variable $best_idx: $old_val → $new_val (error: $min_error)")
        end
        
        # Step 3: Update affected variables using Gauss-Seidel
        # Find all variables that depend on the just-rounded variable
        # (i.e., rows where A[i, best_idx] ≠ 0)
        rows = rowvals(A)
        vals = nonzeros(A)
        
        # Get column best_idx to find affected rows
        affected_vars = Set{Int}()
        col_range = nzrange(A, best_idx)
        for i in col_range
            row_idx = rows[i]
            if row_idx ∉ rounded  # Only update non-rounded variables
                push!(affected_vars, row_idx)
            end
        end
        
        # Priority queue for Gauss-Seidel iterations
        Q = Queue{Int}()
        for var in affected_vars
            enqueue!(Q, var)
        end
        
        # Gauss-Seidel iterations
        iter_count = 0
        gs_updates = 0
        
        while !isempty(Q) && iter_count < max_iterations
            iter_count += 1
            m = dequeue!(Q)
            
            # Compute Ax_m = Σ_n A_mn * x_n
            row_range = nzrange(A, m)
            Ax_m = 0.0
            
            for idx in row_range
                col = rows[idx]
                val = vals[idx]
                Ax_m += val * x[col]
            end
            
            # Residual: r_m = b_m - Ax_m
            residual = b[m] - Ax_m
            
            # Check if update is needed
            if abs(residual) >= tau
                # Get diagonal element A[m,m]
                A_mm = A[m, m]
                if abs(A_mm) < 1e-12
                    continue  # Skip if diagonal is too small
                end
                
                # Update variable: x_m ← x_m + r_m / A_mm
                x[m] = x[m] + residual / A_mm
                gs_updates += 1
                
                # Add dependent variables to queue
                # Find all variables n where A[n,m] ≠ 0
                col_range_m = nzrange(A, m)
                for idx in col_range_m
                    n = rows[idx]
                    if n ∉ rounded && n != m
                        enqueue!(Q, n)
                    end
                end
            end
        end
        
        if verbose && (round_step <= 5 || round_step % 10 == 0 || round_step == n_integer)
            println("    Gauss-Seidel: $gs_updates updates in $iter_count iterations")
        end
        
        if iter_count >= max_iterations
            @warn "Gauss-Seidel did not converge in $max_iterations iterations at round step $round_step"
        end
    end
    
    # Verify integer constraints
    for i in 1:n_integer
        if abs(x[i] - round(x[i])) > 1e-9
            @warn "Variable $i is not properly rounded: $(x[i])"
        end
    end
    
    if verbose
        final_energy = 0.5 * dot(x, A * x) - dot(b, x)
        println("\n[Complete]")
        println("  Final energy: $final_energy")
        println("  Integer variables rounded: $n_integer")
    end
    
    return x
end

"""
    extract_integer_variables(x::Vector{Float64}, n_integer::Int) -> Vector{Int}

Extract and convert the first n_integer variables to integers.
"""
function extract_integer_variables(x::Vector{Float64}, n_integer::Int)
    return round.(Int, x[1:n_integer])
end

"""
    extract_continuous_variables(x::Vector{Float64}, n_integer::Int) -> Vector{Float64}

Extract the continuous variables (after the integer variables).
"""
function extract_continuous_variables(x::Vector{Float64}, n_integer::Int)
    return x[n_integer+1:end]
end
