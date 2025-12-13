using LinearAlgebra
using SparseArrays
using DataStructures

"""
    greedy_mixed_integer_solve(
        A::SparseMatrixCSC,
        b::Vector{Float64},
        k::Int;
        max_gauss_seidel_iters::Int=1000,
        tau::Float64=1e-6,
        fallback_to_full_solve::Bool=true,
        use_queue::Bool=true,
        verbose::Bool=false
    ) -> Vector{Float64}

Greedy mixed integer solver for linear systems where the first k variables must be integers.

Rather than rounding all integer variables at once, this method rounds one integer variable 
at a time, choosing the variable that causes the smallest absolute error when rounded. After 
each rounding, the continuous variables are updated using a local Gauss-Seidel update.

# Algorithm
1. Solve the continuous linear system Ax = b to get x⁰
2. For each step:
   - Find the integer variable xᵢ (i ≤ k) that causes smallest error when rounded
   - Round xᵢ to nearest integer and fix it
   - Update dependent variables using queue-based or priority-queue Gauss-Seidel
3. Continue until all k integer variables are fixed

# Update Methods
- `use_queue=true`: Simple queue-based Gauss-Seidel (CrossField.jl style, faster and more robust)
- `use_queue=false`: Priority queue-based Gauss-Seidel (original implementation)

# Arguments
- `A::SparseMatrixCSC` - Sparse system matrix (n × n)
- `b::Vector{Float64}` - Right-hand side vector (length n)
- `k::Int` - Number of integer variables (first k variables must be integers)
- `max_gauss_seidel_iters::Int` - Maximum Gauss-Seidel iterations (default: 1000)
- `tau::Float64` - Convergence threshold for residual (default: 1e-6)
- `fallback_to_full_solve::Bool` - Fall back to full solve if GS doesn't converge (default: true)
- `use_queue::Bool` - Use simple queue (true) or priority queue (false) for updates (default: true)
- `mod4::Bool` - Apply mod 4 constraint to keep period jumps in [0,3] (default: true)
- `callback::Union{Function,Nothing}` - Optional callback(x, step, k) called after each rounding step (default: nothing)
- `verbose::Bool` - Print progress information (default: false)

# Returns
- `x::Vector{Float64}` - Solution vector with first k entries as integers

# References
Bommes et al., "Mixed-integer quadrangulation", ACM TOG 2009
"""
function greedy_mixed_integer_solve(
    A::SparseMatrixCSC,
    b::Vector{Float64},
    k::Int;
    max_gauss_seidel_iters::Int=1000,
    tau::Float64=1e-6,
    fallback_to_full_solve::Bool=true,
    use_queue::Bool=true,
    mod4::Bool=true,
    callback::Union{Function,Nothing}=nothing,
    verbose::Bool=true
)
    n = length(b)
    @assert size(A) == (n, n) "Matrix A must be square"
    @assert k <= n "Number of integer variables k must be ≤ n"
    @assert k >= 0 "Number of integer variables k must be non-negative"
    
    if k == 0
        # No integer variables, just solve the linear system
        return qr(Matrix(A)) \ b
    end
    
    # Step 1: Solve continuous relaxation
    if verbose
        println("Solving continuous relaxation...")
    end
    
    # Use least-squares solve to handle potentially singular/rectangular systems
    # Julia's \ operator automatically chooses least-squares for rectangular systems
    # For square but rank-deficient systems, we need to be more careful
    if size(A, 1) == size(A, 2)
        # Square system - try direct solve, fall back to least-squares if singular
        try
            x = A \ b
        catch e
            if verbose
                println("  System is singular, using least-squares solve...")
            end
            # Use QR-based least-squares: minimize ||Ax - b||₂
            x = qr(Matrix(A)) \ b
        end
    else
        # Rectangular system - use least-squares
        x = A \ b
    end
    
    if verbose
        println("Initial solution computed")
        println("  First k values: ", x[1:min(k, 10)])
        residual_norm = norm(A * x - b)
        println("  Initial residual norm: $residual_norm")
    end
    
    # Track which variables are already fixed to integers
    fixed = falses(n)
    n_fixed = 0
    
    # Precompute row indices for efficient access
    rows = rowvals(A)
    vals = nonzeros(A)
    
    # Precompute diagonal for queue-based Gauss-Seidel
    diag_A = Vector{Float64}(undef, n)
    for i in 1:n
        diag_A[i] = A[i, i]
    end
    
    # Step 2: Greedy rounding loop
    while n_fixed < k
        # Find the integer variable (among unfixed i ≤ k) that causes smallest rounding error
        best_idx = -1
        best_error = Inf
        best_rounded_val = 0.0
        
        for i in 1:k
            if !fixed[i]
                if mod4
                    # Find the best value in {0, 1, 2, 3} that is closest to x[i]
                    # This respects the π/2 rotational symmetry of frame fields
                    min_error = Inf
                    best_candidate = 0.0
                    
                    for candidate in 0:3
                        # Check values offset by ±4 to handle wrapping
                        for offset in [-4, 0, 4]
                            test_val = Float64(candidate + offset)
                            error = abs(x[i] - test_val)
                            if error < min_error
                                min_error = error
                                best_candidate = Float64(candidate)
                            end
                        end
                    end
                    
                    rounded = best_candidate
                    error = min_error
                else
                    # Standard rounding to nearest integer
                    rounded = round(x[i])
                    error = abs(x[i] - rounded)
                end
                
                if error < best_error
                    best_error = error
                    best_idx = i
                    best_rounded_val = rounded
                end
            end
        end
        
        @assert best_idx > 0 "Should have found a variable to round"
        
        # Round the best variable
        old_val = x[best_idx]
        x[best_idx] = best_rounded_val
        fixed[best_idx] = true
        n_fixed += 1
        
        if verbose
            println("\nRounding step $n_fixed/$k:")
            println("  Variable x[$best_idx]: $old_val → $best_rounded_val (error: $best_error)")
        end
        
        # Step 3: Local update with fixed variables
        if use_queue
            # Use simple queue-based Gauss-Seidel (CrossField.jl style)
            n_iters = gauss_seidel_queue_update!(
                A, b, x, fixed, best_idx, diag_A,
                max_gauss_seidel_iters, tau
            )
            
            if verbose && n_iters > 100
                println("  Queue-based relaxation: $n_iters iterations")
            end
        else
            # Use priority queue-based Gauss-Seidel (original)
            converged = gauss_seidel_update!(
                A, b, x, best_idx, fixed,
                max_gauss_seidel_iters, tau, verbose
            )
            
            if !converged && fallback_to_full_solve
                if verbose
                    println("  Gauss-Seidel did not converge, falling back to full solve")
                end
                x = solve_with_fixed_variables(A, b, x, fixed)
            end
        end
        
        # Call callback if provided
        if callback !== nothing
            callback(copy(x), n_fixed, k)
        end
    end
    
    if verbose
        println("\nGreedy rounding complete!")
        println("  All $k integer variables fixed")
        println("  Final integer values: ", Int.(x[1:k]))
    end
    
    return x
end


"""
    gauss_seidel_queue_update!(
        A, b, x, fixed, start_node, diag_A,
        max_iters, tau
    ) -> Int

Queue-based local Gauss-Seidel update (CrossField.jl style).

After fixing a variable, propagates updates through neighbors using a simple FIFO queue.
More robust and often faster than priority queue approach.

Returns the number of iterations performed.
"""
function gauss_seidel_queue_update!(
    A::SparseMatrixCSC,
    b::Vector{Float64},
    x::Vector{Float64},
    fixed::BitVector,
    start_node::Int,
    diag_A::Vector{Float64},
    max_iters::Int,
    tau::Float64
)
    n = length(b)
    rows = rowvals(A)
    vals = nonzeros(A)
    
    # Simple FIFO queue
    q = Queue{Int}()
    in_queue = Set{Int}()
    
    # Add neighbors of start_node to queue
    for idx in nzrange(A, start_node)
        neighbor = rows[idx]
        if !fixed[neighbor] && !(neighbor in in_queue)
            enqueue!(q, neighbor)
            push!(in_queue, neighbor)
        end
    end
    
    iter = 0
    while !isempty(q) && iter < max_iters
        iter += 1
        k = dequeue!(q)
        delete!(in_queue, k)
        
        if !fixed[k]
            # Compute residual: r_k = b[k] - Σ A_kj * x_j
            Ax_k = 0.0
            for idx in nzrange(A, k)
                Ax_k += vals[idx] * x[rows[idx]]
            end
            r_k = b[k] - Ax_k
            
            if abs(r_k) > tau
                # Update variable
                x[k] += r_k / diag_A[k]
                
                # Add neighbors to queue
                for idx in nzrange(A, k)
                    neighbor = rows[idx]
                    if !fixed[neighbor] && !(neighbor in in_queue)
                        enqueue!(q, neighbor)
                        push!(in_queue, neighbor)
                    end
                end
            end
        end
    end
    
    return iter
end


"""
    gauss_seidel_update!(
        A, b, x, changed_idx, fixed,
        max_iters, tau, verbose
    ) -> Bool

Perform local Gauss-Seidel updates after fixing variable `changed_idx`.

Uses a priority queue to propagate updates to dependent variables.
Returns true if converged, false if max iterations reached.
"""
function gauss_seidel_update!(
    A::SparseMatrixCSC,
    b::Vector{Float64},
    x::Vector{Float64},
    changed_idx::Int,
    fixed::BitVector,
    max_iters::Int,
    tau::Float64,
    verbose::Bool
)
    n = length(b)
    rows = rowvals(A)
    vals = nonzeros(A)
    
    # Priority queue: higher priority = larger absolute residual
    # Store tuples of (-abs(residual), var_index) so larger residuals have higher priority
    queue = PriorityQueue{Int, Float64}()
    in_queue = falses(n)
    
    # Add all variables that depend on changed_idx to the queue
    # These are variables k where A[k, changed_idx] ≠ 0
    for i in nzrange(A, changed_idx)
        row_idx = rows[i]
        if !fixed[row_idx]
            # Compute initial residual for this variable
            r = compute_residual(A, b, x, row_idx)
            println("Residual for variable $row_idx after changing $changed_idx: $r")
            if abs(r) >= tau
                queue[row_idx] = -abs(r)  # Negative for max-heap behavior
                in_queue[row_idx] = true
            end
        end
    end
    
    iter = 0
    n_updates = 0
    
    while !isempty(queue) && iter < max_iters
        iter += 1
        
        # Pop variable with largest residual
        var_idx = dequeue!(queue)
        in_queue[var_idx] = false
        
        # Compute current residual: r_k = b_k - Σⱼ A_kj x_j
        r = compute_residual(A, b, x, var_idx)
        
        # Check if update is needed
        if abs(r) >= tau
            # Get diagonal element
            A_kk = 0.0
            for i in nzrange(A, var_idx)
                if rows[i] == var_idx
                    A_kk = vals[i]
                    break
                end
            end
            
            @assert abs(A_kk) > 1e-12 "Diagonal element A[$var_idx,$var_idx] is too small: $A_kk"
            
            # Update: x_k ← x_k - r_k / A_kk
            old_val = x[var_idx]
            x[var_idx] = x[var_idx] - r / A_kk
            n_updates += 1
            
            # Add dependent variables to queue
            for i in nzrange(A, var_idx)
                dep_idx = rows[i]
                if !fixed[dep_idx] && dep_idx != var_idx && !in_queue[dep_idx]
                    r_dep = compute_residual(A, b, x, dep_idx)
                    if abs(r_dep) >= tau
                        queue[dep_idx] = -abs(r_dep)
                        in_queue[dep_idx] = true
                    end
                end
            end
        end
    end
    
    converged = isempty(queue)
    
    if verbose
        if converged
            println("  Gauss-Seidel converged in $iter iterations ($n_updates updates)")
        else
            println("  Gauss-Seidel reached max iterations ($iter), $n_updates updates performed")
        end
    end
    
    return converged
end


"""
    compute_residual(A, b, x, row_idx) -> Float64

Compute the residual for row `row_idx`: r = b[row_idx] - (A[row_idx,:] * x)
"""
function compute_residual(
    A::SparseMatrixCSC,
    b::Vector{Float64},
    x::Vector{Float64},
    row_idx::Int
)
    # Compute: r = b[row_idx] - Σⱼ A[row_idx, j] * x[j]
    r = b[row_idx]
    
    rows = rowvals(A)
    vals = nonzeros(A)
    
    # A is stored in CSC format (column-major), so we need to iterate through columns
    # to find entries in a specific row
    n = size(A, 2)
    for col in 1:n
        for i in nzrange(A, col)
            if rows[i] == row_idx
                r -= vals[i] * x[col]
                break
            end
        end
    end
    
    return r
end


"""
    solve_with_fixed_variables(A, b, x, fixed) -> Vector{Float64}

Solve the linear system with fixed variables set to their current values in x.

For each fixed variable at index i, modifies the system by:
1. Setting row i of A to have only a 1 on the diagonal (A[i,i] = 1)
2. Setting b[i] to the fixed value x[i]
3. This constrains x[i] to equal the rounded/fixed value

Then solves the modified system for all variables.
"""
function solve_with_fixed_variables(
    A::SparseMatrixCSC,
    b::Vector{Float64},
    x::Vector{Float64},
    fixed::BitVector
)
    n = length(b)
    
    # Identify fixed variables
    fixed_vars = findall(fixed)
    
    if isempty(fixed_vars)
        # No variables are fixed, nothing to modify
        return copy(x)
    end
    
    # Create modified system by converting sparse matrix to dense for easier manipulation
    # (For large systems, could maintain sparsity with more careful bookkeeping)
    A_modified = Matrix(A)
    b_modified = copy(b)
    
    # For each fixed variable, zero out its row and set diagonal to 1
    for idx in fixed_vars
        # Zero out entire row
        A_modified[idx, :] .= 0.0
        
        # Set diagonal to 1
        A_modified[idx, idx] = 1.0
        
        # Set RHS to the fixed value
        b_modified[idx] = x[idx]
    end
    
    # Solve the modified system
    # Convert back to sparse for efficiency
    A_modified_sparse = sparse(A_modified)
    
    # Solve (use least-squares to handle any potential singularities)
    x_new = 0
    try
        x_new = A_modified_sparse \ b_modified
    catch e
        # Fall back to QR-based least-squares if singular
        x_new = qr(A_modified) \ b_modified
    end
    
    return x_new
end


"""
    round_solution(x::Vector{Float64}, k::Int) -> Vector{Float64}

Simple rounding: round the first k entries to nearest integers.
This is a baseline comparison method.
"""
function round_solution(x::Vector{Float64}, k::Int)
    x_rounded = copy(x)
    x_rounded[1:k] = round.(x[1:k])
    return x_rounded
end
