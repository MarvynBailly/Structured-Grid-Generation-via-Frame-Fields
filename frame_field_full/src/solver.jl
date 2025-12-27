module Solver

using ..Types
using ..Topology
using ..Analysis
using LinearAlgebra
using SparseArrays
using Printf
using DataStructures
using Random

export initialize_field, build_spanning_tree!, assemble_system_weighted, local_gauss_seidel_queue!, solve_greedy!, optimize_singularities!

function initialize_field(topo, constraints)
    kappa = compute_transport_angles(topo)
    n = length(topo.faces); theta = zeros(n); fixed = Set{Tuple{Int,Int}}()
    p_jumps = Dict{Tuple{Int,Int}, Int}()
    for (f, a) in constraints; theta[f] = a; end
    return CrossField(topo, theta, p_jumps, kappa, constraints, fixed)
end

function build_spanning_tree!(field; verbose=false)
    visited = falses(length(field.topology.faces))
    queue = Int[]
    
    if verbose
        println("\n--- Building Spanning Forest ---")
    end
    
    n_roots = 0
    if !isempty(field.constrained_faces)
        for f in keys(field.constrained_faces)
            push!(queue, f)
            visited[f] = true
            n_roots += 1
        end
        if verbose
            println("Starting from $n_roots constrained faces (roots)")
        end
    else
        push!(queue, 1)
        visited[1] = true
        n_roots = 1
        if verbose
            println("No constrained faces - starting from face 1")
        end
    end
    
    while !isempty(queue)
        u = popfirst!(queue)
        for (v, _) in field.topology.dual_adj[u]
            if !visited[v]
                visited[v] = true
                push!(queue, v)
                key = minmax(u, v)
                push!(field.fixed_edges, key)
                field.period_jumps[key] = 0
            end
        end
    end
end

function assemble_system_weighted(field)
    topo = field.topology
    n_faces = length(topo.faces)
    
    theta_map = Dict{Int,Int}()
    p_map = Dict{Tuple{Int,Int}, Int}()
    curr_idx = 0
    
    for f in 1:n_faces
        if !haskey(field.constrained_faces, f)
            curr_idx += 1; theta_map[f] = curr_idx
        end
    end
    
    for i in 1:n_faces
        for (j, _) in topo.dual_adj[i]
            if i < j
                edge = (i, j)
                is_constrained_edge = haskey(field.constrained_faces, i) && haskey(field.constrained_faces, j)
                if !(edge in field.fixed_edges) && !is_constrained_edge
                    curr_idx += 1; p_map[edge] = curr_idx
                end
            end
        end
    end
    
    n_vars = curr_idx
    I = Int[]; J = Int[]; V = Float64[]; b = zeros(n_vars)
    
    function add!(r, c, val)
        push!(I, r); push!(J, c); push!(V, val)
    end
    
    for i in 1:n_faces
        for (j, _) in topo.dual_adj[i]
            if i >= j continue end
            edge = (i, j)
            kappa = field.transport_angles[edge]
            
            idx_i = get(theta_map, i, 0)
            idx_j = get(theta_map, j, 0)
            idx_p = get(p_map, edge, 0)
            
            val_i = get(field.constrained_faces, i, 0.0)
            val_j = get(field.constrained_faces, j, 0.0)
            
            val_p = 0.0
            if idx_p == 0
                if edge in field.fixed_edges
                    val_p = Float64(field.period_jumps[edge])
                elseif haskey(field.constrained_faces, i) && haskey(field.constrained_faces, j)
                    term = (2.0/π) * (val_j - val_i - kappa)
                    val_p = round(term)
                    field.period_jumps[edge] = Int(val_p)
                end
            end
            
            vi = topo.faces[i]; vj = topo.faces[j]; shared = intersect(vi, vj)
            p1 = topo.vertices[shared[1]]; p2 = topo.vertices[shared[2]]
            len_sq = (p1.x-p2.x)^2 + (p1.y-p2.y)^2 + (p1.z-p2.z)^2
            weight = 1.0
            
            rhs_const = -kappa
            if idx_i == 0; rhs_const -= val_i; end
            if idx_j == 0; rhs_const += val_j; end
            if idx_p == 0; rhs_const -= (π/2) * val_p; end
            
            if idx_i > 0
                add!(idx_i, idx_i, 2.0 * weight)
                if idx_j > 0; add!(idx_i, idx_j, -2.0 * weight); end
                if idx_p > 0; add!(idx_i, idx_p, π * weight); end
                b[idx_i] += 2.0 * rhs_const * weight
            end
            
            if idx_j > 0
                if idx_i > 0; add!(idx_j, idx_i, -2.0 * weight); end
                add!(idx_j, idx_j, 2.0 * weight)
                if idx_p > 0; add!(idx_j, idx_p, -π * weight); end
                b[idx_j] += -2.0 * rhs_const * weight
            end
            
            if idx_p > 0
                if idx_i > 0; add!(idx_p, idx_i, π * weight); end
                if idx_j > 0; add!(idx_p, idx_j, -π * weight); end
                add!(idx_p, idx_p, (π^2/2.0) * weight)
                b[idx_p] += π * rhs_const * weight
            end
        end
    end
    
    A = sparse(I, J, V, n_vars, n_vars)
    return A, b, theta_map, p_map
end

function local_gauss_seidel_queue!(A, b, x, fixed_mask, start_node, diag_A)
    q = Queue{Int}(); in_q = Set{Int}()
    rows = rowvals(A); vals = nonzeros(A)
    
    for idx in nzrange(A, start_node)
        n = rows[idx]
        if !fixed_mask[n] && !(n in in_q)
            enqueue!(q, n); push!(in_q, n)
        end
    end
    
    iter = 0; max_iter = 50000; tol = 1e-5
    while !isempty(q) && iter < max_iter
        iter += 1
        k = dequeue!(q); delete!(in_q, k)
        
        Ax_k = 0.0
        for idx in nzrange(A, k)
            Ax_k += vals[idx] * x[rows[idx]]
        end
        r_k = b[k] - Ax_k
        
        if abs(r_k) > tol
            x[k] += r_k / diag_A[k]
            for idx in nzrange(A, k)
                n = rows[idx]
                if !fixed_mask[n] && !(n in in_q)
                    enqueue!(q, n); push!(in_q, n)
                end
            end
        end
    end
    return iter
end

function solve_greedy!(field; verbose=false, callback=nothing)
    if verbose
        println("--- Starting Unified MIQ Solver ---")
    end
    build_spanning_tree!(field; verbose=verbose)
    
    if verbose
        println("Assembling system...")
    end
    A, b, theta_map, p_map = assemble_system_weighted(field)
    n_vars = length(b)
    
    if verbose
        println("Initial solve...")
    end
    x = A \ b
    
    for (f, idx) in theta_map
        field.theta[f] = x[idx]
    end
    for (edge, idx) in p_map
        raw_p = Int(round(x[idx]))
        field.period_jumps[edge] = mod(raw_p + 2, 4) - 2
    end
    
    !isnothing(callback) && callback(field, 0, 0, 0.0)

    int_indices = collect(values(p_map))
    fixed_mask = falses(n_vars)
    diag_A = [A[i,i] for i in 1:n_vars]
    
    num_p = length(int_indices)
    if verbose
        println("Greedy rounding $num_p integer variables...")
    end
    
    count = 0
    while true
        best_idx = -1; min_err = Inf; best_val = 0.0
        
        for idx in int_indices
            if !fixed_mask[idx]
                v = x[idx]; r = round(v); err = abs(v-r)
                if err < min_err
                    min_err = err; best_idx = idx; best_val = r
                end
            end
        end
        
        if best_idx == -1; break; end
        
        x[best_idx] = best_val
        fixed_mask[best_idx] = true
        count += 1
        
        local_gauss_seidel_queue!(A, b, x, fixed_mask, best_idx, diag_A)
        
        !isnothing(callback) && callback(field, count, num_p, min_err)

        if verbose && count % 50 == 0
            @printf("  Rounded %d/%d (err: %.4f)\n", count, num_p, min_err)
        end
    end
    
    for (f, idx) in theta_map
        field.theta[f] = x[idx]
    end
    for (edge, idx) in p_map
        raw_p = Int(round(x[idx]))
        field.period_jumps[edge] = mod(raw_p + 2, 4) - 2
    end
    
    if verbose
        println("--- Solver Complete ---")
    end
end

function assemble_laplacian_fixed_p(field)
    topo = field.topology
    n_faces = length(topo.faces)
    
    # We solve for theta only.
    # We map face indices to system indices (excluding constrained faces).
    theta_map = Dict{Int,Int}()
    curr_idx = 0
    for f in 1:n_faces
        if !haskey(field.constrained_faces, f)
            curr_idx += 1
            theta_map[f] = curr_idx
        end
    end
    n_vars = curr_idx
    
    I = Int[]; J = Int[]; V = Float64[]; 
    b = zeros(n_vars)
    
    # Pre-calculate energy offset (constant part) for return if needed, 
    # but strictly we only need A and b for the solve.
    
    for i in 1:n_faces
        for (j, _) in topo.dual_adj[i]
            if i >= j; continue; end
            edge = (i, j)
            
            # Get fixed values
            kappa = field.transport_angles[edge]
            p = get(field.period_jumps, edge, 0)
            
            # Target difference: theta_j - theta_i = kappa + p*pi/2
            # Energy term: (theta_j - theta_i - target)^2
            target = kappa + p * (π/2.0)
            
            idx_i = get(theta_map, i, 0)
            idx_j = get(theta_map, j, 0)
            val_i = get(field.constrained_faces, i, 0.0)
            val_j = get(field.constrained_faces, j, 0.0)
            
            # Weights (simple uniform for now, can be area-weighted)
            weight = 1.0
            
            # Formulate: r = (theta_j - theta_i) - target
            # minimize r^2 -> derivative linear system
            
            # RHS contributions from constraints and target
            # term for i: -2 * (theta_j - theta_i - target) * (-1) = 0
            # -> 2*(theta_j - theta_i) = 2*target
            
            rhs_term = target
            if idx_i == 0; rhs_term += val_i; end
            if idx_j == 0; rhs_term -= val_j; end
            
            if idx_i > 0
                push!(I, idx_i); push!(J, idx_i); push!(V, 2.0 * weight)
                b[idx_i] -= 2.0 * rhs_term * weight
                if idx_j > 0
                    push!(I, idx_i); push!(J, idx_j); push!(V, -2.0 * weight)
                end
            end
            
            if idx_j > 0
                push!(I, idx_j); push!(J, idx_j); push!(V, 2.0 * weight)
                b[idx_j] += 2.0 * rhs_term * weight
                if idx_i > 0
                    push!(I, idx_j); push!(J, idx_i); push!(V, -2.0 * weight)
                end
            end
        end
    end
    
    A = sparse(I, J, V, n_vars, n_vars)
    
    # Add small regularization to diagonal to ensure positive definite
    # (Fixes floating islands or gauge freedom if no Dirichlet constraints exist)
    if isempty(field.constrained_faces)
        # Fix one variable arbitrarily or add regularization
        # Adding epsilon to diagonal is safer for "soft" fixing
        A += sparse(I, I, 1e-6 * ones(length(I)), n_vars, n_vars)
    end

    return A, b, theta_map
end

function compute_energy(field)
    energy = 0.0
    for i in 1:length(field.topology.faces)
        for (j, _) in field.topology.dual_adj[i]
            if i >= j; continue; end
            
            theta_i = field.theta[i]
            theta_j = field.theta[j]
            edge = (i, j)
            kappa = field.transport_angles[edge]
            p = get(field.period_jumps, edge, 0)
            
            diff = theta_j - theta_i - kappa - p * (π/2.0)
            energy += diff^2
        end
    end
    return energy
end

function optimize_singularities!(field; max_iters=100, verbose=false)
    if verbose
        println("\n--- Starting Local Singularity Optimization ---")
    end
    
    # 1. Precompute Factorization
    A, b_base, theta_map = assemble_laplacian_fixed_p(field)
    
    # Use Cholesky (LDLt) if possible, or LU if not perfectly symmetric/PD
    # A should be symmetric.
    F = cholesky(Hermitian(A)) 
    
    current_energy = compute_energy(field)
    if verbose
        println("Initial Energy: $current_energy")
    end
    
    topo = field.topology
    
    # Map faces to system b indices for fast updates
    # We need to know which entries in 'b' correspond to an edge (i,j)
    # Re-building the sparse map is slow, so we do it on the fly or lookup.
    
    for iter in 1:max_iters
        changed = false
        
        # Identify singularities
        sings = compute_singularities(field; verbose=false)
        shuffle!(sings) # Randomize order to prevent cycles
        
        for (v_idx, idx_val) in sings
            # Get star of faces
            # (Re-using logic from analysis, or just iterating adjacent faces)
            
            # Simple approach: Iterate over edges connected to v_idx
            # We need the faces around v_idx to find the edges radiating from it.
            # But "moving" a singularity effectively means changing p on one of the radiating edges.
            
            # Let's iterate over all edges connected to vertex v
            # We need edges incident to v.
            
            incident_faces = Int[] 
            # Building v2f is fast enough.
            v2f = build_vertex_to_faces(topo)
            faces_v = v2f[v_idx]
            
            # Find edges connected to v in this one-ring
            # An edge (i, j) radiates from v if v is one of the shared vertices.
            candidate_edges = Tuple{Int,Int}[]
            for f in faces_v
                for (n, _) in topo.dual_adj[f]
                    if n > f # Unique edges
                        # Check if shared vertices include v_idx
                        v_shared = intersect(topo.faces[f], topo.faces[n])
                        if v_idx in v_shared
                            push!(candidate_edges, (f, n))
                        end
                    end
                end
            end
            
            for edge in candidate_edges
                # Try increasing and decreasing p
                original_p = get(field.period_jumps, edge, 0)
                
                # Movements: +1 and -1 change to p
                for delta in [-1, 1]
                    new_p = original_p + delta
                    
                    # 1. Update p
                    field.period_jumps[edge] = new_p
                    
                    # 2. Update b (RHS) efficiently
                    # Re-calculating full b is safer for implementation speed vs complexity right now
                    # optimizing this requires delta-updates to b, which is easy but verbose.
                    _, b_new, _ = assemble_laplacian_fixed_p(field)
                    
                    # 3. Solve
                    x_new = F \ b_new
                    
                    # 4. Update field temporarily
                    old_thetas = copy(field.theta)
                    for (f, idx) in theta_map
                        field.theta[f] = x_new[idx]
                    end
                    
                    # 5. Check Energy
                    new_energy = compute_energy(field)
                    
                    if new_energy < current_energy - 1e-6 # Epsilon improvement
                        current_energy = new_energy
                        changed = true
                        if verbose
                            @printf("  Moved singularity at v%d: E %.4f -> %.4f\n", v_idx, current_energy + 1e-6, current_energy)
                        end
                        # Keep change, break to next singularity (greedy)
                        break 
                    else
                        # Revert
                        field.period_jumps[edge] = original_p
                        field.theta = old_thetas
                    end
                end
                
                if changed; break; end
            end
        end
        
        if !changed
            if verbose 
                println("Converged after $iter iterations.")
                println("Optimization Singularities found:")
                for (v_idx, I) in sings
                    @printf("  Vertex %d: Index = %.4f\n", v_idx, I)
                end
                println("Final Energy: $current_energy\n")
            end

            break
        end
    end
end




end
