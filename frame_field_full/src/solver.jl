module Solver

using ..Types
using ..Topology
using LinearAlgebra
using SparseArrays
using Printf
using DataStructures

export initialize_field, build_spanning_tree!, assemble_system_weighted, local_gauss_seidel_queue!, solve_greedy!

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

end
