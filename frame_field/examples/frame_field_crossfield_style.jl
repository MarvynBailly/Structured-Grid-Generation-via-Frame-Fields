"""
    frame_field_crossfield_style.jl

Frame field generation using the CrossField.jl assembly and solver approach.
This mimics the solver strategy from the MIQ_old implementation.
"""

using Pkg
Pkg.activate(".")

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using GeometryBasics
using CairoMakie
using LinearAlgebra
using SparseArrays
using DataStructures

include("../src/meshio.jl")
include("../src/mesh_types.jl")
include("../src/dijkstra_forest.jl")
include("../src/plotters.jl")

folder_name = "frame_field_crossfield"
case_name = "simplest-square"

"""
    assemble_system_crossfield_style(field, topology, fixed_edges_per_face, constrained_faces, constrained_angles)

Assemble system matrix using the CrossField.jl approach:
- Only free variables appear in the system
- Fixed edges and constrained faces are handled via RHS adjustments
- Uses the gradient equations directly (matching Equations 2 & 3 from CrossField.jl)
"""
function assemble_system_crossfield_style(field, topology, fixed_edges_per_face, 
                                          constrained_faces=Int[], constrained_angles=Float64[])
    n_faces = topology.n_faces
    
    # Build edge index mapping
    edge_to_idx = Dict{Tuple{Int,Int}, Int}()
    edge_idx = 0
    for face_i in 1:n_faces
        for (face_j, edge) in topology.dual_adj[face_i]
            if face_i < face_j
                edge_idx += 1
                edge_to_idx[(face_i, face_j)] = edge_idx
            end
        end
    end
    n_edges = edge_idx
    
    # Determine fixed edges
    fixed_edges_set = Set(values(fixed_edges_per_face))
    fixed_edge_flags = falses(n_edges)
    for face_i in 1:n_faces
        for (face_j, v_edge) in topology.dual_adj[face_i]
            if face_i < face_j && v_edge in fixed_edges_set
                e_idx = edge_to_idx[(face_i, face_j)]
                fixed_edge_flags[e_idx] = true
            end
        end
    end
    
    # Create constraint angle dictionary
    constrained_dict = Dict{Int, Float64}()
    for (idx, f_idx) in enumerate(constrained_faces)
        constrained_dict[f_idx] = constrained_angles[idx]
    end
    
    # Map Variables (only FREE variables)
    theta_var_map = Dict{Int, Int}()
    p_var_map = Dict{Tuple{Int,Int}, Int}()
    current_idx = 0
    
    # Theta Variables (Free Faces only)
    for f_idx in 1:n_faces
        if !haskey(constrained_dict, f_idx)
            current_idx += 1
            theta_var_map[f_idx] = current_idx
        end
    end
    
    # P Variables (Free Edges only)
    for face_i in 1:n_faces
        for (face_j, v_edge) in topology.dual_adj[face_i]
            if face_i < face_j
                e_idx = edge_to_idx[(face_i, face_j)]
                
                # Skip if: boundary, fixed tree edge, or between two constrained faces
                is_boundary = length(topology.dual_adj[face_i]) < 3 || length(topology.dual_adj[face_j]) < 3
                is_between_constraints = haskey(constrained_dict, face_i) && haskey(constrained_dict, face_j)
                
                if !fixed_edge_flags[e_idx] && !is_boundary && !is_between_constraints
                    current_idx += 1
                    p_var_map[(face_i, face_j)] = current_idx
                end
            end
        end
    end
    
    total_vars = current_idx
    
    println("\nVariable mapping (CrossField style):")
    println("  Free theta variables: $(length(theta_var_map))")
    println("  Free p variables: $(length(p_var_map))")
    println("  Total variables: $total_vars")
    
    # Build system using triplet format
    I, J, V = Int[], Int[], Float64[]
    b = zeros(total_vars)
    
    function add!(r, c, v)
        push!(I, r); push!(J, c); push!(V, v)
    end
    
    # Build Energy Gradient (matching CrossField.jl equations)
    for face_i in 1:n_faces
        for (face_j, v_edge) in topology.dual_adj[face_i]
            if face_i >= face_j; continue; end  # Process each edge once
            
            kappa_ij = get_transport_angle(field, face_i, face_j)
            
            # Get variable indices (0 if fixed/constrained)
            row_i = get(theta_var_map, face_i, 0)
            row_j = get(theta_var_map, face_j, 0)
            row_p = get(p_var_map, (face_i, face_j), 0)
            
            # Get fixed/constrained values
            theta_i_val = get(constrained_dict, face_i, 0.0)
            theta_j_val = get(constrained_dict, face_j, 0.0)
            
            # Determine P value for fixed edges
            p_val = 0.0
            if row_p == 0
                e_idx = edge_to_idx[(face_i, face_j)]
                if fixed_edge_flags[e_idx]
                    p_val = 0.0
                elseif haskey(constrained_dict, face_i) && haskey(constrained_dict, face_j)
                    # Round between constraints
                    term = (2.0/π) * (theta_j_val - theta_i_val - kappa_ij)
                    p_val = round(term)
                end
            end
            
            # RHS accumulation
            const_RHS = -kappa_ij
            if row_i == 0; const_RHS -= theta_i_val; end
            if row_j == 0; const_RHS += theta_j_val; end
            if row_p == 0; const_RHS -= (π/2) * p_val; end
            
            # Add derivatives (from ∂E/∂θ_i, ∂E/∂θ_j, ∂E/∂p)
            # These follow the CrossField.jl equations exactly
            
            # Row i (∂E/∂θ_i contribution)
            if row_i > 0
                add!(row_i, row_i, 2.0)
                if row_j > 0; add!(row_i, row_j, -2.0); end
                if row_p > 0; add!(row_i, row_p, π); end
                b[row_i] += 2.0 * const_RHS
            end
            
            # Row j (∂E/∂θ_j contribution)
            if row_j > 0
                if row_i > 0; add!(row_j, row_i, -2.0); end
                add!(row_j, row_j, 2.0)
                if row_p > 0; add!(row_j, row_p, -π); end
                b[row_j] += -2.0 * const_RHS
            end
            
            # Row p (∂E/∂p contribution)
            if row_p > 0
                if row_i > 0; add!(row_p, row_i, π); end
                if row_j > 0; add!(row_p, row_j, -π); end
                add!(row_p, row_p, π^2/2.0)
                b[row_p] += π * const_RHS
            end
        end
    end
    
    A = sparse(I, J, V, total_vars, total_vars)
    return A, b, theta_var_map, p_var_map
end

"""
    solve_greedy_with_queue(A, b, p_var_map; verbose=true)

Greedy solver using queue-based Gauss-Seidel relaxation (CrossField.jl style).
"""
function solve_greedy_with_queue(A, b, p_var_map; verbose=true)
    if verbose
        println("\nSolving continuous relaxation...")
    end
    
    # Initial solve
    x = A \ b
    
    if verbose
        println("  Initial solution computed")
        p_indices = collect(values(p_var_map))
        if !isempty(p_indices)
            println("  First few p values: ", x[p_indices[1:min(5, length(p_indices))]])
        end
    end
    
    int_indices = collect(values(p_var_map))
    if isempty(int_indices)
        return x
    end
    
    fixed_mask = falses(length(b))
    n_vars = length(b)
    
    # Precompute diagonal for Gauss-Seidel
    diag_A = Vector{Float64}(undef, n_vars)
    for i in 1:n_vars
        diag_A[i] = A[i, i]
    end
    
    if verbose
        println("\nGreedy rounding with queue-based relaxation...")
    end
    
    n_rounded = 0
    while true
        # Find best variable to round
        best_idx = -1
        min_err = Inf
        best_val = 0.0
        
        for idx in int_indices
            if !fixed_mask[idx]
                v = x[idx]
                r = round(v)
                err = abs(v - r)
                if err < min_err
                    min_err = err
                    best_idx = idx
                    best_val = r
                end
            end
        end
        
        if best_idx == -1; break; end
        
        # Fix variable
        x[best_idx] = best_val
        fixed_mask[best_idx] = true
        n_rounded += 1
        
        if verbose && (n_rounded % 5 == 0 || n_rounded <= 5)
            println("  Rounded $n_rounded/$(length(int_indices)): x[$best_idx] = $best_val (error: $(round(min_err, digits=8)))")
        end
        
        # Queue-based Gauss-Seidel relaxation
        n_iters = local_gauss_seidel_queue!(A, b, x, fixed_mask, best_idx, diag_A)
        
        if verbose && n_iters > 100
            println("    Relaxed with $n_iters iterations")
        end
    end
    
    if verbose
        println("\nGreedy rounding complete!")
        println("  All $(length(int_indices)) integer variables fixed")
    end
    
    return x
end

"""
    local_gauss_seidel_queue!(A, b, x, fixed_mask, start_node, diag_A; max_iter=50000, tol=1e-5)

Queue-based local Gauss-Seidel update starting from a changed variable.
"""
function local_gauss_seidel_queue!(A, b, x, fixed_mask, start_node, diag_A; max_iter=50000, tol=1e-5)
    q = Queue{Int}()
    in_queue = Set{Int}()
    
    rows = rowvals(A)
    vals = nonzeros(A)
    
    # Add neighbors of start_node to queue
    for idx in nzrange(A, start_node)
        neighbor = rows[idx]
        if !fixed_mask[neighbor] && !(neighbor in in_queue)
            enqueue!(q, neighbor)
            push!(in_queue, neighbor)
        end
    end
    
    iter = 0
    while !isempty(q) && iter < max_iter
        iter += 1
        k = dequeue!(q)
        delete!(in_queue, k)
        
        if !fixed_mask[k]
            # Compute residual: r_k = b[k] - Σ A_kj * x_j
            Ax_k = 0.0
            for idx in nzrange(A, k)
                Ax_k += vals[idx] * x[rows[idx]]
            end
            r_k = b[k] - Ax_k
            
            if abs(r_k) > tol
                # Update variable
                x[k] += r_k / diag_A[k]
                
                # Add neighbors to queue
                for idx in nzrange(A, k)
                    neighbor = rows[idx]
                    if !fixed_mask[neighbor] && !(neighbor in in_queue)
                        enqueue!(q, neighbor)
                        push!(in_queue, neighbor)
                    end
                end
            end
        end
    end
    
    return iter
end

function main()
    println("="^60)
    println("Frame Field Generation - CrossField.jl Style Solver")
    println("="^60)
    
    println("\nLoading $case_name mesh...")
    mesh_file = joinpath(@__DIR__, "..", "..", "triangulations", "$case_name.msh")
    mesh = load_triangulation(mesh_file)
    
    println("Building mesh topology...")
    topology = build_mesh_topology(mesh)
    
    println("Creating frame field...")
    field = CrossField(topology)
    
    # Constraints
    constrained_faces = [1]
    constrained_angles = [0.0]
    
    println("\nConstrained faces: $constrained_faces")
    println("Constrained angles: $(rad2deg.(constrained_angles)) degrees")
    
    # Compute spanning forest
    println("\nComputing spanning forest...")
    potential_fixed_edges = compute_spanning_forest(mesh; constrained_faces=constrained_faces)
    println("Spanning forest has $(length(potential_fixed_edges)) edges")
    
    # Fix suitable edges
    println("\nFixing suitable edges...")
    fixed_edges_per_face = potential_fixed_edges#fix_suitable_edges(mesh, potential_fixed_edges)
    println("Fixed $(length(fixed_edges_per_face)) edges (one per face)")
    
    # Assemble system using CrossField.jl style
    println("\n" * "="^60)
    println("ASSEMBLING SYSTEM (CrossField.jl Style)")
    println("="^60)
    
    A, b, theta_var_map, p_var_map = assemble_system_crossfield_style(
        field, topology, fixed_edges_per_face, constrained_faces, constrained_angles
    )
    
    println("\nSystem dimensions:")
    println("  Matrix A: $(size(A))")
    println("  Vector b: $(length(b))")
    
    # Solve
    println("\n" * "="^60)
    println("SOLVING")
    println("="^60)
    
    x_solution = solve_greedy_with_queue(A, b, p_var_map; verbose=true)
    
    # Extract results
    println("\n" * "="^60)
    println("EXTRACTING RESULTS")
    println("="^60)
    
    n_faces = topology.n_faces
    angles = zeros(Float64, n_faces)
    
    # Fill free theta variables
    for (f_idx, var_idx) in theta_var_map
        angles[f_idx] = x_solution[var_idx]
    end
    
    # Fill constrained angles
    for (i, f_idx) in enumerate(constrained_faces)
        angles[f_idx] = constrained_angles[i]
    end
    
    # Extract period jumps
    for ((face_i, face_j), var_idx) in p_var_map
        p_val = Int(round(x_solution[var_idx]))
        set_period_jump!(field, face_i, face_j, p_val)
    end
    
    # Update field
    field.theta = angles
    
    # Compute energy
    println("\nComputing smoothness energy...")
    final_energy = compute_smoothness_energy(field)
    println("  Final smoothness energy: $final_energy")
    
    # Verify solution
    println("\nVerifying solution...")
    residual = A * x_solution - b
    max_residual = maximum(abs.(residual))
    println("  Maximum residual: $max_residual")
    
    # Print period jumps
    println("\nPeriod jumps:")
    for ((face_i, face_j), var_idx) in sort(collect(p_var_map), by=x->x[2])
        p_val = Int(round(x_solution[var_idx]))
        println("  p_($face_i,$face_j) = $p_val")
    end
    
    # Save visualizations
    println("\n" * "="^60)
    println("CREATING VISUALIZATIONS")
    println("="^60)
    
    output_dir = joinpath(@__DIR__, "..", "output", folder_name)
    mkpath(output_dir)
    
    # Plot forest
    forest_output = joinpath(output_dir, "forest_$(case_name).png")
    plot_forest(mesh; constrained_faces=constrained_faces, 
                potential_fixed_edges=potential_fixed_edges,
                savepath=forest_output)
    println("\nForest visualization saved to: $forest_output")
    
    # Plot frame field
    field_output = joinpath(output_dir, "frame_field_$(case_name).png")
    plot_frame_field(mesh, angles; savepath=field_output, 
                     constrained_faces=constrained_faces,
                     arrow_scale=0.4)
    println("Frame field visualization saved to: $field_output")
    
    println("\n" * "="^60)
    println("DONE!")
    println("="^60)
end

main()
