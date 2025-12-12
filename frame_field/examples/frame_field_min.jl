using Pkg
Pkg.activate(".")


function assemble_energy_system(field, topology, fixed_edges_per_face, constrained_faces=Int[], constrained_angles=Float64[])
    n_faces = topology.n_faces
    
    # Build mapping from face to constrained angle
    constrained_theta = Dict{Int, Float64}()
    for (i, face_idx) in enumerate(constrained_faces)
        constrained_theta[face_idx] = constrained_angles[i]
    end
    
    # Collect all dual edges
    all_dual_edges = Tuple{Int, Int}[]
    for face_i in 1:n_faces
        for (face_j, edge) in topology.dual_adj[face_i]
            if face_i < face_j
                push!(all_dual_edges, (face_i, face_j))
            end
        end
    end
    
    # Build edge_to_dual mapping (vertex edge -> dual edge)
    edge_to_dual = Dict{Tuple{Int,Int}, Tuple{Int,Int}}()
    for face_i in 1:n_faces
        for (face_j, edge) in topology.dual_adj[face_i]
            if face_i < face_j
                edge_to_dual[edge] = (face_i, face_j)
            end
        end
    end
    
    # Determine which dual edges are fixed
    fixed_edges = Set(values(fixed_edges_per_face))
    fixed_dual_edges = Set{Tuple{Int,Int}}()
    for v_edge in fixed_edges
        if haskey(edge_to_dual, v_edge)
            push!(fixed_dual_edges, edge_to_dual[v_edge])
        end
    end
    
    # Separate free and fixed dual edges
    free_dual_edges = [e for e in all_dual_edges if !(e in fixed_dual_edges)]
    fixed_dual_edges_list = [e for e in all_dual_edges if e in fixed_dual_edges]
    
    n_free_p = length(free_dual_edges)
    n_fixed_p = length(fixed_dual_edges_list)
    n_all_p = n_free_p + n_fixed_p
    n_theta_vars = n_faces
    n_vars = n_all_p + n_theta_vars
    n_equations = n_all_p + n_faces  # One equation per edge + one per face
    
    A = spzeros(Float64, n_equations, n_vars)
    b = zeros(Float64, n_equations)
    
    # Map (face_i, face_j) to p variable index
    # Order: [free p vars, fixed p vars, theta vars]
    p_var_idx = Dict{Tuple{Int,Int}, Int}()
    var_to_p = Dict{Int, Tuple{Int,Int}}()
    
    # Free p variables (indices 1:n_free_p)
    for (idx, (face_i, face_j)) in enumerate(free_dual_edges)
        p_var_idx[(face_i, face_j)] = idx
        var_to_p[idx] = (face_i, face_j)
    end
    
    # Fixed p variables (indices n_free_p+1:n_all_p)
    for (idx, (face_i, face_j)) in enumerate(fixed_dual_edges_list)
        var_idx = n_free_p + idx
        p_var_idx[(face_i, face_j)] = var_idx
        var_to_p[var_idx] = (face_i, face_j)
    end
    
    # Map theta variables (indices n_all_p+1:n_vars)
    # Separate free and constrained theta variables
    free_theta_faces = [i for i in 1:n_faces if !haskey(constrained_theta, i)]
    fixed_theta_faces = [i for i in 1:n_faces if haskey(constrained_theta, i)]
    
    n_free_theta = length(free_theta_faces)
    n_fixed_theta = length(fixed_theta_faces)
    
    # Map: variable index -> face index
    var_to_theta = Dict{Int, Int}()
    # Map: face index -> variable index
    face_to_theta_var = Dict{Int, Int}()
    
    # Free theta variables come first
    for (idx, face_i) in enumerate(free_theta_faces)
        var_idx = n_all_p + idx
        var_to_theta[var_idx] = face_i
        face_to_theta_var[face_i] = var_idx
    end
    
    # Fixed theta variables come after
    for (idx, face_i) in enumerate(fixed_theta_faces)
        var_idx = n_all_p + n_free_theta + idx
        var_to_theta[var_idx] = face_i
        face_to_theta_var[face_i] = var_idx
    end
    
    row = 0
    
    # Equations for FREE edges: ∂E/∂p_ij = π(θ_i + κ_ij + (π/2)p_ij - θ_j) = 0
    # Simplifies to: (π/2)p_ij + θ_i - θ_j = -κ_ij
    for (face_i, face_j) in free_dual_edges
        row += 1
        
        p_idx = p_var_idx[(face_i, face_j)]
        theta_i_idx = face_to_theta_var[face_i]
        theta_j_idx = face_to_theta_var[face_j]
        
        A[row, p_idx] = π/2
        A[row, theta_i_idx] = 1.0
        A[row, theta_j_idx] = -1.0
        
        kappa_ij = get_transport_angle(field, face_i, face_j)
        b[row] = -kappa_ij
    end
    
    # Equations for FIXED edges: p_ij = 0
    for (face_i, face_j) in fixed_dual_edges_list
        row += 1
        
        p_idx = p_var_idx[(face_i, face_j)]
        A[row, p_idx] = 1.0
        b[row] = 0.0
    end
    
    # Equations for FREE face angles: ∂E/∂θ_i = Σ_j 2(θ_i + κ_ij + (π/2)p_ij - θ_j) = 0
    # For each free face i, sum over all neighbors j
    for face_i in free_theta_faces
        row += 1
        theta_i_idx = face_to_theta_var[face_i]
        
        # Degree of face i
        degree = length(topology.dual_adj[face_i])
        A[row, theta_i_idx] = degree
        
        rhs = 0.0
        for (face_j, edge) in topology.dual_adj[face_i]
            # Determine p_ij index and sign
            if face_i < face_j
                p_idx = p_var_idx[(face_i, face_j)]
                p_sign = 1.0
                kappa_ij = get_transport_angle(field, face_i, face_j)
            else
                p_idx = p_var_idx[(face_j, face_i)]
                p_sign = -1.0
                kappa_ij = get_transport_angle(field, face_i, face_j)
            end
            
            A[row, p_idx] += (π/2) * p_sign
            
            theta_j_idx = face_to_theta_var[face_j]
            A[row, theta_j_idx] -= 1.0
            
            rhs -= 1.0 * kappa_ij
        end
        b[row] = rhs
    end
    
    # Equations for FIXED face angles: θ_i = constrained_value
    for face_i in fixed_theta_faces
        row += 1
        theta_i_idx = face_to_theta_var[face_i]
        
        A[row, theta_i_idx] = 1.0
        b[row] = constrained_theta[face_i]
    end
    
    return A, b, var_to_p, var_to_theta, n_free_p, n_free_theta
end



push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using GeometryBasics
using CairoMakie
using LinearAlgebra
using SparseArrays
using FileIO

include("../src/meshio.jl")
include("../src/mesh_types.jl")
include("../src/dijkstra_forest.jl")
include("../src/cross_field_energy.jl")
include("../src/plotters.jl")
include("../src/greedy_solver.jl")

folder_name = "frame_field_min"
case_name = "three-triangles"
# case_name = "simplest-square"
# case_name = "simple-square"
# case_name = "regular-square-4x4"
# case_name = "regular-square-8x8"
# case_name = "disk-radial"
# case_name = "disk-radial-fine"
# case_name = "300_polygon_sphere_100mm"

function main()
    println("Loading $case_name mesh...")
    mesh_file = joinpath(@__DIR__, "..", "..", "triangulations", "$case_name.msh")
    mesh = load_triangulation(mesh_file)
    
    println("Building mesh topology...")
    topology = build_mesh_topology(mesh)
    
    println("Creating frame field...")
    field = CrossField(topology)
    
    # Constrain one face to have a specific direction
    constrained_faces = [1]  # Constrain the first face
    constrained_angles = [0.0]  # Fix angle to 0 (horizontal direction)
    
    println("\nConstrained faces: $constrained_faces")
    println("Constrained angles: $(rad2deg.(constrained_angles)) degrees")
    
    # Compute spanning forest
    println("\nComputing spanning forest...")
    potential_fixed_edges = compute_spanning_forest(mesh; constrained_faces=constrained_faces)
    println("Spanning forest has $(length(potential_fixed_edges)) edges")
    
    fixed_edges_per_face = fix_suitable_edges(mesh, potential_fixed_edges)
    println("Fixed $(length(fixed_edges_per_face)) edges (one per face)")
    
    # Save forest visualization
    println("\nSaving forest visualization...")
    output_dir = joinpath(@__DIR__, "..", "output", folder_name)
    mkpath(output_dir)
    forest_output = joinpath(output_dir, "forest_$(case_name).png")
    plot_forest(mesh; constrained_faces=constrained_faces, 
                potential_fixed_edges=potential_fixed_edges,
                savepath=forest_output)
    println("Forest visualization saved to: $forest_output")
    
    println("\n" * "="^60)
    println("PROBLEM SETUP")
    println("="^60)
    
    # Build the system matrix A and vector b using assemble_energy_system
    println("\n" * "="^60)
    println("BUILDING SYSTEM MATRIX Ax = b")
    println("="^60)
    
    A, b, var_to_p, var_to_theta, n_free_p, n_free_theta = assemble_energy_system(
        field, topology, fixed_edges_per_face, constrained_faces, constrained_angles
    )
    
    n_p_vars = length(var_to_p)
    n_fixed_p = n_p_vars - n_free_p
    n_theta_vars = length(var_to_theta)
    n_fixed_theta = n_theta_vars - n_free_theta
    n_vars = size(A, 2)
    n_equations = size(A, 1)
    
    println("\nVariable counts:")
    println("  Free p variables: $n_free_p")
    println("  Fixed p variables: $n_fixed_p")
    println("  Total p variables: $n_p_vars")
    println("  Free θ variables: $n_free_theta")
    println("  Fixed θ variables: $n_fixed_theta")
    println("  Total θ variables: $n_theta_vars")
    println("  Fixed edges: $(length(fixed_edges_per_face))")
    println("  Total variables: $n_vars")
    println("  Total equations: $n_equations")
    
    println("\nMatrix A dimensions: $(size(A))")
    println("Vector b dimensions: $(size(b))")
    println("\nMatrix A (full):")
    display(Matrix(A))
    println("\n\nVector b ($(length(b))):")
    display(rad2deg.(b)) 
    
    # Solve using greedy mixed integer solver
    println("\n" * "="^60)
    println("SOLVING WITH GREEDY MIXED INTEGER SOLVER")
    println("="^60)
    
    # Only the FREE p variables need to be rounded (fixed ones are already constrained to 0)
    k = n_free_p
    println("\nSolving for $k FREE integer period jump variables...")
    println("  (Note: $n_fixed_p period jumps are already fixed to 0)")
    println("  Plus $n_free_theta FREE continuous angle variables")
    println("  (Note: $n_fixed_theta angles are constrained to user-specified values)")
    
    # Create directory for animation frames
    anim_dir = joinpath(output_dir, "animation_frames")
    mkpath(anim_dir)
    
    # Callback to save frame at each step
    frame_paths = String[]
    function save_frame_callback(x, step, total_steps)
        # Extract angles from solution using the mapping
        angles_step = zeros(Float64, n_theta_vars)
        for (var_idx, face_idx) in var_to_theta
            angles_step[face_idx] = x[var_idx]
        end
        
        # Extract period jumps from solution
        p_jumps = Dict{Tuple{Int,Int}, Float64}()
        for (var_idx, (face_i, face_j)) in var_to_p
            p_jumps[(face_i, face_j)] = x[var_idx]
        end
        
        # Save frame
        frame_path = joinpath(anim_dir, "frame_$(lpad(step, 4, '0')).png")
        push!(frame_paths, frame_path)
        
        plot_frame_field(mesh, angles_step; 
                        savepath=frame_path,
                        constrained_faces=constrained_faces,
                        arrow_scale=0.4,
                        title="Solving: Step $step / $total_steps",
                        period_jumps=p_jumps,
                        topology=topology)
    end
    
    x_solution = greedy_mixed_integer_solve(A, b, k; 
        verbose=true,
        tau=1e-8,
        max_gauss_seidel_iters = 100000,
        fallback_to_full_solve=true,
        use_queue=true,  # Use queue-based Gauss-Seidel (CrossField.jl style)
        mod4=true,        # Apply mod 4 constraint to keep period jumps in [0,3]
        callback=save_frame_callback
    )
    
    println("\nSolution computed!")
    println("\nPeriod jumps (p):")
    p_values = x_solution[1:n_p_vars]
    for (var_idx, edge_pair) in sort(collect(var_to_p), by=x->x[1])
        val = x_solution[var_idx]
        is_fixed = var_idx > n_free_p
        status = is_fixed ? " [FIXED]" : " [FREE]"
        println("  p_$(edge_pair): $(Int(round(val)))$status")
    end
    
    println("\nAngles (θ) in degrees:")
    theta_values = x_solution[n_p_vars+1:end]
    for (var_idx, face_idx) in sort(collect(var_to_theta), by=x->x[1])
        val = x_solution[var_idx]
        println("  θ_$face_idx: $(round(rad2deg(val), digits=2))°")
    end
    
    # Verify solution
    println("\nVerifying solution...")
    residual = A * x_solution - b
    max_residual = maximum(abs.(residual))
    println("  Maximum residual: $max_residual")
    
    # Check that period jumps are integers
    p_int_error = maximum(abs.(p_values - round.(p_values)))
    println("  Maximum period jump rounding error: $p_int_error")
    
    # Extract angle values and period jumps for each face
    angles = zeros(Float64, n_theta_vars)
    for (var_idx, face_idx) in var_to_theta
        angles[face_idx] = x_solution[var_idx]
    end
    
    # Update field with solution
    field.theta = angles
    for (var_idx, (face_i, face_j)) in var_to_p
        p_val = Int(round(x_solution[var_idx]))
        set_period_jump!(field, face_i, face_j, p_val)
    end
    
    # Compute final smoothness energy
    println("\nComputing smoothness energy...")
    final_energy = compute_smoothness_energy(field)
    println("  Final smoothness energy: $final_energy")
    
    # Print individual edge contributions
    # println("\nEdge-by-edge energy contributions:")
    # for face_i in 1:topology.n_faces
    #     for (face_j, edge) in topology.dual_adj[face_i]
    #         if face_i < face_j
    #             kappa_ij = get_transport_angle(field, face_i, face_j)
    #             p_ij = get_period_jump(field, face_i, face_j)
    #             diff = field.theta[face_i] + kappa_ij + (π/2) * p_ij - field.theta[face_j]
    #             edge_energy = diff^2
    #             println("  Edge ($face_i, $face_j): θ_$face_i=$(round(rad2deg(field.theta[face_i]), digits=2))°, " *
    #                     "κ=$(round(rad2deg(kappa_ij), digits=2))°, p=$p_ij, " *
    #                     "θ_$face_j=$(round(rad2deg(field.theta[face_j]), digits=2))° → " *
    #                     "diff=$(round(rad2deg(diff), digits=2))°, energy=$(round(edge_energy, digits=6))")
    #         end
    #     end
    # end
    
    # Plot frame field
    println("\nPlotting frame field...")
    field_output = joinpath(output_dir, "frame_field_$(case_name).png")
    plot_frame_field(mesh, angles; savepath=field_output, 
                     constrained_faces=constrained_faces,
                     arrow_scale=0.4)
    println("Frame field visualization saved to: $field_output")
    
    # Create GIF from frames
    if length(frame_paths) > 0
        println("\nCreating animation GIF...")
        gif_path = joinpath(output_dir, "frame_field_evolution_$(case_name).gif")
        
        # Load all frames
        frames = [load(path) for path in frame_paths]
        
        # Save as GIF with 2 fps (500ms per frame)
        save(gif_path, cat(frames..., dims=3); fps=2)
        println("Animation saved to: $gif_path")
        
        # Clean up frame files
        println("Cleaning up temporary frames...")
        for path in frame_paths
            rm(path, force=true)
        end
        try
            rm(anim_dir, force=true, recursive=true)
        catch
            # Directory might not be empty if other files exist
        end
    end
    
    # Create visualization
    println("\n" * "="^60)
    println("CREATING VISUALIZATIONS")
    println("="^60)
    
    fig = Figure(size=(1600, 1200))
    
    # Matrix A structure
    ax1 = Axis(fig[1, 1], aspect=DataAspect(), 
               title="Matrix A Structure ($(n_equations)×$(n_vars))",
               xlabel="Variables", ylabel="Equations")
    
    # Plot non-zero entries
    rows, cols, vals = findnz(A)
    colors = [v > 0 ? :blue : :red for v in vals]
    scatter!(ax1, cols, rows, color=colors, markersize=20)
    
    # Add grid
    hlines!(ax1, 1:n_equations, color=:gray, linewidth=0.5, alpha=0.3)
    vlines!(ax1, 1:n_vars, color=:gray, linewidth=0.5, alpha=0.3)
    
    # Add variable labels
    var_labels = ["p_$i" for i in 1:n_p_vars] ∪ ["θ_$i" for i in 1:n_theta_vars]
    ax1.xticks = (1:n_vars, var_labels)
    ax1.xticklabelrotation = π/4
    
    # Invert y-axis so equation 1 is at top
    ylims!(ax1, n_equations + 0.5, 0.5)
    
    # Matrix A with values
    ax2 = Axis(fig[2, 1], aspect=DataAspect(),
               title="Matrix A with Values",
               xlabel="Variables", ylabel="Equations")
    
    heatmap!(ax2, 1:n_vars, 1:n_equations, Matrix(A), colormap=:RdBu)
    
    # Add text annotations for values
    for i in 1:n_equations
        for j in 1:n_vars
            val = A[i, j]
            if abs(val) > 1e-5
                val_str = abs(val - π/2) < 1e-5 ? "π/2" : 
                         abs(val - 1.0) < 1e-5 ? "1" :
                         abs(val + 1.0) < 1e-5 ? "-1" : 
                         abs(val - 2.0) < 1e-5 ? "2" :
                         abs(val + 2.0) < 1e-5 ? "-2" :
                         abs(val - π) < 1e-5 ? "π" :
                         abs(val + π) < 1e-5 ? "-π" :
                         "$(round(val, digits=2))"
                text!(ax2, j, i, text=val_str, color=:white, 
                      fontsize=10, align=(:center, :center))
            end
        end
    end
    
    ax2.xticks = (1:n_vars, var_labels)
    ax2.xticklabelrotation = π/4
    ylims!(ax2, n_equations + 0.5, 0.5)
    
    # Save figure
    output_file = joinpath(output_dir, "matrix_$(case_name).png")
    save(output_file, fig)
    println("\nMatrix visualization saved to: $output_file")
    
    println("\n" * "="^60)
end

main()
