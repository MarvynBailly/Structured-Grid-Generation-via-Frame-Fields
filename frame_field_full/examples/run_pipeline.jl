using FrameFieldFull
using Printf


include("visualize_cross_field.jl")

function main()
    # --- Configuration ---
    # filename = "../../triangulations/disk-radial-fine.msh"
    filename = "../triangulations/mesh_airfoil_dae11.su2"
    output_dir = "output"
    
    verbose = true
    animate = false
    frame_interval = 5

    # --- Setup ---
    if !isfile(filename)
        println("File not found: $filename")
        return
    end
    mkpath(output_dir)

    # --- Pipeline ---
    println("Reading mesh...")
    verts, faces = read_mesh(filename)
    topo = build_topology(verts, faces)

    println("Detecting boundaries...")
    constraints = compute_boundary_constraints(topo)
    # constraints = Dict{Int, Float64}(1 => 0.0)  # Placeholder: Constrain face 1 to 0.0 radians
    
    field = initialize_field(topo, constraints)

    # --- Animation Callback ---
    frame_count = 0
    function animation_callback(field, count, num_p, error)
        if animate
            if count == 0
                # Initial frame
                frame_path = joinpath(output_dir, "frame_0000.png")
                plot_frame(field, frame_path, "Initial State (Continuous Solution)")
                println("  Saved initial frame")
                frame_count += 1
            elseif count % frame_interval == 0 || count == num_p
                # Intermediate/Final frames
                frame_path = joinpath(output_dir, @sprintf("frame_%04d.png", frame_count))
                progress_pct = round(100 * count / num_p, digits=1)
                title = @sprintf("Greedy Rounding: %d/%d (%.1f%%) - Error: %.4f", count, num_p, progress_pct, error)
                plot_frame(field, frame_path, title)
                
                if verbose
                    println("  Saved frame $frame_count (iteration $count)")
                end
                frame_count += 1
            end
        end
    end

    # --- Solve ---
    println("Solving cross field...")
    solve_greedy!(field; verbose=verbose, callback=animation_callback)
    println("Cross field computed successfully.")

    # --- Analysis & Extraction ---
    println("Computing singularities...")
    sings = compute_singularities(field; verbose=verbose)
    println("Found $(length(sings)) singularities.")

    # println("Optimizing singularity positions...")
    # optimize_singularities!(field; verbose=verbose)


    # cut_edges = compute_cut_graph(topo)
    
    # face_rotations, _ = propagate_orientations(field, cut_edges)
    
    # println("Plotting globally smooth field...")
    # plot_smooth_global_field(field, face_rotations, joinpath(output_dir, "global_smooth_field.png"); 
    #                          cut_edges=cut_edges, verbose=true)
    
    # u_coords, v_coords = compute_parameterization_least_squares(field, face_rotations, cut_edges)
    
    # quads, _ = extract_quad_mesh(topo, u_coords, v_coords; verbose=true)

    # println("Plotting quad mesh...")
    # plot_quad_mesh(topo, u_coords, v_coords, quads, joinpath(output_dir, "quad_mesh_result.png"); verbose=true)
    
    # println("Generating final visualization...")
    # plot_results(field, joinpath(output_dir, "miq_solution_full_nopt.png"); verbose=verbose, cut_edges=nothing, show_period_jumps = true)
    println("Visualizing...")
    visualize_field(field)

    # println("\nComplete! All outputs saved to the '$output_dir' directory.")
end

main()
