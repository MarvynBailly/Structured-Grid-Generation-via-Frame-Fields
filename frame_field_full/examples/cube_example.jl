using FrameFieldFull
using FrameFieldFull.Types
using FrameFieldFull.Topology
using FrameFieldFull.Solver
using FrameFieldFull.Analysis
using FrameFieldFull.MeshIO
using FrameFieldFull.Cutting
using LinearAlgebra

include("visualize_cross_field.jl")



function create_fine_cube_mesh(filename; n=6)
    # n: Number of subdivisions per edge (e.g., n=1 gives the original 2x2 per face)
    # Total triangles = 6 faces * 2 * n^2
    
    vertices = Vector{Tuple{Float64, Float64, Float64}}()
    # Dict to map coordinates to vertex index (for welding edges)
    vert_map = Dict{Tuple{Float64, Float64, Float64}, Int}()
    faces = Vector{Tuple{Int, Int, Int}}()

    # Helper: Gets or creates vertex index
    function get_vert_idx(x, y, z)
        # Round to avoid floating point issues at shared edges
        p_snap = (round(x, digits=5), round(y, digits=5), round(z, digits=5))
        if haskey(vert_map, p_snap)
            return vert_map[p_snap]
        else
            push!(vertices, p_snap)
            idx = length(vertices) - 1 # SU2 uses 0-based indexing
            vert_map[p_snap] = idx
            return idx
        end
    end

    # Helper: Generates grid of triangles for a face defined by origin and two axes
    function make_face(origin, u_vec, v_vec)
        for i in 0:n-1
            for j in 0:n-1
                # Compute 4 corners of the quad cell
                p0 = origin .+ (i / n) .* u_vec .+ (j / n) .* v_vec
                p1 = origin .+ ((i+1) / n) .* u_vec .+ (j / n) .* v_vec
                p2 = origin .+ ((i+1) / n) .* u_vec .+ ((j+1) / n) .* v_vec
                p3 = origin .+ (i / n) .* u_vec .+ ((j+1) / n) .* v_vec

                idx0 = get_vert_idx(p0...)
                idx1 = get_vert_idx(p1...)
                idx2 = get_vert_idx(p2...)
                idx3 = get_vert_idx(p3...)

                # Split Quad into 2 Triangles (CCW winding)
                push!(faces, (idx0, idx1, idx2))
                push!(faces, (idx0, idx2, idx3))
            end
        end
    end

    # --- Generate 6 Faces with correct outward normals ---
    # Origin is usually bottom-left corner of the face relative to its normal
    
    # 1. Front Face (z = 1)
    make_face((-1.0, -1.0, 1.0), (2.0, 0.0, 0.0), (0.0, 2.0, 0.0))
    
    # 2. Back Face (z = -1) - Look from back, x goes right-to-left
    make_face((1.0, -1.0, -1.0), (-2.0, 0.0, 0.0), (0.0, 2.0, 0.0))
    
    # 3. Right Face (x = 1) - y is up, z is back
    make_face((1.0, -1.0, 1.0), (0.0, 0.0, -2.0), (0.0, 2.0, 0.0))
    
    # 4. Left Face (x = -1) - y is up, z is forward
    make_face((-1.0, -1.0, -1.0), (0.0, 0.0, 2.0), (0.0, 2.0, 0.0))
    
    # 5. Top Face (y = 1) - x is right, z is back
    make_face((-1.0, 1.0, 1.0), (2.0, 0.0, 0.0), (0.0, 0.0, -2.0))
    
    # 6. Bottom Face (y = -1) - x is right, z is forward
    make_face((-1.0, -1.0, -1.0), (2.0, 0.0, 0.0), (0.0, 0.0, 2.0))

    # --- Write SU2 File ---
    open(filename, "w") do io
        write(io, "NDIME= 3\n")
        write(io, "NPOIN= $(length(vertices))\n")
        for v in vertices
            write(io, "$(v[1]) $(v[2]) $(v[3])\n")
        end
        write(io, "NELEM= $(length(faces))\n")
        for f in faces
            # SU2 Type 5 = Triangle
            write(io, "5 $(f[1]) $(f[2]) $(f[3])\n")
        end
    end
end



# --- 3. Run Pipeline ---
function main()
    @info("Generating mesh...")
    mesh_file = "test_cube.su2"
    n = 15
    create_fine_cube_mesh(mesh_file; n=n)
    save_file_path = "output/cube$(n)_cross_field.png"

    
    @info("Reading mesh...")
    verts, faces = read_mesh(mesh_file)
    topo = build_topology(verts, faces)
    
    @info("Setting constraints...")
    constraints = Dict{Int, Float64}() # No constraints for closed surface
    
    @info("Initializing field...")
    field = initialize_field(topo, constraints)
    
    @info("Solving...")
    solve_greedy!(field; verbose=false)
    
    @info("Optimizing Singularities...")
    optimize_singularities!(field; verbose=false)
    
    # Check Topology
    sings = compute_singularities(field)
    if isempty(sings)
        @info("Total Singularity Index: 0.0")
    else
        total_idx = sum(s[2] for s in sings)
        @info("Total Singularity Index: $total_idx (Expected: 2.0)")
    end

    # cut the mesh M into a disk topology
    cuts = compute_cut_graph(topo, sings)


    @info("Visualizing...")
    visualize_field(field, title="Cube Flow (Chi=2)", final_cuts=cuts, save_path=save_file_path)
end

main()