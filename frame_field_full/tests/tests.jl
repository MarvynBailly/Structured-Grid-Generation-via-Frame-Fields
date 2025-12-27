using FrameFieldFull
using Test
using LinearAlgebra
using Printf
include("../examples/visualize_cross_field.jl")

# --- Helper: Mesh Generators (SU2 Format) ---

function create_3x3_grid_mesh(filename)
    open(filename, "w") do io
        write(io, """
        NDIME= 2
        NPOIN= 9
        -1 -1
        0 -1
        1 -1
        -1 0
        0 0
        1 0
        -1 1
        0 1
        1 1
        NELEM= 8
        5 0 1 4
        5 0 4 3
        5 1 2 5
        5 1 5 4
        5 3 4 7
        5 3 7 6
        5 4 5 8
        5 4 8 7
        """)
    end
end

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

function create_cone_mesh(filename)
    # A pyramid/cone. 1 center vertex (top), 4 base vertices.
    # Center vertex defect: Sum of angles < 360.
    # Base: (-1,-1), (1,-1), (1,1), (-1,1) at z=0. Top: (0,0,1).
    # Triangle angles at top: 4 isosceles triangles.
    # If height is high, angle is small -> High Defect -> Positive Index.
    open(filename, "w") do io
        write(io, """
        NDIME= 3
        NPOIN= 5
        -1 -1 0
         1 -1 0
         1  1 0
        -1  1 0
         0  0 1
        NELEM= 4
        5 0 1 4
        5 1 2 4
        5 2 3 4
        5 3 0 4
        """)
    end
end


# --- The Test Suite ---

@testset "Verification Tests" begin

    @testset "1. The Perfect Square (Smoothness)" begin
        # Goal: Ensure interior vertex has Index 0.0 when boundary is aligned
        mesh_file = "test_square.su2"
        create_3x3_grid_mesh(mesh_file)
        
        verts, faces = read_mesh(mesh_file)
        topo = build_topology(verts, faces)
        
        # Constrain boundary faces to align with X-axis (0.0 radians)
        # In this grid, faces 1,2 (bottom) and 7,8 (top) etc are boundaries.
        # We'll just constrain Face 1 (bottom left) to 0.0.
        # The solver should propagate this smoothness.
        constraints = Dict(1 => 0.0)
        
        field = initialize_field(topo, constraints)
        solve_greedy!(field; verbose=false)
        
        sings = compute_singularities(field; verbose=false)
        
        # Vertex 5 (index 5) is the only interior vertex.
        # We expect NO singularities in the interior.
        interior_sings = filter(s -> s[1] == 5, sings)
        
        @test isempty(interior_sings)

        
        # Clean up
        rm(mesh_file)
    end

    @testset "2. The Sphere/Cube (Topology Check)" begin
        # Goal: Ensure Sum of Indices = 2.0 (Euler Characteristic)
        # And specifically, a Cube often creates 8 singularities of index 1/4.
        mesh_file = "test_cube.su2"
        create_fine_cube_mesh(mesh_file; n=10)
        
        verts, faces = read_mesh(mesh_file)
        topo = build_topology(verts, faces)
        
        # No boundary constraints for a closed surface
        constraints = Dict{Int, Float64}()
        
        field = initialize_field(topo, constraints)
        solve_greedy!(field; verbose=false)
        optimize_singularities!(field; verbose=true)
        
        sings = compute_singularities(field; verbose=false)
        
        sum_indices = sum(s[2] for s in sings)
        println("  Cube Index Sum: $sum_indices (Expected: 2.0)")
        for (sing) in sings
            println(@sprintf("    Vertex %d: Index %.5f", sing[1], sing[2]))
        end
        
        @test isapprox(sum_indices, 2.0, atol=1e-5)

                # visualize
        println("Visualizing...")
        visualize_field(field)

        
        # For a cube, we typically expect indices of +0.25 at the 8 corners.
        # (Depending on triangulation, some might merge, but sum must be 2).
        @test !isempty(sings)
        
        rm(mesh_file)
    end
    
    # @testset "3. The Cone (Geometric Singularity)" begin
    #     # Goal: Verify that Index comes from Geometry (Gaussian Curvature)
    #     mesh_file = "test_cone.su2"
    #     create_cone_mesh(mesh_file)
        
    #     verts, faces = read_mesh(mesh_file)
    #     topo = build_topology(verts, faces)
        
    #     # Constrain the base face (Face 1) to 0.0
    #     constraints = Dict(1 => 0.0)
        
    #     field = initialize_field(topo, constraints)
    #     solve_greedy!(field; verbose=false)
        
    #     sings = compute_singularities(field; verbose=false)
        
    #     # Find singularity at the top vertex (Vertex 5, index 5)
    #     top_sing = filter(s -> s[1] == 5, sings)
        
    #     if !isempty(top_sing)
    #         idx = top_sing[1][2]
    #         println("  Cone Top Index: $idx")
    #         # A cone tip has positive Gaussian curvature -> Positive Index.
    #         @test idx > 0.0 
    #     end
        
    #     rm(mesh_file)
    # end

end