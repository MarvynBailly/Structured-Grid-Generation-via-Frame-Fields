"""
    test_frame_field.jl

Unit tests for frame field generation module.
"""

using Test
using LinearAlgebra
using GeometryBasics

# Add parent directory to load path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using FrameField

@testset "FrameField Module Tests" begin
    
    @testset "Mesh Topology Construction" begin
        # Create simple square mesh
        vs = [
            Point3f(0, 0, 0), Point3f(1, 0, 0),
            Point3f(1, 1, 0), Point3f(0, 1, 0)
        ]
        fs = [
            TriangleFace(1, 2, 3),
            TriangleFace(1, 3, 4)
        ]
        mesh = Mesh(vs, fs)
        
        topology = build_mesh_topology(mesh)
        
        @test topology.n_edges == 5  # 4 boundary + 1 internal
        @test length(topology.dual_adj) == 2  # 2 faces
        @test length(topology.dual_adj[1]) == 1  # Face 1 has 1 neighbor
        @test length(topology.dual_adj[2]) == 1  # Face 2 has 1 neighbor
        @test length(topology.reference_edges) == 2
        
        # Check edge map
        @test length(topology.edge_map) == 5
        
        println("✓ Mesh topology construction test passed")
    end
    
    @testset "Dijkstra Forest" begin
        # Simple mesh
        vs = [
            Point3f(0, 0, 0), Point3f(1, 0, 0), Point3f(2, 0, 0),
            Point3f(0, 1, 0), Point3f(1, 1, 0), Point3f(2, 1, 0)
        ]
        fs = [
            TriangleFace(1, 2, 4),
            TriangleFace(2, 5, 4),
            TriangleFace(2, 3, 5),
            TriangleFace(3, 6, 5)
        ]
        mesh = Mesh(vs, fs)
        topology = build_mesh_topology(mesh)
        
        # Build forest from face 1
        forest = build_dijkstra_forest(mesh, topology, [1])
        
        @test !isempty(forest.tree_edges)
        @test length(forest.tree_edges) <= topology.n_edges
        @test !isempty(forest.parent)
        
        # All faces should be reachable
        @test length(forest.parent) == length(fs)
        
        println("✓ Dijkstra forest test passed")
    end
    
    @testset "Energy Assembly" begin
        # Simple 2-triangle mesh
        vs = [
            Point3f(0, 0, 0), Point3f(1, 0, 0),
            Point3f(1, 1, 0), Point3f(0, 1, 0)
        ]
        fs = [
            TriangleFace(1, 2, 3),
            TriangleFace(1, 3, 4)
        ]
        mesh = Mesh(vs, fs)
        topology = build_mesh_topology(mesh)
        forest = build_dijkstra_forest(mesh, topology, Int[])
        
        A, b, var_to_p, var_to_theta = assemble_system_matrix(
            mesh, topology, forest, Int[], Float64[]
        )
        
        @test size(A, 1) == size(A, 2)  # Square matrix
        @test length(b) == size(A, 1)
        @test !isempty(var_to_theta) || !isempty(var_to_p)
        
        println("✓ Energy assembly test passed")
    end
    
    @testset "Simple Frame Field Generation" begin
        # Create small square mesh
        vs = [
            Point3f(0, 0, 0), Point3f(0.5, 0, 0), Point3f(1, 0, 0),
            Point3f(0, 0.5, 0), Point3f(0.5, 0.5, 0), Point3f(1, 0.5, 0),
            Point3f(0, 1, 0), Point3f(0.5, 1, 0), Point3f(1, 1, 0)
        ]
        
        # 8 triangles forming a 2x2 quad grid
        fs = [
            TriangleFace(1, 2, 4), TriangleFace(2, 5, 4),
            TriangleFace(2, 3, 5), TriangleFace(3, 6, 5),
            TriangleFace(4, 5, 7), TriangleFace(5, 8, 7),
            TriangleFace(5, 6, 8), TriangleFace(6, 9, 8)
        ]
        mesh = Mesh(vs, fs)
        
        # Generate frame field without constraints
        result = generate_frame_field(mesh, verbose=false)
        
        @test length(result.angles) == length(fs)
        @test !isempty(result.period_jumps)
        @test result.energy >= 0.0
        @test isa(result.singularities, Vector)
        
        println("✓ Simple frame field generation test passed")
        println("  Energy: $(result.energy)")
        println("  Singularities: $(length(result.singularities))")
    end
    
    @testset "Constrained Frame Field" begin
        # Simple mesh with boundary constraint
        vs = [
            Point3f(0, 0, 0), Point3f(1, 0, 0),
            Point3f(1, 1, 0), Point3f(0, 1, 0)
        ]
        fs = [
            TriangleFace(1, 2, 3),
            TriangleFace(1, 3, 4)
        ]
        mesh = Mesh(vs, fs)
        
        # Constrain first face to angle 0
        result = generate_frame_field(
            mesh,
            constrained_faces=[1],
            constrained_angles=[0.0],
            verbose=false
        )
        
        @test abs(result.angles[1] - 0.0) < 1e-10
        @test result.energy >= 0.0
        
        println("✓ Constrained frame field test passed")
    end
    
    @testset "Frame Direction Extraction" begin
        # Create simple result
        angles = [0.0, π/4, π/2]
        period_jumps = Dict(1 => 0, 2 => 1, 3 => -1)
        singularities = Tuple{Int, Float64}[]
        result = FrameFieldResult(angles, period_jumps, singularities, 0.0, true, 10)
        
        dirs = get_frame_directions(result, 1)
        @test length(dirs) == 4
        @test dirs[1] ≈ 0.0
        @test dirs[2] ≈ π/2
        @test dirs[3] ≈ π
        @test dirs[4] ≈ 3π/2
        
        println("✓ Frame direction extraction test passed")
    end
end

println("\n" * "="^60)
println("All Frame Field Tests Passed!")
println("="^60)
