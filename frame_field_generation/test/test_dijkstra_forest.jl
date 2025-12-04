"""
    test_dijkstra_forest.jl

Test script to verify Dijkstra forest construction and edge fixing logic.
"""

using GeometryBasics

# Add src to load path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using FrameField

"""
    generate_simple_mesh()

Generate a very simple mesh with 2 triangles sharing an edge.
"""
function generate_simple_mesh()
    # 4 vertices forming a square
    vs = [
        Point3f(0, 0, 0),
        Point3f(1, 0, 0),
        Point3f(1, 1, 0),
        Point3f(0, 1, 0)
    ]
    
    # 2 triangles
    fs = [
        TriangleFace(1, 2, 3),  # Bottom-right triangle
        TriangleFace(1, 3, 4)   # Top-left triangle
    ]
    
    return Mesh(vs, fs)
end

"""
    print_edge_info(mesh, topology, forest)

Print detailed information about edges and their status.
"""
function print_edge_info(mesh, topology, forest)
    println("\n" * "="^60)
    println("Edge Information")
    println("="^60)
    
    # Print all edges
    println("\nAll Edges:")
    for (edge_tuple, edge_idx) in topology.edge_map
        v1, v2 = edge_tuple
        is_tree = edge_idx in forest.fixable_edges
        is_fixed = edge_idx in values(forest.fixed_edges_per_face)
        is_free = edge_idx in forest.free_edges
        
        status = []
        is_tree && push!(status, "TREE")
        is_fixed && push!(status, "FIXED")
        is_free && push!(status, "FREE")
        
        println("  Edge $edge_idx: ($v1, $v2) - $(join(status, ", "))")
    end
    
    # Print face information
    fs = faces(mesh)
    println("\nFace Information:")
    for face_idx in 1:length(fs)
        face = fs[face_idx]
        
        # Get all three edges of this face
        edges_of_face = [
            (min(face[1], face[2]), max(face[1], face[2])),
            (min(face[2], face[3]), max(face[2], face[3])),
            (min(face[3], face[1]), max(face[3], face[1]))
        ]
        edge_indices = [topology.edge_map[e] for e in edges_of_face]
        
        # Which edge is fixed for this face?
        fixed_edge = forest.fixed_edges_per_face[face_idx]
        
        println("\n  Face $face_idx: vertices $(face[1]), $(face[2]), $(face[3])")
        println("    Edges: $edge_indices")
        println("    Fixed edge: $fixed_edge")
        
        # Verify fixed edge is one of the face's edges
        if fixed_edge in edge_indices
            println("    ✓ Fixed edge belongs to this face")
        else
            println("    ✗ ERROR: Fixed edge does not belong to this face!")
        end
    end
end

"""
    verify_constraints(mesh, topology, forest)

Verify that the forest satisfies required constraints.
"""
function verify_constraints(mesh, topology, forest)
    println("\n" * "="^60)
    println("Constraint Verification")
    println("="^60)
    
    fs = faces(mesh)
    n_faces = length(fs)
    
    # Constraint 1: Each face should have exactly one fixed edge
    println("\n1. Each face has one fixed edge:")
    all_good = true
    for face_idx in 1:n_faces
        if haskey(forest.fixed_edges_per_face, face_idx)
            println("   ✓ Face $face_idx: edge $(forest.fixed_edges_per_face[face_idx])")
        else
            println("   ✗ Face $face_idx: NO FIXED EDGE")
            all_good = false
        end
    end
    println(all_good ? "   PASS" : "   FAIL")
    
    # Constraint 2: Fixed edge must be one of the face's edges
    println("\n2. Fixed edge belongs to its face:")
    all_good = true
    for face_idx in 1:n_faces
        face = fs[face_idx]
        edges_of_face = [
            (min(face[1], face[2]), max(face[1], face[2])),
            (min(face[2], face[3]), max(face[2], face[3])),
            (min(face[3], face[1]), max(face[3], face[1]))
        ]
        edge_indices = [topology.edge_map[e] for e in edges_of_face]
        
        fixed_edge = forest.fixed_edges_per_face[face_idx]
        if fixed_edge in edge_indices
            println("   ✓ Face $face_idx: edge $fixed_edge is valid")
        else
            println("   ✗ Face $face_idx: edge $fixed_edge NOT in face edges $edge_indices")
            all_good = false
        end
    end
    println(all_good ? "   PASS" : "   FAIL")
    
    # Constraint 3: Free edges + fixed edges = all edges
    println("\n3. Edge accounting (free + fixed = total):")
    n_fixed = length(Set(values(forest.fixed_edges_per_face)))
    n_free = length(forest.free_edges)
    n_total = topology.n_edges
    
    println("   Fixed edges: $n_fixed")
    println("   Free edges: $n_free")
    println("   Total edges: $n_total")
    println("   Sum: $(n_fixed + n_free)")
    
    if n_fixed + n_free == n_total
        println("   ✓ PASS: All edges accounted for")
    else
        println("   ✗ FAIL: Edge count mismatch")
        
        # Debug: show which edges are missing
        all_fixed = Set(values(forest.fixed_edges_per_face))
        all_free = Set(forest.free_edges)
        all_edges = Set(1:n_total)
        
        fixed_and_free = union(all_fixed, all_free)
        missing = setdiff(all_edges, fixed_and_free)
        duplicate = intersect(all_fixed, all_free)
        
        if !isempty(missing)
            println("   Missing edges: $missing")
        end
        if !isempty(duplicate)
            println("   Duplicate edges (both fixed and free): $duplicate")
        end
    end
    
    # Constraint 4: For a mesh with F faces, we should fix F edges
    println("\n4. Number of fixed edges equals number of faces:")
    n_faces_total = length(fs)
    n_fixed_edges = length(forest.fixed_edges_per_face)
    
    println("   Faces: $n_faces_total")
    println("   Fixed edge assignments: $n_fixed_edges")
    
    if n_fixed_edges == n_faces_total
        println("   ✓ PASS: One fixed edge per face")
    else
        println("   ✗ FAIL: Mismatch in fixed edge count")
    end
end

println("="^60)
println("Testing Dijkstra Forest Construction")
println("="^60)

# Test 1: Simple 2-triangle mesh
println("\n" * "="^60)
println("TEST 1: Simple 2-Triangle Mesh")
println("="^60)

mesh = generate_simple_mesh()
println("\nMesh: $(length(coordinates(mesh))) vertices, $(length(faces(mesh))) faces")

topology = build_mesh_topology(mesh)
println("Topology: $(topology.n_edges) edges")

forest = build_dijkstra_forest(mesh, topology, Int[])
println("Forest built successfully")

print_edge_info(mesh, topology, forest)
verify_constraints(mesh, topology, forest)

# Test 2: Larger mesh
println("\n\n" * "="^60)
println("TEST 2: 3×3 Square Mesh")
println("="^60)

function generate_square_mesh(; nx=3, ny=3)
    vs = Point3f[]
    for j in 0:ny
        for i in 0:nx
            push!(vs, Point3f(i/nx, j/ny, 0))
        end
    end
    
    fs = TriangleFace{Int}[]
    row_stride = nx + 1
    
    for j in 1:ny
        for i in 1:nx
            v_bl = (j-1) * row_stride + i
            v_br = v_bl + 1
            v_tl = v_bl + row_stride
            v_tr = v_tl + 1
            
            push!(fs, TriangleFace(v_bl, v_br, v_tr))
            push!(fs, TriangleFace(v_bl, v_tr, v_tl))
        end
    end
    
    return Mesh(vs, fs)
end

mesh2 = generate_square_mesh(nx=3, ny=3)
println("\nMesh: $(length(coordinates(mesh2))) vertices, $(length(faces(mesh2))) faces")

topology2 = build_mesh_topology(mesh2)
println("Topology: $(topology2.n_edges) edges")

forest2 = build_dijkstra_forest(mesh2, topology2, Int[])
println("Forest built successfully")

verify_constraints(mesh2, topology2, forest2)

println("\n" * "="^60)
println("Summary")
println("="^60)
println("Both tests completed. Check output above for any FAIL messages.")
println("="^60)
