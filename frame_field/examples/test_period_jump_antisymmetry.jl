"""
    test_period_jump_antisymmetry.jl

Test that period jumps are stored with antisymmetric property (p_ji = -p_ij).
"""

using Pkg
Pkg.activate(".")

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using GeometryBasics

include("../src/meshio.jl")
include("../src/mesh_types.jl")

function main()
    println("Loading simplest-square mesh...")
    mesh_file = joinpath(@__DIR__, "..", "..", "triangulations", "simplest-square.msh")
    mesh = load_triangulation(mesh_file)
    
    println("Building mesh topology...")
    topology = build_mesh_topology(mesh)
    
    println("Creating frame field...")
    field = CrossField(topology)
    
    println("\n" * "="^60)
    println("TESTING PERIOD JUMP ANTISYMMETRY")
    println("="^60)
    
    # Test 1: Verify initial state (all zeros)
    println("\nTest 1: Initial state (all zeros)")
    println("-"^60)
    
    face_i = 1
    face_j = 4
    
    p_ij = get_period_jump(field, face_i, face_j)
    p_ji = get_period_jump(field, face_j, face_i)
    
    println("p_($face_i,$face_j) = $p_ij")
    println("p_($face_j,$face_i) = $p_ji")
    println("Initial values both zero: ", p_ij == 0 && p_ji == 0)
    
    # Test 2: Set a period jump and verify antisymmetry
    println("\nTest 2: Set p_ij = 2 and verify p_ji = -2")
    println("-"^60)
    
    set_period_jump!(field, face_i, face_j, 2)
    
    p_ij = get_period_jump(field, face_i, face_j)
    p_ji = get_period_jump(field, face_j, face_i)
    
    println("p_($face_i,$face_j) = $p_ij")
    println("p_($face_j,$face_i) = $p_ji")
    println("Antisymmetry satisfied (p_ji = -p_ij): ", p_ji == -p_ij)
    
    # Test 3: Verify energy computation uses the correct sign
    println("\nTest 3: Energy computation with period jump")
    println("-"^60)
    
    # Set theta values
    field.theta[face_i] = 0.0
    field.theta[face_j] = π/4
    
    kappa_ij = get_transport_angle(field, face_i, face_j)
    
    # E_ij = (θ_i + κ_ij + (π/2)p_ij - θ_j)²
    diff_ij = field.theta[face_i] + kappa_ij + (π/2) * p_ij - field.theta[face_j]
    energy_ij = diff_ij^2
    
    println("θ_$face_i = $(field.theta[face_i])")
    println("θ_$face_j = $(field.theta[face_j])")
    println("κ_($face_i,$face_j) = $kappa_ij rad = $(rad2deg(kappa_ij))°")
    println("p_($face_i,$face_j) = $p_ij")
    println("\nEnergy term E_ij = (θ_i + κ_ij + (π/2)p_ij - θ_j)²")
    println("  = ($diff_ij)²")
    println("  = $energy_ij")
    
    # Now check from j to i perspective
    kappa_ji = get_transport_angle(field, face_j, face_i)
    p_ji = get_period_jump(field, face_j, face_i)
    
    # E_ji = (θ_j + κ_ji + (π/2)p_ji - θ_i)²
    diff_ji = field.theta[face_j] + kappa_ji + (π/2) * p_ji - field.theta[face_i]
    energy_ji = diff_ji^2
    
    println("\nFrom reverse direction:")
    println("κ_($face_j,$face_i) = $kappa_ji rad = $(rad2deg(kappa_ji))°")
    println("p_($face_j,$face_i) = $p_ji")
    println("\nEnergy term E_ji = (θ_j + κ_ji + (π/2)p_ji - θ_i)²")
    println("  = ($diff_ji)²")
    println("  = $energy_ji")
    
    println("\nEnergy values match: ", isapprox(energy_ij, energy_ji, atol=1e-10))
    
    # Test 4: Verify all face pairs have antisymmetric property
    println("\nTest 4: Check all face pairs have antisymmetric storage")
    println("-"^60)
    
    n_faces = topology.n_faces
    all_antisymmetric = true
    
    for i in 1:n_faces
        for (j, edge) in topology.dual_adj[i]
            p_ij = get_period_jump(field, i, j)
            p_ji = get_period_jump(field, j, i)
            
            if p_ij != -p_ji
                println("ERROR: p_($i,$j) = $p_ij but p_($j,$i) = $p_ji")
                all_antisymmetric = false
            end
        end
    end
    
    if all_antisymmetric
        println("✓ All face pairs satisfy antisymmetry: p_ji = -p_ij")
    else
        println("✗ Some face pairs do not satisfy antisymmetry!")
    end
    
    # Test 5: Verify only canonical direction is stored
    println("\nTest 5: Verify storage efficiency (only i < j stored)")
    println("-"^60)
    
    # Count stored keys
    stored_keys = collect(keys(field.period_jumps))
    n_stored = length(stored_keys)
    
    # Count dual edges
    n_dual_edges = 0
    for i in 1:n_faces
        for (j, edge) in topology.dual_adj[i]
            if i < j
                n_dual_edges += 1
            end
        end
    end
    
    println("Number of dual edges (i < j): $n_dual_edges")
    println("Number of stored period jumps: $n_stored")
    println("Storage matches edge count: ", n_stored == n_dual_edges)
    
    # Verify all stored keys have i < j
    all_canonical = all(key[1] < key[2] for key in stored_keys)
    println("All stored keys canonical (i < j): $all_canonical")
    
    # Verify we can still access both directions
    test_access = true
    for i in 1:n_faces
        for (j, edge) in topology.dual_adj[i]
            try
                p_ij = get_period_jump(field, i, j)
                p_ji = get_period_jump(field, j, i)
            catch
                println("ERROR: Cannot access p_($i,$j)")
                test_access = false
            end
        end
    end
    println("Can access all directions via helper: $test_access")
    
    # Test 6: Verify transport angles also use canonical storage
    println("\nTest 6: Verify transport angles also use canonical storage")
    println("-"^60)
    
    stored_kappa_keys = collect(keys(field.transport_angles))
    n_stored_kappa = length(stored_kappa_keys)
    
    println("Number of dual edges (i < j): $n_dual_edges")
    println("Number of stored transport angles: $n_stored_kappa")
    println("Storage matches edge count: ", n_stored_kappa == n_dual_edges)
    
    # Verify all stored keys have i < j
    all_kappa_canonical = all(key[1] < key[2] for key in stored_kappa_keys)
    println("All stored keys canonical (i < j): $all_kappa_canonical")
    
    # Verify antisymmetry works for transport angles
    test_kappa_antisym = true
    for i in 1:n_faces
        for (j, edge) in topology.dual_adj[i]
            kappa_ij = get_transport_angle(field, i, j)
            kappa_ji = get_transport_angle(field, j, i)
            
            if !isapprox(kappa_ij, -kappa_ji, atol=1e-10)
                println("ERROR: κ_($i,$j) = $kappa_ij but κ_($j,$i) = $kappa_ji")
                test_kappa_antisym = false
            end
        end
    end
    
    if test_kappa_antisym
        println("✓ All transport angles satisfy antisymmetry: κ_ji = -κ_ij")
    else
        println("✗ Some transport angles do not satisfy antisymmetry!")
    end
    
    println("\n" * "="^60)
    println("TEST COMPLETE")
    println("="^60)
end

main()
