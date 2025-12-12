"""
    test_kappa_computation.jl

Test the transport angle computation with a simple example
to verify it matches the formula from the slides.
"""

using Pkg
Pkg.activate(".")

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using GeometryBasics
using LinearAlgebra

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
    
    vs = coordinates(topology.mesh)
    
    # Pick a specific adjacent pair to verify manually
    face_i = 1
    face_j = 4
    
    println("\n" * "="^60)
    println("MANUAL VERIFICATION OF κ COMPUTATION")
    println("="^60)
    println("\nTesting transport angle between face $face_i and face $face_j")
    
    # Get shared edge
    shared_edge = nothing
    for (neighbor, edge) in topology.dual_adj[face_i]
        if neighbor == face_j
            shared_edge = edge
            break
        end
    end
    
    if shared_edge === nothing
        println("Faces $face_i and $face_j are not adjacent!")
        return
    end
    
    println("\nShared edge: $shared_edge")
    
    # Get reference edges
    ref_i = topology.face_reference_edges[face_i]
    ref_j = topology.face_reference_edges[face_j]
    
    println("Reference edge of face $face_i: $ref_i")
    println("Reference edge of face $face_j: $ref_j")
    
    # Compute angles manually
    ref_i_p1 = vs[ref_i[1]]
    ref_i_p2 = vs[ref_i[2]]
    ref_i_vec = (ref_i_p2[1] - ref_i_p1[1], ref_i_p2[2] - ref_i_p1[2])
    ref_i_angle = atan(ref_i_vec[2], ref_i_vec[1])
    
    shared_p1 = vs[shared_edge[1]]
    shared_p2 = vs[shared_edge[2]]
    shared_vec = (shared_p2[1] - shared_p1[1], shared_p2[2] - shared_p1[2])
    shared_angle = atan(shared_vec[2], shared_vec[1])
    
    phi_i = shared_angle - ref_i_angle
    
    ref_j_p1 = vs[ref_j[1]]
    ref_j_p2 = vs[ref_j[2]]
    ref_j_vec = (ref_j_p2[1] - ref_j_p1[1], ref_j_p2[2] - ref_j_p1[2])
    ref_j_angle = atan(ref_j_vec[2], ref_j_vec[1])
    
    phi_j = shared_angle - ref_j_angle
    
    kappa_ij_manual = phi_j - phi_i + π
    
    # Adjust to (-π, π]
    while kappa_ij_manual > π
        kappa_ij_manual -= 2π
    end
    while kappa_ij_manual <= -π
        kappa_ij_manual += 2π
    end
    
    # Get computed value using helper function
    kappa_ij_computed = get_transport_angle(field, face_i, face_j)
    
    println("\n" * "-"^60)
    println("ANGLES:")
    println("-"^60)
    println("Reference edge angle (face $face_i): $(rad2deg(ref_i_angle))°")
    println("Reference edge angle (face $face_j): $(rad2deg(ref_j_angle))°")
    println("Shared edge angle: $(rad2deg(shared_angle))°")
    println("\nφ_i (ref to shared in face $face_i): $(rad2deg(phi_i))°")
    println("φ_j (ref to shared in face $face_j): $(rad2deg(phi_j))°")
    
    println("\n" * "-"^60)
    println("TRANSPORT ANGLE κ_ij:")
    println("-"^60)
    println("Formula: κ_ij = φ_j - φ_i + π")
    println("Manual computation: $(rad2deg(kappa_ij_manual))°")
    println("Computed value:     $(rad2deg(kappa_ij_computed))°")
    println("Match: ", isapprox(kappa_ij_manual, kappa_ij_computed, atol=1e-10))
    
    # Also verify the reverse direction
    kappa_ji_computed = get_transport_angle(field, face_j, face_i)
    println("\nReverse direction κ_ji: $(rad2deg(kappa_ji_computed))°")
    println("Should be -κ_ij: ", isapprox(kappa_ji_computed, -kappa_ij_computed, atol=1e-10))
    
    println("\n" * "="^60)
    println("VERIFICATION COMPLETE")
    println("="^60)
end

main()
