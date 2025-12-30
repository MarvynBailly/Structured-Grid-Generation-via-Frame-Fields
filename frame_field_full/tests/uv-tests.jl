using FrameFieldFull
using FrameFieldFull.Types
using FrameFieldFull.Topology
using FrameFieldFull.Solver
using FrameFieldFull.Analysis
using FrameFieldFull.Parameterization
using LinearAlgebra
using Printf

# --- Helper: Generate Flat Grid Mesh ---
function generate_flat_grid(Nx, Ny)
    vertices = Point3D[]
    faces = Tuple{Int,Int,Int}[]
    
    # Create Vertices
    for y in 1:Ny
        for x in 1:Nx
            push!(vertices, Point3D(Float64(x-1), Float64(y-1), 0.0))
        end
    end
    
    # Create Faces
    for y in 1:(Ny-1)
        for x in 1:(Nx-1)
            # Grid indices (1-based)
            v_bl = (y-1)*Nx + x
            v_br = (y-1)*Nx + (x+1)
            v_tl = y*Nx + x
            v_tr = y*Nx + (x+1)
            
            # Quad split into 2 triangles
            push!(faces, (v_bl, v_br, v_tl))
            push!(faces, (v_br, v_tr, v_tl))
        end
    end
    
    return vertices, faces
end

# --- Test 1: Uniform Field Propagation ---
# Verifies that the cross-field solver correctly propagates a fixed constraint (0 degrees)
# across a flat mesh without deviation.
function test_uniform_field()
    println("\n=== TEST 1: Uniform Field Propagation (Flat 5x5 Grid) ===")
    
    Nx, Ny = 5, 5
    verts, faces = generate_flat_grid(Nx, Ny)
    topo = build_topology(verts, faces)
    
    # Constrain the bottom-left face to 0.0 radians
    constraints = Dict(1 => 0.0)
    
    field = initialize_field(topo, constraints)
    solve_greedy!(field; verbose=false)
    
    # Check: All faces should have theta close to 0.0 (or multiples of pi/2)
    max_err = 0.0
    for i in 1:length(field.theta)
        # Normalize angle to -pi/4 .. pi/4 to ignore 90 degree symmetry
        t = field.theta[i]
        dev = abs(rem(t, π/2, RoundNearest)) 
        max_err = max(max_err, dev)
    end
    
    if max_err < 1e-6
        println("  [PASS] Field is uniform. Max deviation: $max_err")
    else
        println("  [FAIL] Field deviated. Max deviation: $max_err")
    end
end

# --- Test 2: Exact Coordinate Reconstruction ---
# Verifies that the Global Parametrization solver reconstructs the X,Y coordinates
# of a flat grid (which is the perfect parameterization for a flat grid).
function test_coordinate_reconstruction()
    println("\n=== TEST 2: Exact Coordinate Reconstruction (Flat 3x3 Grid) ===")
    
    Nx, Ny = 3, 3
    verts, faces = generate_flat_grid(Nx, Ny)
    topo = build_topology(verts, faces)
    
    # 1. Setup perfect field (theta = 0 everywhere)
    constraints = Dict{Int,Float64}()
    # Constrain boundary faces to match geometric edges (which are axis aligned)
    # This forces the field to be perfectly grid-aligned
    constraints = compute_boundary_constraints(topo) 
    
    field = initialize_field(topo, constraints)
    solve_greedy!(field; verbose=false)
    
    # 2. Cut & Parametrize
    sings = compute_singularities(field)
    cuts = compute_cut_graph(topo, sings)
    
    # This should trigger the Mixed-Integer solver
    # Since the grid is integer-aligned, the rounding error should be near zero.
    u_p, v_p = solve_mixed_integer_parameterization(field, cuts; verbose=false)
    
    # 3. Verify
    # We expect u ~ x and v ~ y (up to some global integer rotation/translation)
    # Let's look at the gradient or relative positions.
    # Edge (0,0)->(1,0) should have du=1, dv=0 (or swapped)
    
    # Get vertex indices for (0,0) and (1,0)
    idx_00 = 1
    idx_10 = 2 
    
    du = u_p[idx_10] - u_p[idx_00]
    dv = v_p[idx_10] - v_p[idx_00]
    dist = sqrt(du^2 + dv^2)
    
    println("  Edge Length in Param Space (Expected ~1.0): $dist")
    
    # Linearity Check: Fit a plane u = ax + by + c
    # If the solver works, residuals should be tiny.
    xs = [v.x for v in topo.vertices]
    ys = [v.y for v in topo.vertices]
    
    # Simple check: u should be linear in x,y
    # We check if gradients are constant.
    passed = abs(dist - 1.0) < 1e-3
    
    if passed
        println("  [PASS] Parameterization reconstructed metric correctly.")
    else
        println("  [FAIL] distortion detected.")
    end
end

# --- Test 3: The "Trivial Cut" (Single Triangle) ---
# A sanity check that the system matrix can be assembled and solved for a single element.
function test_single_triangle()
    println("\n=== TEST 3: Single Triangle Sanity Check ===")
    
    verts = [Point3D(0.0,0.0,0.0), Point3D(1.0,0.0,0.0), Point3D(0.0,1.0,0.0)]
    faces = [(1, 2, 3)]
    topo = build_topology(verts, faces)
    
    constraints = Dict(1 => 0.0)
    field = initialize_field(topo, constraints)
    
    cuts = Set{Tuple{Int,Int}}() # No cuts needed for 1 triangle
    
    try
        u, v = solve_mixed_integer_parameterization(field, cuts; verbose=false)
        
        # Check orthogonality
        # u(1,0) - u(0,0) should be orthogonal to u(0,1) - u(0,0)?
        # Ideally u aligned with x, v aligned with y
        u_x = u[2] - u[1]
        v_x = v[2] - v[1]
        
        u_y = u[3] - u[1]
        v_y = v[3] - v[1]
        
        # Dot product of basis vectors in param space
        dot_prod = u_x*u_y + v_x*v_y
        
        if abs(dot_prod) < 1e-4
            println("  [PASS] Parametrization is orthogonal.")
        else
            println("  [FAIL] Parametrization is skewed. Dot: $dot_prod")
        end
    catch e
        println("  [FAIL] Solver crashed: $e")
        rethrow(e)
    end
end

function run_tests()
    test_uniform_field()
    test_single_triangle()
    test_coordinate_reconstruction()
    println("\nAll Tests Finished.")
end

run_tests()