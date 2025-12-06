using GeometryBasics
using DataStructures
using LinearAlgebra
using Statistics
using WGLMakie  # Interactive 3D plots in browser
using FileIO
using MeshIO
using SparseArrays
using WriteVTK
using NearestNeighbors

# Resolve Mesh name conflict by explicitly using GeometryBasics version
const Mesh = GeometryBasics.Mesh


# A helper struct to handle mesh topology explicitly
struct MeshTopology
    # Maps (v1, v2) sorted tuple to a unique edge index
    edge_map::Dict{Tuple{Int, Int}, Int} 
    # Maps face index to list of (neighbor_face_index, shared_edge_index)
    dual_adj::Vector{Vector{Tuple{Int, Int}}} 
    # We need this to iterate edges in the solver later
    edge_to_faces::Dict{Int, Vector{Int}}
    n_edges::Int
end



include("meshIO.jl")
include("plotters.jl")


##########
# ISSUES #
##########
# - singularities 
# - input mesh
    

"""
    generate_square_mesh(; nx=5, ny=5)

Generates a clean, connected triangle mesh of a unit square [0,1]x[0,1].
- nx, ny: Number of subdivisions (quads) along X and Y.
"""
function generate_square_mesh(; nx=6, ny=6)
    # 1. Generate Vertices (Row-major order)
    vs = Point3f[]
    for j in 0:ny      # 0 to 6 (7 rows)
        for i in 0:nx  # 0 to 6 (7 cols)
            push!(vs, Point3f(i/nx, j/ny, 0))
        end
    end
    
    # 2. Generate Faces
    fs = TriangleFace{Int}[]
    # Number of vertices per row
    row_stride = nx + 1
    
    for j in 1:ny      # Iterate cells in Y
        for i in 1:nx  # Iterate cells in X
            # Get indices of the 4 corners of the quad
            # Vertices array is 1-based
            v_bl = (j-1) * row_stride + i
            v_br = v_bl + 1
            v_tl = v_bl + row_stride
            v_tr = v_tl + 1
            
            # Add two triangles (forming a quad)
            # Ensure CCW winding order
            push!(fs, TriangleFace(v_bl, v_br, v_tr)) # Bottom-Right Tri
            push!(fs, TriangleFace(v_bl, v_tr, v_tl)) # Top-Left Tri
        end
    end
    
    return Mesh(vs, fs)
end

"""
    weld_mesh(mesh::Mesh; tol=1e-5)

Merges vertices that are closer than `tol` and updates face indices.
Returns a new welded Mesh.
"""
function weld_mesh(mesh::Mesh; tol=1e-5)
    vs = coordinates(mesh)
    fs = faces(mesh)
    
    # 1. Use a KDTree to find duplicate vertices efficiently
    tree = KDTree(vs)
    # Query all points against the tree to find their nearest neighbor index
    # (which will be the "canonical" index for that position)
    idxs, dists = nn(tree, vs)
    
    # 2. Remap faces to use the canonical indices
    new_faces = Vector{TriangleFace{Int}}()
    
    for f in fs
        # Map old vertex indices to new canonical indices
        v1 = idxs[f[1]]
        v2 = idxs[f[2]]
        v3 = idxs[f[3]]
        
        # Skip degenerate faces (where 2+ vertices merged into one)
        if v1 != v2 && v1 != v3 && v2 != v3
            push!(new_faces, TriangleFace(v1, v2, v3))
        end
    end
    
    # 3. (Optional) Clean up unused vertices 
    # For simplicity, we can keep the original vertex list since unused ones won't hurt topology
    # provided the faces only reference the canonical ones.
    
    return Mesh(vs, new_faces)
end



"""
    build_topology(mesh::Mesh)

Constructs the dual graph adjacency and unique edge numbering.
"""
function build_topology(mesh::Mesh)
    faces_list = decompose(TriangleFace{Int}, mesh)
    edge_map = Dict{Tuple{Int, Int}, Int}()
    edge_to_faces = Dict{Int, Vector{Int}}()
    
    edge_counter = 0
    
    # 1. Identify all unique edges and assign IDs
    for (f_idx, face) in enumerate(faces_list)
        vs = [face[1], face[2], face[3]]
        for k in 1:3
            v_a, v_b = minmax(vs[k], vs[mod1(k+1, 3)])
            if !haskey(edge_map, (v_a, v_b))
                edge_counter += 1
                edge_map[(v_a, v_b)] = edge_counter
                edge_to_faces[edge_counter] = Int[]
            end
            e_idx = edge_map[(v_a, v_b)]
            push!(edge_to_faces[e_idx], f_idx)
        end
    end
    
    # 2. Build Dual Adjacency (Face -> Face)
    n_faces = length(faces_list)
    dual_adj = [Vector{Tuple{Int, Int}}() for _ in 1:n_faces]
    
    for (e_idx, connected_faces) in edge_to_faces
        if length(connected_faces) == 2
            f1, f2 = connected_faces
            push!(dual_adj[f1], (f2, e_idx))
            push!(dual_adj[f2], (f1, e_idx))
        end
    end
    
    # Return the struct including edge_to_faces
    return MeshTopology(edge_map, dual_adj, edge_to_faces, edge_counter)
end


"""
    compute_spanning_forest(topology::MeshTopology, n_faces::Int, constrained_faces::Vector{Int})

Returns a BitVector of length `n_edges`. True indicates the edge is part of the 
spanning forest (where period jump p_ij should be fixed to 0).
"""
function compute_spanning_forest(topology::MeshTopology, n_faces::Int, constrained_faces::Vector{Int})
    # Track which faces have been visited (conquered by a tree)
    visited_faces = falses(n_faces)
    
    # The result: Edges where p_ij = 0
    # We use a BitVector for efficient indexing later in the solver
    fixed_edges = falses(topology.n_edges)
    
    # Queue for BFS (Dijkstra with uniform weights)
    # Stores the face index
    q = Queue{Int}()
    
    # Initialize roots (constrained faces) 
    for f_idx in constrained_faces
        enqueue!(q, f_idx)
        visited_faces[f_idx] = true
    end
    
    # If no constraints are provided (rare, but possible edge case), pick random root
    if isempty(constrained_faces) && n_faces > 0
        enqueue!(q, 1)
        visited_faces[1] = true
    end
    
    while !isempty(q)
        current_face = dequeue!(q)
        
        # Check neighbors in the dual graph
        for (neighbor_face, shared_edge) in topology.dual_adj[current_face]
            if !visited_faces[neighbor_face]
                # We found a new face via this edge. 
                # This edge becomes part of the spanning tree.
                fixed_edges[shared_edge] = true
                
                # Mark neighbor as visited and add to queue
                visited_faces[neighbor_face] = true
                enqueue!(q, neighbor_face)
            end
        end
    end
    
    # Validation: Ensure all faces reachable from constraints are visited.
    # Note: Disconnected components in the mesh without constraints in them 
    # won't be visited, but standard meshes are usually one component.
    
    return fixed_edges
end


##########################################################
"""
    get_local_frame(p1, p2, p3)

Returns (center, x_axis, y_axis) for a triangle defined by vertices p1, p2, p3.
X-axis is normalized p2-p1.
"""
function get_local_frame(p1::Point3f, p2::Point3f, p3::Point3f)
    center = (p1 + p2 + p3) / 3f0
    e1 = p2 - p1
    x_axis = normalize(e1)
    
    # Normal
    n = normalize(cross(e1, p3 - p1))
    # Y-axis (orthogonal to X and Normal)
    y_axis = cross(n, x_axis)
    
    return center, x_axis, y_axis
end

"""
    compute_kappas(mesh, topo)

Computes the transport angle kappa_ij for every edge.
Returns Dict: edge_idx -> angle (in radians).
"""
function compute_kappas(mesh::Mesh, topo::MeshTopology)
    kappas = Dict{Int, Float64}()
    vs = coordinates(mesh)
    fs = faces(mesh)
    
    for (e_idx, faces_list) in topo.edge_to_faces
        if length(faces_list) == 2
            i, j = faces_list # Face indices
            
            # Get vertices for Face i
            fi = fs[i]
            pi1, pi2, pi3 = vs[fi[1]], vs[fi[2]], vs[fi[3]]
            _, xi, yi = get_local_frame(pi1, pi2, pi3)
            
            # Get vertices for Face j
            fj = fs[j]
            pj1, pj2, pj3 = vs[fj[1]], vs[fj[2]], vs[fj[3]]
            _, xj, yj = get_local_frame(pj1, pj2, pj3)
            
            # To compute kappa_ij, we must unfold Face j onto Face i.
            # 1. Find shared edge vector
            # The edge e_idx connects two vertices. Let's find them.
            # (Optimization: We could look up edge_map, but geometry is safer)
            
            # We need the rotation that aligns the shared edge normals.
            # Actually, simpler method: 
            # The angle kappa is the rotation needed to align Frame i with Frame j 
            # relative to the shared edge.
            
            # Let's project xj onto the plane of i.
            # Since meshes are not flat, we unfold around the shared edge.
            
            # Robust approach:
            # 1. Measure angle of Shared Edge in Frame i (alpha_i)
            # 2. Measure angle of Shared Edge in Frame j (alpha_j)
            # 3. Ideally, shared edge should overlap. 
            # 4. kappa = alpha_j - alpha_i + pi (because shared edge is reversed in neighbor?)
            
            # Let's determine shared vertices u, v
            # Intersection of indices
            common = intersect(fi, fj)
            if length(common) != 2; continue; end
            p_u = vs[common[1]]
            p_v = vs[common[2]]
            edge_vec = normalize(p_v - p_u)
            
            # Angle of edge_vec in Frame i
            ang_i = atan(dot(edge_vec, yi), dot(edge_vec, xi))
            
            # Angle of edge_vec in Frame j
            ang_j = atan(dot(edge_vec, yj), dot(edge_vec, xj))
            
            # The difference represents the coordinate system mismatch
            # kappa_ij = angle_j - angle_i
            
            k_val = ang_j - ang_i
            
            # Wrap to (-pi, pi]
            k_val = mod2pi(k_val + pi) - pi
            
            kappas[e_idx] = k_val
        end
    end
    return kappas
end

"""
    constraints_to_angles(mesh, constraints_dict)

Converts 3D vector constraints to local 2D angles.
constraints_dict: FaceID -> Vec3f
"""
function constraints_to_angles(mesh::Mesh, constraints::Dict{Int, Vec3f})
    angle_constraints = Dict{Int, Float64}()
    vs = coordinates(mesh)
    fs = faces(mesh)
    
    for (f_idx, dir_vec) in constraints
        face = fs[f_idx]
        p1, p2, p3 = vs[face[1]], vs[face[2]], vs[face[3]]
        
        _, x_axis, y_axis = get_local_frame(p1, p2, p3)
        
        # Project 3D constraint onto 2D local frame
        local_x = dot(dir_vec, x_axis)
        local_y = dot(dir_vec, y_axis)
        
        angle_constraints[f_idx] = atan(local_y, local_x)
    end
    return angle_constraints
end


###########################################################################
# --- Matrix Assembly ---
function assemble_system(mesh, topo, fixed_edges, kappas, constrained_angles)
    n_faces = length(faces(mesh))
    n_edges = topo.n_edges

    # Map Variables
    theta_var_map = Dict{Int, Int}()
    p_var_map = Dict{Int, Int}()
    
    current_idx = 0
    
    # Theta Variables (Free Faces only)
    # here we find all free faces and assign them variable indices
    for f_idx in 1:n_faces
        if !haskey(constrained_angles, f_idx)
            current_idx += 1
            theta_var_map[f_idx] = current_idx
        end
    end
    
    # P Variables (Free Edges only)
    # here we find all free edges and assign them variable indices
    # - skip fixed edges, boundary edges, and edges between two constrained faces
    for e_idx in 1:n_edges
        # Skip boundary and fixed tree edges
        if length(topo.edge_to_faces[e_idx]) < 2; continue; end
        
        f1, f2 = topo.edge_to_faces[e_idx]
        is_cc = haskey(constrained_angles, f1) && haskey(constrained_angles, f2)
        
        if !fixed_edges[e_idx] && !is_cc
            current_idx += 1
            p_var_map[e_idx] = current_idx
        end
    end
    
    total_vars = current_idx
    
    # Triplets
    I, J, V = Int[], Int[], Float64[]
    b = zeros(total_vars)
    
    function add!(r, c, v)
        push!(I, r); push!(J, c); push!(V, v)
    end
    
    # Build Energy Gradient
    for (e_idx, faces_list) in topo.edge_to_faces
        if length(faces_list) != 2; continue; end
        i, j = faces_list
        k_ij = kappas[e_idx]

        # get the variable rows for theta_i, theta_j, p_ij
        row_i = get(theta_var_map, i, 0)
        row_j = get(theta_var_map, j, 0)
        row_p = get(p_var_map, e_idx, 0)
        
        # Get constant values for fixed variables
        # these terms will end up in the RHS
        th_i_val = get(constrained_angles, i, 0.0)
        th_j_val = get(constrained_angles, j, 0.0)
        
        # Determine fixed P value
        # these terms will end up in the RHS
        p_val = 0.0
        if row_p == 0
            # if a fixed edged, set p_val = 0
            if fixed_edges[e_idx]
                p_val = 0.0
                # p_val = round(-k_ij * 2.0/pi)
            # if between two constrained faces, compute rounding
            elseif haskey(constrained_angles, i) && haskey(constrained_angles, j)
                # Rounding between constraints
                term = (2.0/pi) * (th_j_val - th_i_val - k_ij)
                p_val = round(term)
            end
        end
        
        # RHS accumulation
        const_RHS = -k_ij
        if row_i == 0; 
            const_RHS -= th_i_val; 
        end
        if row_j == 0; 
            const_RHS += th_j_val; 
        end # Sign flip: -theta_j
        if row_p == 0; 
            const_RHS -= (pi/2)*p_val; 
        end
        
        # Add derivatives (Equations 2 & 3)
        # Row i
        if row_i > 0
            add!(row_i, row_i, 2.0); 
            
            if row_j>0
                add!(row_i, row_j, -2.0)
            end
            if row_p>0
                add!(row_i, row_p, pi)
            end

            b[row_i] += 2.0 * const_RHS
        end
        # Row j
        if row_j > 0
            if row_i>0; add!(row_j, row_i, -2.0); end
            add!(row_j, row_j, 2.0)
            if row_p>0; add!(row_j, row_p, -pi); end
            b[row_j] += -2.0 * const_RHS
        end
        # Row p
        if row_p > 0
            if row_i>0; add!(row_p, row_i, pi); end
            if row_j>0; add!(row_p, row_j, -pi); end
            add!(row_p, row_p, pi^2/2.0)
            b[row_p] += pi * const_RHS
        end
    end
    
    A = sparse(I, J, V, total_vars, total_vars)
    return A, b, theta_var_map, p_var_map
end

# --- Solver Logic ---
function local_gauss_seidel!(A, b, x, fixed_mask, start_node; max_iter=500000, tol=1e-5)
    q = Queue{Int}()
    in_queue = Set{Int}()
    
    rows = rowvals(A)
    vals = nonzeros(A)
    
    # Iterate over column 'start_node' (which is symmetric to row 'start_node')
    for idx in nzrange(A, start_node)
        neighbor = rows[idx]
        if !fixed_mask[neighbor] && !(neighbor in in_queue)
            enqueue!(q, neighbor)
            push!(in_queue, neighbor)
        end
    end

    iter = 0
    diag_A = diag(A)

    while !isempty(q) && iter < max_iter
        iter += 1
        k = dequeue!(q)
        delete!(in_queue, k)
        
        # We only process 'k' if it is free (redundant check if queue logic is correct, but safe)
        if !fixed_mask[k]
            # Calculate Residual: r_k = b[k] - Sum(A_kj * x_j)
            # Optimization: Only sum over non-zeros
            Ax_k = 0.0
            for idx in nzrange(A, k)
                Ax_k += vals[idx] * x[rows[idx]]
            end
            
            r_k = b[k] - Ax_k
            
            if abs(r_k) > tol
                # Update variable
                delta = r_k / diag_A[k]
                x[k] += delta
                
                # If we changed k significantly, we must notify ITS neighbors
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
    
    # Optional: Print occasionally to verify it's working (e.g., if iter > 100)
    # if iter > 100; println("  Relaxed with $iter iterations"); end
    @info "  Local Gauss-Seidel relaxed $iter iterations"
end


function solve_greedy(A, b, p_var_map)
    println("  Factorizing system...")    
    F = cholesky(Hermitian(A)) # Assumes A is SPD
    x = F \ b
    
    println("  Rounding integers...")
    int_indices = collect(values(p_var_map))
    fixed_mask = falses(length(b))
    
    # Rounding Loop
    while true
        best_idx = -1
        min_err = Inf
        best_val = 0.0
        
        # Find best candidate (closest to integer)
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
        
        # Fix it
        x[best_idx] = best_val
        fixed_mask[best_idx] = true
        
        # Relax
        local_gauss_seidel!(A, b, x, fixed_mask, best_idx)
    end
    
    return x
end





function solve_greedy_global(A, b, p_var_map)
    println("  Factorizing system...")
    # F = cholesky(Hermitian(A))
    F = deepcopy(A)  
    x = F \ b
    
    println("  Rounding integers...")
    int_indices = collect(values(p_var_map))
    fixed_mask = falses(length(b))
    
    # Pre-allocate queue with ALL variables for robust updates
    all_vars = 1:length(b)
    
    while true
        # 1. Find best rounding candidate
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
        
        # 2. Fix variable
        x[best_idx] = best_val
        fixed_mask[best_idx] = true
        
        # 3. ROBUST RELAXATION: Add EVERYONE to the queue
        # This prevents the solver from "sleeping" after 2 iterations
        q = Queue{Int}()
        in_queue = falses(length(b))
        for i in all_vars
            if !fixed_mask[i]
                enqueue!(q, i)
                in_queue[i] = true
            end
        end
        
        # Run relaxation
        iter = run_gauss_seidel_queue!(A, b, x, fixed_mask, q, in_queue)
        
        # Debug print to prove it's working
        println("Rounded var $best_idx. Relaxed $iter iterations.")
    end
    return x
end

# Ensure you have this helper function defined too!
function run_gauss_seidel_queue!(A, b, x, fixed_mask, q, in_queue; max_iter=1_000_000, tol=1e-5)
    iter = 0
    diag_A = diag(A)
    rows = rowvals(A)
    vals = nonzeros(A)
    
    while !isempty(q) && iter < max_iter
        iter += 1
        k = dequeue!(q)
        in_queue[k] = false 
        
        if !fixed_mask[k]
            # Compute Ax_k manually for speed
            Ax_k = 0.0
            for idx in nzrange(A, k)
                Ax_k += vals[idx] * x[rows[idx]]
            end
            
            r_k = b[k] - Ax_k
            if abs(r_k) > tol
                x[k] += r_k / diag_A[k]
                
                # Add neighbors
                for idx in nzrange(A, k)
                    neighbor = rows[idx]
                    if !fixed_mask[neighbor] && !in_queue[neighbor]
                        enqueue!(q, neighbor)
                        in_queue[neighbor] = true
                    end
                end
            end
        end
    end
    return iter
end





function log_euler_char(mesh_info, topology)
    # compute the euler characteristic
    V = length(coordinates(mesh_info))
    E = topology.n_edges
    F = length(faces(mesh_info))
    chi = V - E + F
    println("|V| = $V, |E| = $E, |F| = $F, χ = $chi")
end






#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
# user options:
# Load the mesh (Gmsh .msh file, or .obj, .stl)

# filename = "simple-square.msh"  
# filename = "unit-square.msh"  
# filename = "hole.msh"
# filename = "teapot.msh"
# filename = "300_polygon_sphere_100mm.msh"

filepath = joinpath("SimpleAutoQuad", "output", "meshes", filename)
mesh_info = load_triangulation(filepath)



# println("Generating clean square mesh...")
# mesh_info = generate_square_mesh(nx=6, ny=6) # Adjust resolution as needed

# @info "Welding mesh to remove duplicate vertices..."
# mesh_info = weld_mesh(mesh_info_raw; tol=1)
# @info "Raw Stats:"
# log_euler_char(mesh_info_raw)
# @info "Welded Stats:"
# readline()


# Define constraints (hardcode directions)
# randomly pick n constrained faces
n = 1
total_faces = length(faces(mesh_info))
constrained_faces = rand(1:total_faces, n)
# constraint_vecs = Dict(1 => Vec3f(1, 1, 0), 4 => Vec3f(0, 1, 0)) 
# pick random directions for each constrained face
constraint_vecs = Dict{Int, Vec3f}()
rand_dir = normalize(Vec3f(1, 0, 0))
for f_idx in constrained_faces
    constraint_vecs[f_idx] = rand_dir
end



# 3. Build the topology 
topo = build_topology(mesh_info)
kappas = compute_kappas(mesh_info, topo)

# Log Euler Characteristic
log_euler_char(mesh_info, topo)


# 4. Spanning Forest
fixed_edges = compute_spanning_forest(topo, length(faces(mesh_info)), constrained_faces)

# 4 part 2
# let's pick only one edge per triangle to be fixed
for f_idx in 1:length(faces(mesh_info))
    face = faces(mesh_info)[f_idx]
    vs_face = [face[1], face[2], face[3]]
    # find edges of this face
    face_edges = Tuple{Int, Int}[]
    for k in 1:3
        v_a, v_b = minmax(vs_face[k], vs_face[mod1(k+1, 3)])
        push!(face_edges, (v_a, v_b))
    end
    # find the first edge that is not fixed yet and fix it
    for edge in face_edges
        e_idx = topo.edge_map[edge]
        if !fixed_edges[e_idx]
            fixed_edges[e_idx] = true
            break
        end
    end
end




# 5. Convert Constraints to Angles
constrained_angles = constraints_to_angles(mesh_info, constraint_vecs)

# 6. Assemble System
println("Assembling system...")
A, b, theta_map, p_map = assemble_system(mesh_info, topo, fixed_edges, kappas, constrained_angles)

# 7. Solve
println("Solving...")
x_sol = solve_greedy_global(A, b, p_map)

# 8. Extract Results
final_thetas = zeros(length(faces(mesh_info)))
final_p_values = Dict{Int, Int}()

for (e_idx, mat_idx) in p_map
    # final_p_values[e_idx] = round(Int, x_sol[mat_idx])
    final_p_values[e_idx] = x_sol[mat_idx]
end

# Fill free variables
for (f_idx, mat_idx) in theta_map
    final_thetas[f_idx] = x_sol[mat_idx]
end
# Fill fixed constraints
for (f_idx, val) in constrained_angles
    final_thetas[f_idx] = val
end

println("Final Integer P Values (Edge Index => P_ij):")
for (e_idx, p_val) in final_p_values
    println("  Edge $e_idx => $p_val")
end



# 10. Compute Singularities
include("SingularityFinder.jl")
sings = compute_singularities(mesh_info, topo, final_p_values)
println("Found $(length(sings)) singularities.")
println("Singularities (Vertex Index => Index Sum):")
for (v_idx, sing_type) in sings
    println("  Vertex $v_idx => $sing_type")
end




println("Done! Visualizing...")







# 9. Save to VTU
output_path = joinpath("SimpleAutoQuad", "output", "framefields", "cross_field_results")
save_cross_field_to_vtu(output_path, mesh_info, final_thetas, constraint_vecs, topo, p_map, x_sol)

# 10. Visualize
fig = plot_cross_field(mesh_info, final_thetas, constrained_faces; 
                       constraint_vecs=constraint_vecs, 
                       sings=sings,
                       topo=topo,
                       p_values=final_p_values)
fig

# 11. Visualize Smoothness
fig_smooth = plot_field_smoothness(mesh_info, final_thetas, topo, kappas, final_p_values)
fig_smooth