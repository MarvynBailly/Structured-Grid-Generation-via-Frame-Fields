module Parameterization

using ..Types
using ..Topology
using DataStructures # For Queue
using LinearAlgebra  # Added for dot, cross, normalize

export propagate_orientations

# --- Helper: Compute Global Vector (Duplicated from Plotting to avoid circular dependency) ---
function get_global_vector_local(topo::MeshTopology, f_idx::Int, theta::Float64, r::Int)
    # 1. Get Reference Frame
    tri = topo.faces[f_idx]
    ref_edge = topo.face_ref_edges[f_idx] # (v_start, v_end)
    
    p1 = topo.vertices[ref_edge[1]]
    p2 = topo.vertices[ref_edge[2]]
    
    # Basis X (e1) = Normalized Reference Edge
    e1 = [p2.x - p1.x, p2.y - p1.y, p2.z - p1.z]
    normalize!(e1)
    
    # Normal (n)
    v1 = topo.vertices[tri[1]]
    v2 = topo.vertices[tri[2]]
    v3 = topo.vertices[tri[3]]
    u = [v2.x - v1.x, v2.y - v1.y, v2.z - v1.z]
    v = [v3.x - v1.x, v3.y - v1.y, v3.z - v1.z]
    n = cross(u, v)
    normalize!(n)
    
    # Basis Y (e2) = n x e1
    e2 = cross(n, e1)
    
    # 2. Total Angle = Optimized Theta + Integer Rotation * 90 degrees
    total_angle = theta + (r * π / 2.0)
    
    # 3. Rotate in Tangent Plane
    v_global = cos(total_angle) .* e1 .+ sin(total_angle) .* e2
    return v_global
end

"""
    propagate_orientations(field::CrossField, cut_graph::Set{Tuple{Int,Int}})

Assigns an integer rotation r_i ∈ {0, 1, 2, 3} to each face i.
Includes a geometric check to verify that vectors align across edges.
"""
function propagate_orientations(field::CrossField, cut_graph::Set{Tuple{Int,Int}}; 
                                verbose::Bool=true, 
                                callback::Union{Function, Nothing}=nothing)
    println("\n--- Propagating Orientations (Section 5) ---")
    
    topo = field.topology
    n_faces = length(topo.faces)
    
    # rotations[i] stores the integer rotation (0,1,2,3) for face i
    rotations = fill(-1, n_faces)
    
    q = Queue{Int}()
    iter = 0
    
    # Statistics for the check
    alignment_errors = 0
    min_dot = 1.0
    
    # Iterate to ensure coverage (handle disconnected components if any)
    for root_face in 1:length(topo.faces)
        if rotations[root_face] != -1
            continue
        end
        
        # Start BFS from this root
        enqueue!(q, root_face)
        rotations[root_face] = 0
        
        while !isempty(q)
            u = dequeue!(q)
            iter += 1
            
            if !isnothing(callback)
                callback(rotations, iter)
            end
            
            r_u = rotations[u]
            
            # Pre-compute vector for u (for checking)
            vec_u = get_global_vector_local(topo, u, field.theta[u], r_u)
            
            # Iterate over neighbors
            for (v, edge_key) in topo.dual_adj[u]
                
                # 1. STOP at Cut Graph
                if edge_key in cut_graph
                    continue
                end
                
                # 2. Check if already visited
                if rotations[v] != -1
                    continue
                end
                
                # 3. Retrieve Period Jump (Computed from Theta)
                # The solver may have fixed p=0 on the spanning tree, but theta might have wound around.
                # We recover the effective jump from the geometry of the solution.
                edge_key = minmax(u, v)
                kappa = get(field.transport_angles, edge_key, 0.0)
                
                # Effective kappa for u -> v
                kappa_uv = (u < v) ? kappa : -kappa
                
                # Calculate effective period jump
                # theta_v - theta_u ≈ kappa_uv + p_eff * (pi/2)
                diff = field.theta[v] - field.theta[u]
                p_eff = round(Int, (diff - kappa_uv) / (π/2))
                
                # 4. Calculate neighbor rotation
                # r_v = r_u - p_eff
                rotations[v] = mod(r_u - p_eff, 4)
                
                # --- GEOMETRIC CHECK ---
                # Calculate the actual 3D vector for neighbor v
                vec_v = get_global_vector_local(topo, v, field.theta[v], rotations[v])
                
                # Dot product should be close to 1.0 (Aligned)
                alignment = dot(vec_u, vec_v)
                min_dot = min(min_dot, alignment)
                
                if alignment < 0.7 # Threshold for "roughly aligned"
                    alignment_errors += 1
                    if verbose && alignment_errors <= 10 # Don't spam
                        println("  [MISALIGNMENT] Edge $u -> $v")
                        println("    p_eff: $p_eff")
                        println("    r_u: $r_u, r_v: $(rotations[v])")
                        println("    Dot Product: $(round(alignment, digits=4))")
                    end
                end
                # -----------------------
                # -----------------------
                
                enqueue!(q, v)
            end
        end
    end
    
    if verbose
        println("Propagation Complete.")
        println("  Minimum Alignment (Dot): $(round(min_dot, digits=4))")
        if alignment_errors > 0
            println("  WARNING: Found $alignment_errors edges with poor vector alignment (< 0.7).")
        else
            println("  Success: All vectors align smoothly.")
        end
    end
    
    return rotations
end

end