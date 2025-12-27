module Cutting

using ..Types
using DataStructures

export compute_cut_graph, propagate_orientations

# --- Helper: Build Vertex Adjacency (Primal Graph) ---
function build_primal_adjacency(topo::MeshTopology)
    adj = Dict{Int, Set{Int}}()
    for i in 1:length(topo.vertices)
        adj[i] = Set{Int}()
    end
    for (v1, v2, v3) in topo.faces
        push!(adj[v1], v2); push!(adj[v1], v3)
        push!(adj[v2], v1); push!(adj[v2], v3)
        push!(adj[v3], v1); push!(adj[v3], v2)
    end
    return adj
end

# --- Helper: Dijkstra to Set with Barriers ---
# Added `barrier_edges` argument
function shortest_path_to_set(
    start_node::Int, 
    target_set::Set{Int}, 
    adj::Dict{Int, Set{Int}}, 
    vertices::Vector{Point3D},
    barrier_edges::Set{Tuple{Int,Int}} # <--- NEW ARGUMENT
)
    # If start is already in target, return empty path
    if start_node in target_set
        return Int[]
    end

    q = PriorityQueue{Int, Float64}()
    dist = Dict{Int, Float64}()
    parent = Dict{Int, Int}()
    
    # Initialize
    q[start_node] = 0.0
    dist[start_node] = 0.0
    
    found_target = -1
    
    while !isempty(q)
        u = dequeue!(q)
        
        if u in target_set
            found_target = u
            break
        end
        
        p_u = vertices[u]
        
        for v in adj[u]
            # [FIX] Check for Barrier: Do not cross existing cuts!
            edge_key = minmax(u, v)
            if edge_key in barrier_edges
                continue 
            end

            p_v = vertices[v]
            w = sqrt((p_u.x - p_v.x)^2 + (p_u.y - p_v.y)^2 + (p_u.z - p_v.z)^2)
            
            new_dist = dist[u] + w
            
            if !haskey(dist, v) || new_dist < dist[v]
                dist[v] = new_dist
                parent[v] = u
                q[v] = new_dist
            end
        end
    end
    
    # Reconstruct Path
    path = Int[]
    if found_target != -1
        curr = found_target
        while curr != start_node
            push!(path, curr)
            curr = parent[curr]
        end
        push!(path, start_node)
        reverse!(path)
    end
    return path
end

function compute_cut_graph(topo::MeshTopology, singularities::Vector{Tuple{Int, Float64}})
    println("\n--- Computing Cut Graph ---")
    
    # === STEP 1: Dual Spanning Tree & Initial Cuts ===
    # This generates the rough topological cuts
    n_faces = length(topo.faces)
    visited_face = falses(n_faces)
    tree_edges = Set{Tuple{Int,Int}}()
    
    queue = Int[1]
    visited_face[1] = true
    
    while !isempty(queue)
        f_curr = popfirst!(queue)
        for (f_neigh, edge_key) in topo.dual_adj[f_curr]
            if !visited_face[f_neigh]
                visited_face[f_neigh] = true
                push!(tree_edges, edge_key)
                push!(queue, f_neigh)
            end
        end
    end
    
    raw_cut_edges = Set{Tuple{Int,Int}}()
    for i in 1:n_faces
        for (j, key) in topo.dual_adj[i]
            if i < j && !(key in tree_edges)
                push!(raw_cut_edges, key)
            end
        end
    end
    
    # === Identify Boundary & Detect Open/Closed ===
    boundary_verts = Set{Int}()
    all_half_edges = Dict{Tuple{Int,Int}, Int}()
    for tri in topo.faces
        e1, e2, e3 = minmax(tri[1],tri[2]), minmax(tri[2],tri[3]), minmax(tri[3],tri[1])
        all_half_edges[e1] = get(all_half_edges, e1, 0) + 1
        all_half_edges[e2] = get(all_half_edges, e2, 0) + 1
        all_half_edges[e3] = get(all_half_edges, e3, 0) + 1
    end
    for (e, count) in all_half_edges
        if count == 1
            push!(boundary_verts, e[1])
            push!(boundary_verts, e[2])
        end
    end

    # Explicit Detection
    is_closed_mesh = isempty(boundary_verts)
    println(is_closed_mesh ? "  [Detected Closed Mesh]" : "  [Detected Open Mesh]")

    # === STEP 2: Prune Open Paths ===
    final_cuts = copy(raw_cut_edges)
    
    # We only protect boundaries. 
    # We explicitly UNPROTECT singularities so jagged paths are removed.
    protected_verts = copy(boundary_verts)
    
    changed = true
    while changed
        changed = false
        cut_adj = Dict{Int, Vector{Int}}()
        for (u, v) in final_cuts
            if !haskey(cut_adj, u); cut_adj[u] = Int[]; end
            if !haskey(cut_adj, v); cut_adj[v] = Int[]; end
            push!(cut_adj[u], v)
            push!(cut_adj[v], u)
        end
        
        to_remove = Set{Tuple{Int,Int}}()
        for (v, neighbors) in cut_adj
            # Prune if valence 1 AND not protected
            if length(neighbors) == 1 && !(v in protected_verts)
                neighbor = neighbors[1]
                edge = minmax(v, neighbor)
                push!(to_remove, edge)
                changed = true
            end
        end
        if !isempty(to_remove)
            setdiff!(final_cuts, to_remove)
        end
    end
    
    println("  After Pruning: $(length(final_cuts)) edges")

    # === STEP 3: Connect Singularities (Geodesic Paths) ===
    primal_adj = build_primal_adjacency(topo)
    
    # Define Targets: Boundaries + Any loops that survived pruning (e.g., on Torus)
    cut_verts = copy(protected_verts) 
    for (u, v) in final_cuts
        push!(cut_verts, u)
        push!(cut_verts, v)
    end
    
    # === HANDLE CLOSED MESH (Cube/Sphere) ===
    # If closed, pruning (Step 2) removed everything because there was no boundary.
    # We must manually pick the first singularity as the "Root Anchor".
    if is_closed_mesh && isempty(cut_verts) && !isempty(singularities)
        root_s = singularities[1][1]
        push!(cut_verts, root_s)
        println("  [Closed Mesh] Seeding cut graph with Singularity $root_s")
    end
    
    added_paths = 0
    for (s_idx, s_val) in singularities
        if s_idx in cut_verts
            continue
        end
        
        # Draw straight path to the nearest anchor (Boundary or Root)
        # PASS final_cuts as BARRIER to prevent crossing existing cuts
        path = shortest_path_to_set(s_idx, cut_verts, primal_adj, topo.vertices, final_cuts)
        
        if !isempty(path)
            added_paths += 1
            for i in 1:(length(path)-1)
                u, v = path[i], path[i+1]
                edge = minmax(u, v)
                push!(final_cuts, edge)
                push!(cut_verts, u)
                push!(cut_verts, v)
            end
        else
            println("  Warning: Could not connect Singularity $s_idx to graph!")
        end
    end
    
    println("  Added $added_paths paths to connect singularities.")
    println("Final Cut Graph size: $(length(final_cuts)) edges")
    
    return final_cuts
end



function propagate_orientations(field::CrossField, cut_edges::Set{Tuple{Int,Int}})
    topo = field.topology
    n_faces = length(topo.faces)
    
    r = zeros(Int, n_faces)
    visited = falses(n_faces)
    
    q = Queue{Int}()
    
    root = 1
    enqueue!(q, root)
    visited[root] = true
    r[root] = 0
    
    println("\n--- Propagating Global Orientation ---")
    
    while !isempty(q)
        u = dequeue!(q)
        
        for (v, edge_key) in topo.dual_adj[u]
            if edge_key in cut_edges
                continue
            end
            
            if !visited[v]
                p = get(field.period_jumps, edge_key, 0)
                
                if u < v
                    r[v] = r[u] - p
                else
                    r[v] = r[u] + p
                end
                
                visited[v] = true
                enqueue!(q, v)
            end
        end
    end
    
    cut_rotations = Dict{Tuple{Int,Int}, Int}()
    
    for edge in cut_edges
        u, v = (edge[1] < edge[2]) ? (edge[1], edge[2]) : (edge[2], edge[1])
        
        if visited[u] && visited[v]
            p = get(field.period_jumps, edge, 0)
            k = r[v] - r[u] + p
            cut_rotations[(u, v)] = k
        else
            println("Warning: Cut edge $edge connects unvisited faces (disconnected component?)")
        end
    end
    
    println("Propagated orientation to $(count(visited)) faces.")
    println("Computed rotations for $(length(cut_rotations)) cut edges.")
    
    return r, cut_rotations
end

end
