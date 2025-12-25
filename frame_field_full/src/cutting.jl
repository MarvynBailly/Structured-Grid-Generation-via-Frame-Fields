module Cutting

using ..Types
using DataStructures

export compute_cut_graph, propagate_orientations

function compute_cut_graph(topo::MeshTopology)
    println("\n--- Computing Cut Graph (Dual Spanning Tree) ---")
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
            if i < j
                if !(key in tree_edges)
                    push!(raw_cut_edges, key)
                end
            end
        end
    end
    
    println("Initial Cut Graph size: $(length(raw_cut_edges)) edges")
    
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
    
    final_cuts = copy(raw_cut_edges)
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
            if length(neighbors) == 1 && !(v in boundary_verts)
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
