"""
    cut_graph.jl

Compute cut graphs for mesh topology handling in global parameterization.

For closed surfaces (sphere, torus), we must cut the mesh to create a topological
disk that can be mapped to 2D. The cut graph connects singularities and turns
the closed surface into a flat sheet.

Key concepts:
- Cut graph: A tree-like structure connecting all singularities
- Seam edges: The edges that are "cut" to unfold the mesh
- Cross-boundary compatibility: Ensures grid lines match across cuts
"""

using GeometryBasics
using DataStructures

"""
    CutGraph

Represents a cut graph for mesh parameterization.

# Fields
- `seam_edges::Set{Tuple{Int,Int}}` - Edges that form the cut (canonical form: v1 < v2)
- `singularity_tree::Dict{Int,Vector{Int}}` - Tree connecting singularities (vertex → neighbors)
- `boundary_loops::Vector{Vector{Int}}` - Boundary loops after cutting
"""
struct CutGraph
    seam_edges::Set{Tuple{Int,Int}}
    singularity_tree::Dict{Int,Vector{Int}}
    boundary_loops::Vector{Vector{Int}}
end


"""
    compute_cut_graph(
        mesh::Mesh,
        singularities::Dict{Int,Float64};
        boundary_vertices::Set{Int}=Set{Int}()
    ) -> CutGraph

Compute a cut graph that connects all singularities and makes the mesh
a topological disk.

# Algorithm
1. If mesh has boundary, connect singularities to boundary
2. If closed surface, build spanning tree connecting all singularities
3. Add additional cuts to ensure topological disk (Euler characteristic)

# Arguments
- `mesh`: Triangular mesh
- `singularities`: Dict mapping vertex index → singularity index
- `boundary_vertices`: Set of vertices on mesh boundary (empty for closed surfaces)

# Returns
- `CutGraph` structure containing seam edges and topology information
"""
function compute_cut_graph(
    mesh::GeometryBasics.Mesh,
    singularities::Dict{Int,Float64};
    boundary_vertices::Set{Int}=Set{Int}()
)
    vs = coordinates(mesh)
    fs = faces(mesh)
    n_verts = length(vs)
    
    # Build adjacency graph (vertex-to-vertex connectivity)
    adj = Dict{Int,Vector{Int}}()
    for i in 1:n_verts
        adj[i] = Int[]
    end
    
    for f in fs
        idxs = Tuple(f)
        for i in 1:3
            v1 = idxs[i]
            v2 = idxs[mod1(i+1, 3)]
            push!(adj[v1], v2)
            push!(adj[v2], v1)
        end
    end
    
    # Remove duplicates in adjacency
    for v in 1:n_verts
        adj[v] = unique(adj[v])
    end
    
    seam_edges = Set{Tuple{Int,Int}}()
    singularity_tree = Dict{Int,Vector{Int}}()
    
    sing_vertices = collect(keys(singularities))
    
    if isempty(sing_vertices)
        # No singularities - for closed surface, still need to cut to make disk
        if isempty(boundary_vertices)
            println("  No singularities on closed surface - creating minimal cut")
            # TODO: Add minimal cut for closed surface without singularities
        end
        return CutGraph(seam_edges, singularity_tree, Vector{Int}[])
    end
    
    println("  Computing cut graph connecting $(length(sing_vertices)) singularities...")
    
    if !isempty(boundary_vertices)
        # Mesh has boundary: connect each singularity to nearest boundary
        seam_edges = connect_singularities_to_boundary(
            mesh, sing_vertices, boundary_vertices, adj
        )
    else
        # Closed surface: build tree connecting all singularities
        seam_edges, singularity_tree = build_singularity_tree(
            mesh, sing_vertices, adj
        )
    end
    
    println("  Cut graph computed: $(length(seam_edges)) seam edges")
    
    # TODO: Compute boundary loops after cutting
    boundary_loops = Vector{Int}[]
    
    return CutGraph(seam_edges, singularity_tree, boundary_loops)
end


"""
    connect_singularities_to_boundary(
        mesh, singularities, boundary_vertices, adj
    ) -> Set{Tuple{Int,Int}}

For meshes with boundary: connect each singularity to nearest boundary point
using shortest path (Dijkstra).
"""
function connect_singularities_to_boundary(
    mesh::GeometryBasics.Mesh,
    singularities::Vector{Int},
    boundary_vertices::Set{Int},
    adj::Dict{Int,Vector{Int}}
)
    vs = coordinates(mesh)
    seam_edges = Set{Tuple{Int,Int}}()
    
    for sing_v in singularities
        # Find shortest path from singularity to boundary using BFS
        path = shortest_path_to_boundary(sing_v, boundary_vertices, adj, vs)
        
        println("    Singularity $sing_v → boundary: path length = $(length(path))")
        
        # Add edges along path to seam
        for i in 1:(length(path)-1)
            v1, v2 = path[i], path[i+1]
            edge = v1 < v2 ? (v1, v2) : (v2, v1)
            push!(seam_edges, edge)
        end
    end
    
    return seam_edges
end


"""
    build_singularity_tree(mesh, singularities, adj)
    -> (Set{Tuple{Int,Int}}, Dict{Int,Vector{Int}})

For closed surfaces: build a tree connecting all singularities.
This is a spanning tree in the dual graph restricted to singular vertices.
"""
function build_singularity_tree(
    mesh::GeometryBasics.Mesh,
    singularities::Vector{Int},
    adj::Dict{Int,Vector{Int}}
)
    vs = coordinates(mesh)
    seam_edges = Set{Tuple{Int,Int}}()
    tree = Dict{Int,Vector{Int}}()
    
    if length(singularities) < 2
        # Single singularity - no tree needed
        return seam_edges, tree
    end
    
    # Initialize tree with first singularity as root
    for s in singularities
        tree[s] = Int[]
    end
    
    visited = Set{Int}([singularities[1]])
    
    # Build tree by connecting remaining singularities
    for i in 2:length(singularities)
        target = singularities[i]
        
        # Find shortest path from any visited singularity to target
        best_path = nothing
        best_length = Inf
        
        for source in visited
            path = shortest_path_between_vertices(source, target, adj, vs)
            if length(path) < best_length
                best_path = path
                best_length = length(path)
            end
        end
        
        if !isnothing(best_path) && length(best_path) >= 2
            # Add path edges to seam
            for j in 1:(length(best_path)-1)
                v1, v2 = best_path[j], best_path[j+1]
                edge = v1 < v2 ? (v1, v2) : (v2, v1)
                push!(seam_edges, edge)
            end
            
            # Update tree connectivity
            source = best_path[1]
            push!(tree[source], target)
            push!(tree[target], source)
            push!(visited, target)
        end
    end
    
    return seam_edges, tree
end


"""
    shortest_path_to_boundary(start, boundary, adj, coords) -> Vector{Int}

Find shortest path from start vertex to any boundary vertex using BFS.
Returns path as list of vertex indices.
"""
function shortest_path_to_boundary(
    start::Int,
    boundary::Set{Int},
    adj::Dict{Int,Vector{Int}},
    coords::Vector
)
    if start in boundary
        return [start]
    end
    
    # BFS to find shortest path
    queue = Queue{Int}()
    enqueue!(queue, start)
    parent = Dict{Int,Int}()
    parent[start] = -1
    
    target = -1
    while !isempty(queue)
        v = dequeue!(queue)
        
        if v in boundary
            target = v
            break
        end
        
        for neighbor in adj[v]
            if !haskey(parent, neighbor)
                parent[neighbor] = v
                enqueue!(queue, neighbor)
            end
        end
    end
    
    if target == -1
        # No path found (shouldn't happen)
        return [start]
    end
    
    # Reconstruct path
    path = Int[]
    current = target
    while current != -1
        pushfirst!(path, current)
        current = get(parent, current, -1)
    end
    
    return path
end


"""
    shortest_path_between_vertices(start, target, adj, coords) -> Vector{Int}

Find shortest path between two vertices using BFS.
"""
function shortest_path_between_vertices(
    start::Int,
    target::Int,
    adj::Dict{Int,Vector{Int}},
    coords::Vector
)
    if start == target
        return [start]
    end
    
    # BFS
    queue = Queue{Int}()
    enqueue!(queue, start)
    parent = Dict{Int,Int}()
    parent[start] = -1
    
    found = false
    while !isempty(queue)
        v = dequeue!(queue)
        
        if v == target
            found = true
            break
        end
        
        for neighbor in adj[v]
            if !haskey(parent, neighbor)
                parent[neighbor] = v
                enqueue!(queue, neighbor)
            end
        end
    end
    
    if !found
        # No path (disconnected graph?)
        return [start, target]
    end
    
    # Reconstruct path
    path = Int[]
    current = target
    while current != -1
        pushfirst!(path, current)
        current = get(parent, current, -1)
    end
    
    return path
end


"""
    identify_boundary_vertices(mesh::Mesh) -> Set{Int}

Identify vertices on the boundary of the mesh (vertices on edges that
belong to only one face).
"""
function identify_boundary_vertices(mesh::GeometryBasics.Mesh)
    fs = faces(mesh)
    
    # Count how many faces each edge belongs to
    edge_count = Dict{Tuple{Int,Int}, Int}()
    
    for f in fs
        idxs = Tuple(f)
        for i in 1:3
            v1 = idxs[i]
            v2 = idxs[mod1(i+1, 3)]
            edge = v1 < v2 ? (v1, v2) : (v2, v1)
            edge_count[edge] = get(edge_count, edge, 0) + 1
        end
    end
    
    # Boundary edges have count == 1
    boundary_verts = Set{Int}()
    for (edge, count) in edge_count
        if count == 1
            push!(boundary_verts, edge[1])
            push!(boundary_verts, edge[2])
        end
    end
    
    return boundary_verts
end


"""
    visualize_cut_graph(mesh, cut_graph, savepath=nothing)

Visualize the cut graph on the mesh (for debugging).
"""
function visualize_cut_graph(
    mesh::GeometryBasics.Mesh,
    cut_graph::CutGraph;
    savepath=nothing
)
    # TODO: Add visualization using plotters.jl
    println("Cut graph visualization:")
    println("  Seam edges: $(length(cut_graph.seam_edges))")
    println("  Singularity tree nodes: $(length(cut_graph.singularity_tree))")
end
