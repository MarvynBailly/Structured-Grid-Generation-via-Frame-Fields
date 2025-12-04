"""
    dijkstra_forest.jl

Implementation of Dijkstra forest construction for fixing period jump variables.

The Dijkstra forest is used to identify edges whose period jumps can be fixed to zero
without changing the minimization energy. This reduces the dimensionality of the
mixed-integer problem.
"""

using DataStructures
using GeometryBasics

"""
    DijkstraForest

Structure representing a Dijkstra forest over the dual graph.

# Fields
- `tree_edges::Set{Int}` - Set of edge indices in the spanning forest
- `fixable_edges::Set{Int}` - Edges that CAN be fixed to zero (tree edges)
- `fixed_edges_per_face::Dict{Int, Int}` - For each face, which edge is fixed to zero
- `free_edges::Vector{Int}` - Edges whose period jumps are optimization variables
- `parent::Dict{Int, Int}` - parent[face] = parent face in tree (or -1 for roots)
"""
struct DijkstraForest
    tree_edges::Set{Int}
    fixable_edges::Set{Int}
    fixed_edges_per_face::Dict{Int, Int}
    free_edges::Vector{Int}
    parent::Dict{Int, Int}
end

"""
    build_dijkstra_forest(
        mesh::Mesh,
        topology::MeshTopology,
        constrained_faces::Vector{Int}=Int[]
    ) -> DijkstraForest

Build a Dijkstra spanning forest over the dual graph of the mesh, starting from
constrained faces. Then for each triangle, fix one of its edges to have p=0.

# Arguments
- `mesh::Mesh` - Triangular mesh
- `topology::MeshTopology` - Mesh topology structure
- `constrained_faces::Vector{Int}` - Indices of faces with fixed angles (roots of forest)

# Returns
- `DijkstraForest` - Forest structure with fixable edges and per-face fixed edges

# Algorithm
1. Run Dijkstra's algorithm to identify tree edges (edges that CAN be fixed)
2. For each triangle, select one of its edges to fix to zero:
   - Prefer tree edges (fixable without changing energy)
   - This ensures one edge per triangle is fixed, reducing problem size

# Note
If no constrained faces are provided, starts from face 1 as root.
"""
function build_dijkstra_forest(
    mesh::Mesh,
    topology::MeshTopology,
    constrained_faces::Vector{Int}=Int[]
)
    fs = faces(mesh)
    n_faces = length(fs)
    
    # Initialize
    visited = falses(n_faces)
    parent = Dict{Int, Int}()
    tree_edges = Set{Int}()
    
    # Priority queue: (distance, face_idx)
    pq = PriorityQueue{Int, Float64}()
    
    # If no constrained faces, start from face 1
    if isempty(constrained_faces)
        constrained_faces = [1]
    end
    
    # Initialize roots
    for root in constrained_faces
        pq[root] = 0.0
        parent[root] = -1  # Root has no parent
    end
    
    # Dijkstra's algorithm on dual graph
    while !isempty(pq)
        current_face = dequeue!(pq)
        
        if visited[current_face]
            continue
        end
        
        visited[current_face] = true
        
        # Process neighbors in dual graph
        for (neighbor_face, edge_idx) in topology.dual_adj[current_face]
            if !visited[neighbor_face]
                # Add edge to tree
                push!(tree_edges, edge_idx)
                parent[neighbor_face] = current_face
                
                # Add neighbor to queue with distance (can use uniform weights)
                if !haskey(pq, neighbor_face)
                    pq[neighbor_face] = length(tree_edges)  # Simple distance metric
                end
            end
        end
    end
    
    # Fixable edges are tree edges
    fixable_edges = tree_edges
    
    # For each triangle, fix one edge to zero
    # Prefer edges in the tree (fixable without changing energy)
    # Ensure each edge is fixed by at most one face
    fixed_edges_per_face = Dict{Int, Int}()
    already_fixed = Set{Int}()
    
    for face_idx in 1:n_faces
        face = fs[face_idx]
        
        # Get all three edges of this triangle
        edges_of_face = [
            (min(face[1], face[2]), max(face[1], face[2])),
            (min(face[2], face[3]), max(face[2], face[3])),
            (min(face[3], face[1]), max(face[3], face[1]))
        ]
        
        # Convert to edge indices
        edge_indices = [topology.edge_map[e] for e in edges_of_face]
        
        # Find edges that are available (not already fixed by another face)
        available_edges = [e_idx for e_idx in edge_indices if !(e_idx in already_fixed)]
        
        # Among available edges, prefer fixable (tree) edges
        fixable_available = [e_idx for e_idx in available_edges if e_idx in fixable_edges]
        
        # Pick one edge to fix
        chosen_edge = if !isempty(fixable_available)
            fixable_available[1]  # Prefer fixable edges
        elseif !isempty(available_edges)
            available_edges[1]  # Otherwise pick any available edge
        else
            # All edges already fixed (shouldn't happen in valid mesh)
            edge_indices[1]  # Fallback
        end
        
        fixed_edges_per_face[face_idx] = chosen_edge
        push!(already_fixed, chosen_edge)
    end
    
    # Now determine which edges are truly free (not fixed by any face)
    all_fixed = Set(values(fixed_edges_per_face))
    all_edges = Set(1:topology.n_edges)
    free_edges = collect(setdiff(all_edges, all_fixed))
    
    return DijkstraForest(tree_edges, fixable_edges, fixed_edges_per_face, free_edges, parent)
end

"""
    is_tree_edge(forest::DijkstraForest, edge_idx::Int) -> Bool

Check if an edge is in the spanning forest.
"""
function is_tree_edge(forest::DijkstraForest, edge_idx::Int)
    return edge_idx in forest.tree_edges
end

"""
    is_fixed_edge(forest::DijkstraForest, edge_idx::Int) -> Bool

Check if an edge is fixed to p=0 by some face.
"""
function is_fixed_edge(forest::DijkstraForest, edge_idx::Int)
    return edge_idx in values(forest.fixed_edges_per_face)
end

"""
    get_fixed_edge_for_face(forest::DijkstraForest, face_idx::Int) -> Int

Get the edge that is fixed to p=0 for a given face.
"""
function get_fixed_edge_for_face(forest::DijkstraForest, face_idx::Int)
    return forest.fixed_edges_per_face[face_idx]
end

"""
    get_free_variables(forest::DijkstraForest) -> Vector{Int}

Get the list of edge indices whose period jumps are free variables.
"""
function get_free_variables(forest::DijkstraForest)
    return forest.free_edges
end

"""
    count_free_edges(forest::DijkstraForest) -> Int

Count the number of free edge variables.
"""
function count_free_edges(forest::DijkstraForest)
    return length(forest.free_edges)
end
