using GeometryBasics
using DataStructures
using LinearAlgebra
using Statistics
using WGLMakie  # Interactive 3D plots in browser
using FileIO
using MeshIO

# Resolve Mesh name conflict by explicitly using GeometryBasics version
const Mesh = GeometryBasics.Mesh

include("meshIO.jl")
include("plotters.jl")


"""
    build_topology(mesh::Mesh)

Constructs the dual graph adjacency and unique edge numbering.
"""
function build_topology(mesh::Mesh)
    faces_list = decompose(TriangleFace{Int}, mesh)
    edge_map = Dict{Tuple{Int, Int}, Int}()
    # Temporary storage to build dual adjacency: Edge Index -> List of Faces
    edge_to_faces = Dict{Int, Vector{Int}}()
    
    edge_counter = 0
    
    # 1. Identify all unique edges and assign IDs
    for (f_idx, face) in enumerate(faces_list)
        # Edges of the triangle: (1,2), (2,3), (3,1)
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
        # Interior edges connect exactly 2 faces. Boundary edges connect 1.
        # We only care about traversing interior edges for the dual graph.
        if length(connected_faces) == 2
            f1, f2 = connected_faces
            push!(dual_adj[f1], (f2, e_idx))
            push!(dual_adj[f2], (f1, e_idx))
        end
    end
    
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




#####################################################
# 1. Load the mesh (Gmsh .msh file, or .obj, .stl)
filename = "simple-square.msh"  
# filename = "unit-square.msh"  
# filename = "teapot.msh"

filepath = joinpath("SimpleAutoQuad", "output", "meshes", filename)

# This handles all the dirty work of file parsing and type conversion
mesh_info = load_triangulation(filepath)

# 2. Define constraints (e.g., hardcode specific faces or use a selection heuristic)
# for now, let's just pick face 1.
# constrained_faces = [1,5, 300 ,500, 240, 400] 
# randomly pick n constrained faces
n = 2
total_faces = length(faces(mesh_info))
constrained_faces = rand(1:total_faces, n)

# 3. Build the topology (Uses the code from Step 1)
topo = build_topology(mesh_info)

# 4. Compute the Cross Field Spanning Forest
fixed_edges = compute_spanning_forest(topo, length(faces(mesh_info)), constrained_faces)

# 5. Visualize 
fig = visualize_forest(mesh_info, topo, fixed_edges, constrained_faces)