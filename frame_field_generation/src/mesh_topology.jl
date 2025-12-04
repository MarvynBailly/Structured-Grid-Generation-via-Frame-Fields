"""
    mesh_topology.jl

Utilities for building and managing triangular mesh topology, including:
- Edge maps
- Dual graph adjacency
- Face neighbors
- Reference edge computation
"""

using GeometryBasics
using LinearAlgebra

"""
    MeshTopology

Structure containing mesh topology information for frame field generation.

# Fields
- `edge_map::Dict{Tuple{Int,Int}, Int}` - Maps sorted (v1,v2) pairs to edge indices
- `dual_adj::Vector{Vector{Tuple{Int,Int}}}` - dual_adj[face_i] = [(neighbor_j, edge_k), ...]
- `edge_to_faces::Dict{Int, Vector{Int}}` - Maps edge index to incident face indices
- `n_edges::Int` - Total number of edges
- `reference_edges::Vector{Int}` - reference_edges[face_i] = edge index used as reference
- `kappa::Dict{Tuple{Int,Int}, Float64}` - Angle between reference edges across each edge
"""
struct MeshTopology
    edge_map::Dict{Tuple{Int,Int}, Int}
    dual_adj::Vector{Vector{Tuple{Int,Int}}}
    edge_to_faces::Dict{Int, Vector{Int}}
    n_edges::Int
    reference_edges::Vector{Int}
    kappa::Dict{Tuple{Int,Int}, Float64}
end

"""
    build_mesh_topology(mesh::Mesh) -> MeshTopology

Construct complete mesh topology from a triangular mesh.

# Arguments
- `mesh::Mesh` - Triangular mesh (GeometryBasics.Mesh)

# Returns
- `MeshTopology` - Complete topology structure

# Algorithm
1. Build edge map: enumerate all unique edges
2. Build dual adjacency: find face neighbors across each edge
3. Assign reference edges: choose one edge per face as x-axis reference
4. Compute κ_ij: angle between reference edges of adjacent faces
"""
function build_mesh_topology(mesh::Mesh)
    vs = coordinates(mesh)
    fs = faces(mesh)
    n_faces = length(fs)
    
    # Step 1: Build edge map
    edge_map = Dict{Tuple{Int,Int}, Int}()
    edge_counter = 0
    
    for face_idx in 1:n_faces
        face = fs[face_idx]
        # Get three edges of the triangle
        edges = [
            (min(face[1], face[2]), max(face[1], face[2])),
            (min(face[2], face[3]), max(face[2], face[3])),
            (min(face[3], face[1]), max(face[3], face[1]))
        ]
        
        for edge in edges
            if !haskey(edge_map, edge)
                edge_counter += 1
                edge_map[edge] = edge_counter
            end
        end
    end
    
    n_edges = edge_counter
    
    # Step 2: Build edge_to_faces map
    edge_to_faces = Dict{Int, Vector{Int}}()
    for e_idx in 1:n_edges
        edge_to_faces[e_idx] = Int[]
    end
    
    for face_idx in 1:n_faces
        face = fs[face_idx]
        edges = [
            (min(face[1], face[2]), max(face[1], face[2])),
            (min(face[2], face[3]), max(face[2], face[3])),
            (min(face[3], face[1]), max(face[3], face[1]))
        ]
        
        for edge in edges
            e_idx = edge_map[edge]
            push!(edge_to_faces[e_idx], face_idx)
        end
    end
    
    # Step 3: Build dual adjacency (face neighbors)
    dual_adj = [Vector{Tuple{Int,Int}}() for _ in 1:n_faces]
    
    for (edge, e_idx) in edge_map
        incident_faces = edge_to_faces[e_idx]
        if length(incident_faces) == 2
            f1, f2 = incident_faces
            push!(dual_adj[f1], (f2, e_idx))
            push!(dual_adj[f2], (f1, e_idx))
        end
    end
    
    # Step 4: Assign reference edges (choose first edge of each face)
    reference_edges = zeros(Int, n_faces)
    
    for face_idx in 1:n_faces
        face = fs[face_idx]
        # Use first edge (v1, v2) as reference
        edge = (min(face[1], face[2]), max(face[1], face[2]))
        reference_edges[face_idx] = edge_map[edge]
    end
    
    # Step 5: Compute κ_ij (angle between reference edges)
    kappa = Dict{Tuple{Int,Int}, Float64}()
    
    # Create reverse map from edge index to edge tuple
    edge_idx_to_tuple = Dict{Int, Tuple{Int,Int}}()
    for (edge_tuple, idx) in edge_map
        edge_idx_to_tuple[idx] = edge_tuple
    end
    
    for face_idx in 1:n_faces
        for (neighbor_idx, edge_idx) in dual_adj[face_idx]
            # Get reference edge directions
            ref_edge_i = reference_edges[face_idx]
            ref_edge_j = reference_edges[neighbor_idx]
            
            # Find edge tuples
            edge_i = edge_idx_to_tuple[ref_edge_i]
            edge_j = edge_idx_to_tuple[ref_edge_j]
            
            # Get edge directions (as 2D vectors on the mesh)
            v1_i, v2_i = edge_i
            dir_i = vs[v2_i] - vs[v1_i]
            dir_i_2d = [dir_i[1], dir_i[2]]  # Project to xy-plane
            dir_i_2d = dir_i_2d / norm(dir_i_2d)
            
            v1_j, v2_j = edge_j
            dir_j = vs[v2_j] - vs[v1_j]
            dir_j_2d = [dir_j[1], dir_j[2]]
            dir_j_2d = dir_j_2d / norm(dir_j_2d)
            
            # Compute angle between directions using atan2
            # κ_ij ∈ (-π, π]
            angle = atan(dir_j_2d[2], dir_j_2d[1]) - atan(dir_i_2d[2], dir_i_2d[1])
            
            # Normalize to (-π, π]
            while angle > π
                angle -= 2π
            end
            while angle <= -π
                angle += 2π
            end
            
            kappa[(face_idx, neighbor_idx)] = angle
        end
    end
    
    return MeshTopology(edge_map, dual_adj, edge_to_faces, n_edges, reference_edges, kappa)
end

"""
    get_face_neighbors(topology::MeshTopology, face_idx::Int) -> Vector{Int}

Get all neighboring face indices for a given face.
"""
function get_face_neighbors(topology::MeshTopology, face_idx::Int)
    return [neighbor for (neighbor, _) in topology.dual_adj[face_idx]]
end

"""
    get_shared_edge(topology::MeshTopology, face_i::Int, face_j::Int) -> Union{Int, Nothing}

Get the edge index shared between two faces, or nothing if not adjacent.
"""
function get_shared_edge(topology::MeshTopology, face_i::Int, face_j::Int)
    for (neighbor, edge_idx) in topology.dual_adj[face_i]
        if neighbor == face_j
            return edge_idx
        end
    end
    return nothing
end

"""
    compute_angle_defect(mesh::Mesh, vertex_idx::Int) -> Float64

Compute the angle defect at a vertex: 2π - sum of incident face angles.

For flat meshes, this should be approximately zero for interior vertices.
"""
function compute_angle_defect(mesh::Mesh, vertex_idx::Int)
    vs = coordinates(mesh)
    fs = faces(mesh)
    
    # Find all faces incident to this vertex
    incident_faces = Int[]
    for (f_idx, face) in enumerate(fs)
        if vertex_idx in face
            push!(incident_faces, f_idx)
        end
    end
    
    # Sum angles at this vertex in each incident face
    total_angle = 0.0
    
    for f_idx in incident_faces
        face = fs[f_idx]
        # Find position of vertex in face
        local_idx = findfirst(==(vertex_idx), face)
        
        # Get three vertices of triangle
        v1 = vs[face[local_idx]]
        v2 = vs[face[mod1(local_idx + 1, 3)]]
        v3 = vs[face[mod1(local_idx - 1, 3)]]
        
        # Compute angle at v1
        edge1 = v2 - v1
        edge2 = v3 - v1
        
        cos_angle = dot(edge1, edge2) / (norm(edge1) * norm(edge2))
        cos_angle = clamp(cos_angle, -1.0, 1.0)
        angle = acos(cos_angle)
        
        total_angle += angle
    end
    
    # Angle defect = 2π - total_angle
    return 2π - total_angle
end
