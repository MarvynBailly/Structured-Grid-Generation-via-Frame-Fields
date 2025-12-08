


using DataStructures
using GeometryBasics
using LinearAlgebra

"""
    identify_boundary_faces_and_angles(mesh::GeometryBasics.Mesh) -> (Vector{Int}, Vector{Float64})

Identify all boundary faces and compute their target angles based on boundary tangent directions.

For each face on the boundary, computes the angle of the boundary edge tangent.
This angle will be used as a hard constraint to align the frame field with the boundary.

# Returns
- `boundary_faces::Vector{Int}` - Indices of faces incident to boundary edges
- `boundary_angles::Vector{Float64}` - Target angles for each boundary face (in radians)
"""
function identify_boundary_faces_and_angles(mesh::GeometryBasics.Mesh)
    faces_list = decompose(TriangleFace{Int}, mesh)
    vs = coordinates(mesh)
    n_faces = length(faces_list)
    
    # Build edge to faces mapping
    edge_to_faces = Dict{Tuple{Int,Int}, Vector{Int}}()
    face_to_edges = Dict{Int, Vector{Tuple{Int,Int}}}()
    
    for (f_idx, face) in enumerate(faces_list)
        verts = (face[1], face[2], face[3])
        face_edges = Tuple{Int,Int}[]
        
        for k in 1:3
            v1, v2 = verts[k], verts[mod1(k+1, 3)]
            edge = v1 < v2 ? (v1, v2) : (v2, v1)
            push!(face_edges, edge)
            
            if !haskey(edge_to_faces, edge)
                edge_to_faces[edge] = Int[]
            end
            push!(edge_to_faces[edge], f_idx)
        end
        
        face_to_edges[f_idx] = face_edges
    end
    
    # Find boundary edges (edges with only one incident face)
    boundary_edges = Set{Tuple{Int,Int}}()
    for (edge, incident_faces) in edge_to_faces
        if length(incident_faces) == 1
            push!(boundary_edges, edge)
        end
    end
    
    if isempty(boundary_edges)
        return Int[], Float64[]
    end
    
    # For each boundary face, compute the angle of its boundary edge
    boundary_face_angles = Dict{Int, Float64}()
    
    for (f_idx, edges) in face_to_edges
        # Find which edge(s) of this face are on the boundary
        boundary_edge_of_face = nothing
        for edge in edges
            if edge in boundary_edges
                boundary_edge_of_face = edge
                break
            end
        end
        
        if boundary_edge_of_face !== nothing
            # Compute tangent direction of boundary edge
            v1_idx, v2_idx = boundary_edge_of_face
            p1 = vs[v1_idx]
            p2 = vs[v2_idx]
            
            # Tangent vector (in 2D projection)
            tangent = Point2f(p2[1] - p1[1], p2[2] - p1[2])
            tangent_normalized = normalize(tangent)
            
            # Compute angle of tangent
            angle = atan(tangent_normalized[2], tangent_normalized[1])
            
            # Store this as the constraint for this boundary face
            boundary_face_angles[f_idx] = angle
        end
    end
    
    # Convert to vectors
    boundary_faces = collect(keys(boundary_face_angles))
    boundary_angles = [boundary_face_angles[f] for f in boundary_faces]
    
    return boundary_faces, boundary_angles
end

"""
    compute_spanning_forest(mesh::GeometryBasics.Mesh; constrained_faces=nothing)

Builds the mesh dual adjacency and computes a spanning forest (BFS/Dijkstra with unit
weights) starting from `constrained_faces`. Returns a Set of edge tuples (vmin, vmax)
that form the forest (these are the edges that can be fixed).

If `mesh` is omitted, the function will try to use a global `mesh` variable from `Main`.
If `constrained_faces` is `nothing`, the function will try to use `Main.constrained_faces`.
"""
function compute_spanning_forest(mesh::GeometryBasics.Mesh; constrained_faces=nothing)
    # Build edges and dual adjacency
    faces_list = decompose(TriangleFace{Int}, mesh)
    edge_map = Dict{Tuple{Int,Int}, Int}()
    edge_to_faces = Dict{Int, Vector{Int}}()
    edge_counter = 0

    for (f_idx, face) in enumerate(faces_list)
        vs = (face[1], face[2], face[3])
        for k in 1:3
            a, b = minmax(vs[k], vs[mod1(k+1, 3)])
            key = (a,b)
            if !haskey(edge_map, key)
                edge_counter += 1
                edge_map[key] = edge_counter
                edge_to_faces[edge_counter] = Int[]
            end
            push!(edge_to_faces[edge_map[key]], f_idx)
        end
    end

    # Build dual adjacency: face -> [(neighbor_face, edge_index), ...]
    n_faces = length(faces_list)
    dual_adj = [Vector{Tuple{Int,Int}}() for _ in 1:n_faces]
    for (e_idx, conn) in edge_to_faces
        if length(conn) == 2
            f1, f2 = conn[1], conn[2]
            push!(dual_adj[f1], (f2, e_idx))
            push!(dual_adj[f2], (f1, e_idx))
        end
    end

    # BFS queue
    visited = falses(n_faces)
    forest_edges_idx = Set{Int}()
    q = Queue{Int}()

    # initialize roots
    for f in constrained_faces
        if 1 <= f <= n_faces
            enqueue!(q, f)
            visited[f] = true
        end
    end
    if isempty(constrained_faces) && n_faces > 0
        enqueue!(q, 1)
        visited[1] = true
    end

    while !isempty(q)
        cur = dequeue!(q)
        for (nbr, eidx) in dual_adj[cur]
            if !visited[nbr]
                push!(forest_edges_idx, eidx)
                visited[nbr] = true
                enqueue!(q, nbr)
            end
        end
    end

    # Convert edge indices back to vertex tuple keys
    rev_map = Dict{Int, Tuple{Int,Int}}()
    for (k,v) in edge_map
        rev_map[v] = k
    end

    forest_edges = Set{Tuple{Int,Int}}()
    for eidx in forest_edges_idx
        push!(forest_edges, rev_map[eidx])
    end

    return forest_edges
end


"""
    fix_suitable_edges(mesh::GeometryBasics.Mesh, potential_fixed_edges)
From the set of `potential_fixed_edges`, selects one edge per face to be fixed. 
Returns a Dict{Int, Tuple{Int,Int}} mapping face_idx → edge.
"""
function fix_suitable_edges(mesh::GeometryBasics.Mesh, potential_fixed_edges)
    faces_list = decompose(TriangleFace{Int}, mesh)
    fixed_edges_per_face = Dict{Int, Tuple{Int,Int}}()

    for (f_idx, face) in enumerate(faces_list)
        vs = (face[1], face[2], face[3])
        # find edges of this face
        face_edges = Tuple{Int, Int}[]
        for k in 1:3
            v_a, v_b = minmax(vs[k], vs[mod1(k+1, 3)])
            push!(face_edges, (v_a, v_b))
        end
        # Randomize selection: skip some faces to ensure proper subset
        if rand() < 0.5
            continue
        end
        # find the first edge that is in potential_fixed_edges and fix it
        for edge in face_edges
            if edge in potential_fixed_edges
                fixed_edges_per_face[f_idx] = edge
                break
            end
        end
    end
    return fixed_edges_per_face
end
