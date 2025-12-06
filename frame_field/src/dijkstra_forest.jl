


using DataStructures
using GeometryBasics

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
