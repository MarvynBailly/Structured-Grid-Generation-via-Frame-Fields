"""
    mesh_types.jl

Data structures for triangular meshes and frame field representations.
"""

using GeometryBasics
using SparseArrays

"""
    TriangularMeshTopology

Stores topological information about a triangular mesh.

# Fields
- `mesh::GeometryBasics.Mesh` - GeometryBasics mesh with vertices and faces
- `n_faces::Int` - Number of triangle faces
- `n_edges::Int` - Number of edges
- `edge_map::Dict{Tuple{Int,Int}, Int}` - Maps vertex pair (v1,v2) to edge index
- `edge_to_faces::Dict{Int, Vector{Int}}` - Maps edge index to adjacent face indices
- `dual_adj::Dict{Int, Vector{Tuple{Int, Tuple{Int,Int}}}}` - Dual graph: face → [(neighbor_face, edge)]
- `face_reference_edges::Dict{Int, Tuple{Int,Int}}` - Maps face index to its reference edge (v1,v2)
"""
struct TriangularMeshTopology
    mesh::GeometryBasics.Mesh
    n_faces::Int
    n_edges::Int
    edge_map::Dict{Tuple{Int,Int}, Int}
    edge_to_faces::Dict{Int, Vector{Int}}
    dual_adj::Dict{Int, Vector{Tuple{Int, Tuple{Int,Int}}}}
    face_reference_edges::Dict{Int, Tuple{Int,Int}}
end

"""
    build_mesh_topology(mesh::GeometryBasics.Mesh) -> TriangularMeshTopology

Construct topological information from a triangular mesh.

# Arguments
- `mesh` - GeometryBasics mesh containing vertices and triangle faces

# Returns
- `TriangularMeshTopology` with all topological data structures built
"""
function build_mesh_topology(mesh::GeometryBasics.Mesh)
    fs = faces(mesh)
    n_faces = length(fs)
    
    # Initialize data structures
    edge_map = Dict{Tuple{Int,Int}, Int}()
    edge_to_faces = Dict{Int, Vector{Int}}()
    dual_adj = Dict{Int, Vector{Tuple{Int, Tuple{Int,Int}}}}()
    face_reference_edges = Dict{Int, Tuple{Int,Int}}()
    
    for i in 1:n_faces
        dual_adj[i] = []
    end
    
    edge_idx = 0
    
    # Build edge map and face adjacency
    for (i, face_i) in enumerate(fs)
        verts_i = [face_i[1], face_i[2], face_i[3]]
        
        # Set reference edge as first edge of triangle (v1, v2)
        ref_edge = (verts_i[1], verts_i[2])#verts_i[1] < verts_i[2] ? (verts_i[1], verts_i[2]) : (verts_i[2], verts_i[1])
        face_reference_edges[i] = ref_edge
        
        for (j, face_j) in enumerate(fs)
            if i >= j
                continue
            end
            
            verts_j = [face_j[1], face_j[2], face_j[3]]
            shared = intersect(verts_i, verts_j)
            
            if length(shared) == 2
                v1, v2 = shared[1], shared[2]
                edge = (v1, v2)  # Use vertices as given, no sorting
                
                # Check both orientations since edge might be stored either way
                edge_canonical = haskey(edge_map, edge) ? edge : (v2, v1)
                
                if !haskey(edge_map, edge) && !haskey(edge_map, (v2, v1))
                    edge_idx += 1
                    edge_map[edge] = edge_idx
                    edge_to_faces[edge_idx] = [i, j]
                    
                    # Add to dual adjacency
                    push!(dual_adj[i], (j, edge))
                    push!(dual_adj[j], (i, edge))
                end
            end
        end
    end
    
    n_edges = edge_idx
    
    return TriangularMeshTopology(
        mesh,
        n_faces,
        n_edges,
        edge_map,
        edge_to_faces,
        dual_adj,
        face_reference_edges
    )
end

"""
    CrossField

Stores a frame field solution on a triangular mesh.

# Fields
- `topology::TriangularMeshTopology` - Underlying mesh topology
- `theta::Vector{Float64}` - Frame angles for each face (in radians)
- `period_jumps::Dict{Tuple{Int,Int}, Int}` - Period jumps p_ij for each (face_i, face_j) pair
    * Only canonical direction (i < j) stored as independent variable
    * Reverse direction computed as p_ji = -p_ij for convenience
- `transport_angles::Dict{Tuple{Int,Int}, Float64}` - Transport angles κ_ij between adjacent faces
    * Only canonical direction (i < j) stored
    * Reverse direction computed as κ_ji = -κ_ij for convenience
- `constrained_faces::Vector{Int}` - Indices of faces with fixed angles
- `constrained_angles::Vector{Float64}` - Fixed angle values for constrained faces
- `fixed_edges::Set{Tuple{Int,Int}}` - Edges where period jumps are fixed to zero
"""
mutable struct CrossField
    topology::TriangularMeshTopology
    theta::Vector{Float64}
    period_jumps::Dict{Tuple{Int,Int}, Int}
    transport_angles::Dict{Tuple{Int,Int}, Float64}
    constrained_faces::Vector{Int}
    constrained_angles::Vector{Float64}
    fixed_edges::Set{Tuple{Int,Int}}
end

"""
    CrossField(
        topology::TriangularMeshTopology;
        constrained_faces::Vector{Int}=Int[],
        constrained_angles::Vector{Float64}=Float64[],
        fixed_edges::Set{Tuple{Int,Int}}=Set{Tuple{Int,Int}}()
    ) -> CrossField

Create a new CrossField with initialized values.

# Arguments
- `topology` - Triangular mesh topology
- `constrained_faces` - (Optional) Face indices with fixed angles
- `constrained_angles` - (Optional) Fixed angle values
- `fixed_edges` - (Optional) Set of edges with p_ij = 0

# Returns
- `CrossField` with theta initialized to zeros and period jumps to zeros
"""
function CrossField(
    topology::TriangularMeshTopology;
    constrained_faces::Vector{Int}=Int[],
    constrained_angles::Vector{Float64}=Float64[],
    fixed_edges::Set{Tuple{Int,Int}}=Set{Tuple{Int,Int}}()
)
    n_faces = topology.n_faces
    
    # Initialize theta to zeros (or constrained values where applicable)
    theta = zeros(Float64, n_faces)
    for (i, face_idx) in enumerate(constrained_faces)
        theta[face_idx] = constrained_angles[i]
    end
    
    # Initialize period jumps to zeros (only canonical direction i < j stored)
    period_jumps = Dict{Tuple{Int,Int}, Int}()
    for face_i in 1:n_faces
        for (face_j, edge) in topology.dual_adj[face_i]
            if face_i < face_j
                period_jumps[(face_i, face_j)] = 0
            end
        end
    end
    
    # Compute transport angles
    transport_angles = compute_transport_angles(topology)
    
    return CrossField(
        topology,
        theta,
        period_jumps,
        transport_angles,
        constrained_faces,
        constrained_angles,
        fixed_edges
    )
end

"""
    compute_transport_angles(topology::TriangularMeshTopology) -> Dict{Tuple{Int,Int}, Float64}

Compute transport angles κ_ij for all adjacent face pairs.

Transport angle represents the rotation between local coordinate systems of adjacent faces.
We measure the angle of the shared edge in each face's local coordinate system, 
and the difference gives the rotation needed to align the frames.

Algorithm:
1. For each pair of adjacent faces (i, j) sharing an edge
2. Get the shared edge vector (v1, v2)
3. Measure angle of edge in face i's local frame (using reference edge as x-axis)
4. Measure angle of edge in face j's local frame (using reference edge as x-axis)
5. κ_ij = angle_j - angle_i (adjusted to (-π, π])

# Arguments
- `topology` - Triangular mesh topology

# Returns
- Dictionary mapping (face_i, face_j) → κ_ij (transport angle in radians)
"""
function compute_transport_angles(topology::TriangularMeshTopology)
    fs = faces(topology.mesh)
    vs = coordinates(topology.mesh)
    n_faces = topology.n_faces
    
    kappa = Dict{Tuple{Int,Int}, Float64}()
    
    for face_i in 1:n_faces
        for (face_j, shared_edge) in topology.dual_adj[face_i]
            if face_i < face_j
                # Get reference edges for both faces
                ref_i_v1, ref_i_v2 = topology.face_reference_edges[face_i]
                ref_j_v1, ref_j_v2 = topology.face_reference_edges[face_j]
                
                # Get shared edge vector
                shared_p1 = vs[shared_edge[1]]
                shared_p2 = vs[shared_edge[2]]
                edge_vec = (shared_p2[1] - shared_p1[1], shared_p2[2] - shared_p1[2])
                
                # Compute local coordinate systems (x-axis along reference edge)
                # For face i:
                ref_i_p1 = vs[ref_i_v1]
                ref_i_p2 = vs[ref_i_v2]
                xi = (ref_i_p2[1] - ref_i_p1[1], ref_i_p2[2] - ref_i_p1[2])
                xi_norm = sqrt(xi[1]^2 + xi[2]^2)
                xi = (xi[1] / xi_norm, xi[2] / xi_norm)  # Normalized x-axis
                yi = (-xi[2], xi[1])  # Perpendicular y-axis (90° CCW rotation)
                
                # For face j:
                ref_j_p1 = vs[ref_j_v1]
                ref_j_p2 = vs[ref_j_v2]
                xj = (ref_j_p2[1] - ref_j_p1[1], ref_j_p2[2] - ref_j_p1[2])
                xj_norm = sqrt(xj[1]^2 + xj[2]^2)
                xj = (xj[1] / xj_norm, xj[2] / xj_norm)  # Normalized x-axis
                yj = (-xj[2], xj[1])  # Perpendicular y-axis (90° CCW rotation)
                
                # Normalize edge vector
                edge_norm = sqrt(edge_vec[1]^2 + edge_vec[2]^2)
                edge_vec = (edge_vec[1] / edge_norm, edge_vec[2] / edge_norm)
                
                # Angle of edge in frame i
                angle_i = atan(edge_vec[1] * yi[1] + edge_vec[2] * yi[2],  # dot(edge, yi)
                              edge_vec[1] * xi[1] + edge_vec[2] * xi[2])   # dot(edge, xi)
                
                # Angle of edge in frame j
                angle_j = atan(edge_vec[1] * yj[1] + edge_vec[2] * yj[2],  # dot(edge, yj)
                              edge_vec[1] * xj[1] + edge_vec[2] * xj[2])   # dot(edge, xj)
                
                # Transport angle: rotation from frame i to frame j
                # Since we measure the same edge in both frames, the rotation is angle_i - angle_j
                kappa_ij = angle_i - angle_j
                
                # Wrap to (-π, π]
                kappa_ij = mod2pi(kappa_ij + π) - π
                
                println("For faces $face_i and $face_j sharing edge $shared_edge:")
                println("  Angle of shared edge in frame i: $(rad2deg(angle_i))°")
                println("  Angle of shared edge in frame j: $(rad2deg(angle_j))°")
                println("  κ_ij = angle_i - angle_j = $(rad2deg(kappa_ij))°\n")
                
                # Only store canonical direction (i < j)
                kappa[(face_i, face_j)] = kappa_ij
            end
        end
    end
    
    return kappa
end

"""
    get_transport_angle(field::CrossField, face_i::Int, face_j::Int) -> Float64

Get the transport angle between two adjacent faces.
Handles antisymmetry automatically: κ_ji = -κ_ij

# Arguments
- `field` - CrossField containing transport angles
- `face_i` - First face index
- `face_j` - Second face index

# Returns
- Transport angle value κ_ij in radians
"""
function get_transport_angle(field::CrossField, face_i::Int, face_j::Int)
    if face_i < face_j
        return field.transport_angles[(face_i, face_j)]
    else
        return -field.transport_angles[(face_j, face_i)]
    end
end

"""
    get_period_jump(field::CrossField, face_i::Int, face_j::Int) -> Int

Get the period jump between two adjacent faces.
Handles antisymmetry automatically: p_ji = -p_ij

# Arguments
- `field` - CrossField containing period jumps
- `face_i` - First face index
- `face_j` - Second face index

# Returns
- Period jump value p_ij
"""
function get_period_jump(field::CrossField, face_i::Int, face_j::Int)
    if face_i < face_j
        return field.period_jumps[(face_i, face_j)]
    else
        return -field.period_jumps[(face_j, face_i)]
    end
end

"""
    set_period_jump!(field::CrossField, face_i::Int, face_j::Int, value::Int)

Set the period jump between two adjacent faces.
Only stores canonical direction (i < j).

# Arguments
- `field` - CrossField to modify
- `face_i` - First face index
- `face_j` - Second face index
- `value` - Period jump value p_ij
"""
function set_period_jump!(field::CrossField, face_i::Int, face_j::Int, value::Int)
    if face_i < face_j
        field.period_jumps[(face_i, face_j)] = value
    else
        field.period_jumps[(face_j, face_i)] = -value
    end
end

"""
    compute_smoothness_energy(field::CrossField) -> Float64

Compute the frame field smoothness energy.

E = Σ (θ_i + κ_ij + (π/2)p_ij - θ_j)²

# Arguments
- `field` - CrossField containing theta values and period jumps

# Returns
- Energy value (sum of squared misalignments)
"""
function compute_smoothness_energy(field::CrossField)
    energy = 0.0
    topology = field.topology
    n_faces = topology.n_faces
    
    # Sum over all edges (count each once)
    for face_i in 1:n_faces
        for (face_j, edge) in topology.dual_adj[face_i]
            if face_i < face_j  # Count each edge once
                kappa_ij = get_transport_angle(field, face_i, face_j)
                p_ij = get_period_jump(field, face_i, face_j)
                
                diff = field.theta[face_i] + kappa_ij + (π/2) * p_ij - field.theta[face_j]
                energy += diff^2
            end
        end
    end
    
    return energy
end

"""
    get_reference_angle(field::CrossField, face_idx::Int) -> Float64

Compute the angle of the reference edge for a given face.

# Arguments
- `field` - CrossField containing the topology
- `face_idx` - Index of the face

# Returns
- Angle in radians of the reference edge direction
"""
function get_reference_angle(field::CrossField, face_idx::Int)
    vs = coordinates(field.topology.mesh)
    ref_edge = field.topology.face_reference_edges[face_idx]
    v1, v2 = ref_edge
    
    # Get 2D positions (ignore z-coordinate if present)
    p1 = vs[v1]
    p2 = vs[v2]
    
    # Compute angle of reference edge
    ref_vec = (p2[1] - p1[1], p2[2] - p1[2])
    return atan(ref_vec[2], ref_vec[1])
end

"""
    get_absolute_angle(field::CrossField, face_idx::Int) -> Float64

Get the absolute frame angle for a given face.
This is the reference edge angle plus theta.

# Arguments
- `field` - CrossField containing the solution
- `face_idx` - Index of the face

# Returns
- Absolute angle in radians (reference_angle + theta)
"""
function get_absolute_angle(field::CrossField, face_idx::Int)
    ref_angle = get_reference_angle(field, face_idx)
    return ref_angle + field.theta[face_idx]
end

"""
    get_frame_directions(field::CrossField, face_idx::Int) -> Vector{Vector{Float64}}

Get the 4 orthogonal frame directions for a given face in absolute coordinates.
When theta=0, the first direction aligns with the reference edge.

# Arguments
- `field` - CrossField containing the solution
- `face_idx` - Index of the face

# Returns
- Vector of 4 direction vectors (2D), each rotated by π/2 from the previous
"""
function get_frame_directions(field::CrossField, face_idx::Int)
    # Get absolute angle (reference angle + theta)
    absolute_angle = get_absolute_angle(field, face_idx)
    
    # Generate 4 directions at angles: angle, angle + π/2, angle + π, angle + 3π/2
    directions = Vector{Float64}[]
    for k in 0:3
        angle = absolute_angle + k * π/2
        push!(directions, [cos(angle), sin(angle)])
    end
    
    return directions
end

"""
    compute_vertex_to_faces(topology::TriangularMeshTopology) -> Dict{Int, Vector{Int}}

Build a mapping from vertex indices to the faces that contain them.

# Arguments
- `topology` - Triangular mesh topology

# Returns
- Dictionary mapping vertex index → list of face indices containing that vertex
"""
function compute_vertex_to_faces(topology::TriangularMeshTopology)
    fs = faces(topology.mesh)
    vertex_to_faces = Dict{Int, Vector{Int}}()
    
    for (face_idx, face) in enumerate(fs)
        for v in face
            if !haskey(vertex_to_faces, v)
                vertex_to_faces[v] = Int[]
            end
            push!(vertex_to_faces[v], face_idx)
        end
    end
    
    return vertex_to_faces
end

"""
    compute_angle_defect(topology::TriangularMeshTopology, vertex_idx::Int) -> Float64

Compute the angle defect at a vertex.
For interior vertices: A_d(v) = 2π - sum of corner angles
For boundary vertices: A_d(v) = π - sum of corner angles

# Arguments
- `topology` - Triangular mesh topology
- `vertex_idx` - Index of the vertex

# Returns
- Angle defect in radians
"""
function compute_angle_defect(topology::TriangularMeshTopology, vertex_idx::Int)
    fs = faces(topology.mesh)
    vs = coordinates(topology.mesh)
    
    # Get all faces containing this vertex
    vertex_to_faces = compute_vertex_to_faces(topology)
    incident_faces = get(vertex_to_faces, vertex_idx, Int[])
    
    if isempty(incident_faces)
        return 0.0
    end
    
    # Sum the corner angles at this vertex
    total_angle = 0.0
    for face_idx in incident_faces
        face = fs[face_idx]
        verts = [face[1], face[2], face[3]]
        
        # Find position of vertex_idx in the face
        local_idx = findfirst(==(vertex_idx), verts)
        if isnothing(local_idx)
            continue
        end
        
        # Get the three vertices of the triangle
        v0 = vs[verts[local_idx]]        # The vertex we're computing for
        v1 = vs[verts[mod1(local_idx + 1, 3)]]  # Next vertex
        v2 = vs[verts[mod1(local_idx + 2, 3)]]  # Previous vertex
        
        # Compute vectors from v0 to v1 and v0 to v2
        e1 = [v1[1] - v0[1], v1[2] - v0[2]]
        e2 = [v2[1] - v0[1], v2[2] - v0[2]]
        
        # Compute angle between edges using dot product
        len1 = sqrt(e1[1]^2 + e1[2]^2)
        len2 = sqrt(e2[1]^2 + e2[2]^2)
        if len1 > 1e-10 && len2 > 1e-10
            cos_angle = (e1[1]*e2[1] + e1[2]*e2[2]) / (len1 * len2)
            cos_angle = clamp(cos_angle, -1.0, 1.0)
            total_angle += acos(cos_angle)
        end
    end
    
    # Check if vertex is on boundary (incomplete star)
    # For a 2D planar mesh, interior vertices should have angle sum = 2π
    # We use a simple heuristic: if angle sum is significantly less than 2π, it's a boundary vertex
    is_boundary = total_angle < 1.9 * π
    
    if is_boundary
        return π - total_angle
    else
        return 2π - total_angle
    end
end

"""
    compute_singularities(field::CrossField) -> Vector{Tuple{Int, Float64, Point, Bool}}

Compute the singularity index for all vertices and return non-zero singularities.

For a flat 2D domain, the singularity index simplifies to:
    I(v) = Σ_{e_ij ∈ star(v)} p_ij / 4

Where p_ij are the period jumps on edges in the star of vertex v.

Vertices with non-zero index are singularities:
- Index +1/4 → valence 3 vertex in quad mesh
- Index -1/4 → valence 5 vertex in quad mesh
- Index +1/2 → valence 2 vertex
- Index -1/2 → valence 6 vertex

# Arguments
- `field` - CrossField with solved theta and period jumps

# Returns
- Vector of tuples: (vertex_index, index_value, position, is_boundary)
"""
function compute_singularities(field::CrossField)
    topology = field.topology
    fs = faces(topology.mesh)
    vs = coordinates(topology.mesh)
    n_verts = length(vs)
    
    # Build vertex to faces mapping
    vertex_to_faces = compute_vertex_to_faces(topology)
    
    # Build vertex to edges mapping (edges in star of vertex)
    # An edge is in star(v) if it connects two faces that share v
    vertex_star_edges = Dict{Int, Vector{Tuple{Int, Int}}}()
    
    for v in 1:n_verts
        vertex_star_edges[v] = Tuple{Int, Int}[]
        incident_faces = get(vertex_to_faces, v, Int[])
        
        # Find all pairs of adjacent faces sharing this vertex
        for i in 1:length(incident_faces)
            face_i = incident_faces[i]
            for j in i+1:length(incident_faces)
                face_j = incident_faces[j]
                
                # Check if these faces share an edge (not just the vertex)
                verts_i = Set([fs[face_i][1], fs[face_i][2], fs[face_i][3]])
                verts_j = Set([fs[face_j][1], fs[face_j][2], fs[face_j][3]])
                shared = intersect(verts_i, verts_j)
                
                if length(shared) == 2 && v in shared
                    # These faces share an edge containing v
                    # This edge contributes to the star of v
                    canonical_edge = face_i < face_j ? (face_i, face_j) : (face_j, face_i)
                    if !(canonical_edge in vertex_star_edges[v])
                        push!(vertex_star_edges[v], canonical_edge)
                    end
                end
            end
        end
    end
    
    # Identify boundary vertices (appear in fewer than their face count would suggest)
    # A vertex is on boundary if the total angle around it is less than 2π
    boundary_vertices = Set{Int}()
    for v in 1:n_verts
        incident_faces = get(vertex_to_faces, v, Int[])
        if isempty(incident_faces)
            continue
        end
        
        # Sum angles around vertex
        total_angle = 0.0
        for face_idx in incident_faces
            face = fs[face_idx]
            verts = [face[1], face[2], face[3]]
            local_idx = findfirst(==(v), verts)
            if isnothing(local_idx)
                continue
            end
            v0 = vs[verts[local_idx]]
            v1 = vs[verts[mod1(local_idx + 1, 3)]]
            v2 = vs[verts[mod1(local_idx + 2, 3)]]
            e1 = [v1[1] - v0[1], v1[2] - v0[2]]
            e2 = [v2[1] - v0[1], v2[2] - v0[2]]
            len1 = sqrt(e1[1]^2 + e1[2]^2)
            len2 = sqrt(e2[1]^2 + e2[2]^2)
            if len1 > 1e-10 && len2 > 1e-10
                cos_angle = clamp((e1[1]*e2[1] + e1[2]*e2[2]) / (len1 * len2), -1.0, 1.0)
                total_angle += acos(cos_angle)
            end
        end
        
        # Boundary vertices have total angle < 2π
        if total_angle < 2π - 0.1
            push!(boundary_vertices, v)
        end
    end
    
    # Compute singularity index for each INTERIOR vertex
    singularities = Tuple{Int, Float64, typeof(vs[1]), Bool}[]  # Added boundary flag
    
    for v in 1:n_verts
        is_boundary = v in boundary_vertices
        
        # For interior vertices: angle defect should be 0 on flat domain
        # For boundary vertices: we compute them but mark them separately
        
        # Sum of period jumps over star edges only (for flat 2D domain, kappa sum cancels)
        p_sum = 0
        n_star_edges = length(vertex_star_edges[v])
        
        for (face_i, face_j) in vertex_star_edges[v]
            # Get period jump (oriented)
            p_ij = get_period_jump(field, face_i, face_j)
            p_sum += p_ij
        end
        
        # For a flat 2D domain: I_0 = 0 for interior vertices
        # Singularity index = p_sum / 4
        I_v = p_sum / 4.0
        
        # Only record non-zero singularities (threshold for numerical noise)
        if abs(I_v) > 1e-6
            push!(singularities, (v, I_v, vs[v], is_boundary))
        end
    end
    
    return singularities
end

"""
    singularity_valence(index::Float64) -> Int

Convert singularity index to quad mesh valence.

# Arguments  
- `index` - Singularity index value

# Returns
- Expected valence in quad mesh (4 for regular, 3 for +1/4, 5 for -1/4, etc.)
"""
function singularity_valence(index::Float64)
    # Regular vertex has valence 4
    # Each +1/4 index decreases valence by 1
    # Each -1/4 index increases valence by 1
    return 4 - round(Int, 4 * index)
end
