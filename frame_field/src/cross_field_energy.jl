"""
    cross_field_energy.jl

Energy function computation and matrix assembly for frame field optimization.

Implements the smoothness energy:
    E_smooth = Σ (θ_i + κ_ij + (π/2)p_ij - θ_j)²

where the sum is over all interior edges, θ_i are face angles, 
p_ij are period jumps (integers), and κ_ij are transport angles.
"""

using LinearAlgebra
using SparseArrays
using GeometryBasics
using Statistics

"""
    compute_smoothness_energy(
        mesh::Mesh,
        edge_to_faces::Dict{Int, Vector{Int}},
        dual_adj::Dict{Int, Vector{Tuple{Int, Tuple{Int, Int}}}},
        kappa::Dict{Tuple{Int, Int}, Float64},
        theta::Vector{Float64},
        p::Dict{Tuple{Int, Int}, Int}
    ) -> Float64

Compute the frame field smoothness energy.

# Arguments
- `mesh` - Triangular mesh
- `edge_to_faces` - Map from edge index to adjacent face indices
- `dual_adj` - Dual graph adjacency: face → [(neighbor_face, edge)]
- `kappa` - Transport angles κ_ij for each face pair (i,j)
- `theta` - Face angles (length = n_faces)
- `p` - Period jumps indexed by edge (vertex pair tuple)

# Returns
- Energy value E = Σ (θ_i + κ_ij + (π/2)p_ij - θ_j)²
"""
function compute_smoothness_energy(
    mesh::Mesh,
    edge_to_faces::Dict{Int, Vector{Int}},
    dual_adj::Dict{Int, Vector{Tuple{Int, Tuple{Int, Int}}}},
    kappa::Dict{Tuple{Int, Int}, Float64},
    theta::Vector{Float64},
    p::Dict{Tuple{Int, Int}, Int}
)
    energy = 0.0
    fs = faces(mesh)
    n_faces = length(fs)
    
    # Sum over all edges (count each once)
    for face_i in 1:n_faces
        for (face_j, edge) in dual_adj[face_i]
            if face_i < face_j  # Count each edge once
                # Kappa stored for edge (use + for i→j per edge_to_faces)
                kappa_ij = get(kappa, edge, 0.0)
                p_ij = get(p, edge, 0)
                
                diff = theta[face_i] + kappa_ij + (π/2) * p_ij - theta[face_j]
                energy += diff^2
            end
        end
    end
    
    return energy
end

"""
    compute_transport_angle(p1::Point, p2::Point, p3::Point, 
                           q1::Point, q2::Point, q3::Point) -> Float64

Compute the transport angle κ_ij between two adjacent triangles.

The transport angle is the rotation needed to align the local coordinate
frames of the two triangles across their shared edge.

# Arguments
- `p1, p2, p3` - Vertices of first triangle
- `q1, q2, q3` - Vertices of second triangle (q1, q2 is shared edge)

# Returns
- Transport angle κ_ij in radians
"""
function compute_transport_angle(
    p1, p2, p3,
    q1, q2, q3
)
    # Compute normal for first triangle
    v1 = p2 - p1
    v2 = p3 - p1
    n1 = normalize(cross(v1, v2))
    
    # Compute normal for second triangle
    u1 = q2 - q1
    u2 = q3 - q1
    n2 = normalize(cross(u1, u2))
    
    # Find shared edge direction
    shared_edge = p2 - p1  # Assuming p1, p2 are on shared edge
    e = normalize(shared_edge)
    
    # Project normals onto plane perpendicular to edge
    n1_proj = normalize(n1 - dot(n1, e) * e)
    n2_proj = normalize(n2 - dot(n2, e) * e)
    
    # Compute angle between projected normals
    cos_angle = clamp(dot(n1_proj, n2_proj), -1.0, 1.0)
    sin_angle = dot(cross(n1_proj, n2_proj), e)
    
    return atan(sin_angle, cos_angle)
end

"""
    assemble_system_matrix(
        mesh::Mesh,
        potential_fixed_edges::Set{Tuple{Int, Int}},
        fixed_edges_per_face::Dict{Int, Tuple{Int, Int}},
        constrained_faces::Vector{Int}=Int[],
        constrained_angles::Vector{Float64}=Float64[]
    ) -> (SparseMatrixCSC, Vector{Float64}, Dict, Dict)

Assemble the linear system for frame field optimization.

Creates equations for all variables (both free and constrained):
- For free edges: ∂E/∂p_ij = π(θ_i + κ_ij + (π/2)p_ij - θ_j) = 0
- For fixed edges: p_ij = 0
- For free faces: ∂E/∂θ_k = Σ_j 2(θ_k + κ_kj + (π/2)p_kj - θ_j) = 0
- For constrained faces: θ_k = fixed_value

# Arguments
- `mesh` - Triangular mesh
- `potential_fixed_edges` - Edges in spanning forest (from compute_spanning_forest)
- `fixed_edges_per_face` - Dict mapping face_idx → edge to fix (from fixed_suitable_edges)
- `constrained_faces` - Face indices with fixed angles
- `constrained_angles` - Corresponding angle values

# Returns
- `A::SparseMatrixCSC` - Sparse system matrix (n_vars × n_vars)
- `b::Vector{Float64}` - Right-hand side vector
- `var_to_p::Dict{Int, Tuple{Int,Int}}` - Maps variable index → edge (for period jumps)
- `var_to_theta::Dict{Int, Int}` - Maps variable index → face index (for angles)

# Variable Ordering
Variables are ordered as: [p_free; p_fixed; θ_free; θ_constrained]
- First n_free_p variables: period jumps on free edges
- Next n_fixed_p variables: period jumps on fixed edges (constrained to 0)
- Next n_free_theta variables: angles on free faces
- Last n_constrained_theta variables: angles on constrained faces
"""
function assemble_system_matrix(
    mesh::Mesh,
    potential_fixed_edges::Set{Tuple{Int, Int}},
    fixed_edges_per_face::Dict{Int, Tuple{Int, Int}},
    constrained_faces::Vector{Int}=Int[],
    constrained_angles::Vector{Float64}=Float64[],
    debug::Bool=false
)
    fs = faces(mesh)
    vs = coordinates(mesh)
    n_faces = length(fs)
    
    # Build topology
    edge_map = Dict{Tuple{Int, Int}, Int}()
    edge_to_faces = Dict{Int, Vector{Int}}()
    dual_adj = Dict{Int, Vector{Tuple{Int, Tuple{Int, Int}}}}()
    
    edge_idx = 0
    for i in 1:n_faces
        dual_adj[i] = []
    end
    
    # Build edge map and face adjacency
    for (i, face_i) in enumerate(fs)
        verts_i = [face_i[1], face_i[2], face_i[3]]
        
        for (j, face_j) in enumerate(fs)
            if i >= j
                continue
            end
            
            verts_j = [face_j[1], face_j[2], face_j[3]]
            shared = intersect(verts_i, verts_j)
            
            if length(shared) == 2
                v1, v2 = shared[1], shared[2]
                edge = v1 < v2 ? (v1, v2) : (v2, v1)
                
                if !haskey(edge_map, edge)
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
    
    # Compute transport angles κ_ij
    kappa = Dict{Tuple{Int, Int}, Float64}()
    
    for face_i in 1:n_faces
        for (face_j, edge) in dual_adj[face_i]
            if face_i < face_j
                # Get vertices for both triangles
                f_i = fs[face_i]
                f_j = fs[face_j]
                
                pi1, pi2, pi3 = vs[f_i[1]], vs[f_i[2]], vs[f_i[3]]
                pj1, pj2, pj3 = vs[f_j[1]], vs[f_j[2]], vs[f_j[3]]
                
                # Compute local frames for both faces
                # X-axis is normalized edge from p1 to p2
                ei1 = pi2 - pi1
                xi = normalize(ei1)
                ni = normalize(cross(ei1, pi3 - pi1))
                yi = cross(ni, xi)
                
                ej1 = pj2 - pj1
                xj = normalize(ej1)
                nj = normalize(cross(ej1, pj3 - pj1))
                yj = cross(nj, xj)
                
                # Find shared edge vertices
                common = intersect([f_i[1], f_i[2], f_i[3]], [f_j[1], f_j[2], f_j[3]])
                if length(common) == 2
                    p_u = vs[common[1]]
                    p_v = vs[common[2]]
                    edge_vec = normalize(p_v - p_u)
                    
                    # Angle of edge in frame i
                    ang_i = atan(dot(edge_vec, yi), dot(edge_vec, xi))
                    
                    # Angle of edge in frame j  
                    ang_j = atan(dot(edge_vec, yj), dot(edge_vec, xj))
                    
                    # Transport angle
                    kappa_ij = ang_j - ang_i
                    
                    # Wrap to (-π, π]
                    kappa_ij = mod2pi(kappa_ij + π) - π
                else
                    kappa_ij = 0.0
                end
                
                kappa[(face_i, face_j)] = kappa_ij
                kappa[(face_j, face_i)] = -kappa_ij
            end
        end
    end
    
    # Debug kappa values
    nonzero_kappas = filter(p -> abs(p.second) > 1e-10, kappa)
    println("  Kappa values computed: $(length(kappa))")
    println("  Non-zero kappas: $(length(nonzero_kappas))")
    if length(nonzero_kappas) > 0
        kappa_vals = [abs(v) for v in values(kappa)]
        println("  Kappa range: [$(minimum(kappa_vals)), $(maximum(kappa_vals))]")
        println("  First 5 non-zero kappas:")
        for (i, ((fi, fj), k)) in enumerate(collect(nonzero_kappas)[1:min(5, length(nonzero_kappas))])
            println("    κ[($fi, $fj)] = $k")
        end
    end
    
    # Determine free and fixed variables
    fixed_edges = Set(values(fixed_edges_per_face))
    free_faces = setdiff(1:n_faces, constrained_faces)
    
    # All edges are variables (both free and fixed)
    all_edges = collect(keys(edge_map))
    free_edges = [e for e in all_edges if !(e in fixed_edges)]
    
    n_all_p = length(all_edges)
    n_free_p = length(free_edges)
    n_fixed_p = n_all_p - n_free_p
    n_all_theta = n_faces
    n_free_theta = length(free_faces)
    n_fixed_theta = n_all_theta - n_free_theta
    n_vars = n_all_p + n_all_theta

    debug && println("Total number of variables: $n_vars")
    debug && println("- Period jump variables: $n_all_p (free: $n_free_p, fixed: $n_fixed_p)")
    debug && println("- Angle variables: $n_all_theta (free: $n_free_theta, fixed: $n_fixed_theta)")
    
    # Create variable mappings
    # Period jumps: variables 1 to n_all_p (free edges first, then fixed edges)
    var_to_p = Dict{Int, Tuple{Int, Int}}()
    p_to_var = Dict{Tuple{Int, Int}, Int}()
    p_is_fixed = Dict{Int, Bool}()
    
    var_idx = 0
    # First add free edges
    for edge in free_edges
        var_idx += 1
        var_to_p[var_idx] = edge
        p_to_var[edge] = var_idx
        p_is_fixed[var_idx] = false
    end
    # Then add fixed edges
    for edge in fixed_edges
        var_idx += 1
        var_to_p[var_idx] = edge
        p_to_var[edge] = var_idx
        p_is_fixed[var_idx] = true
    end
    
    # Angles: variables (n_all_p + 1) to n_vars (all faces: free first, then constrained)
    var_to_theta = Dict{Int, Int}()
    theta_to_var = Dict{Int, Int}()
    theta_is_fixed = Dict{Int, Bool}()
    
    var_idx = n_all_p
    # First add free faces
    for face_idx in sort(collect(free_faces))
        var_idx += 1
        var_to_theta[var_idx] = face_idx
        theta_to_var[face_idx] = var_idx
        theta_is_fixed[var_idx] = false
    end
    # Then add constrained faces
    for face_idx in sort(constrained_faces)
        var_idx += 1
        var_to_theta[var_idx] = face_idx
        theta_to_var[face_idx] = var_idx
        theta_is_fixed[var_idx] = true
    end
    
    # Build sparse system using coordinate format
    I_idx = Int[]
    J_idx = Int[]
    vals = Float64[]
    b = zeros(n_vars)
    
    row = 0
    
    # Equations from ∂E/∂p_ij = 0 for each period jump variable
    for (var_idx, edge) in sort(collect(var_to_p))
        row += 1
        
        # Check if this is a fixed edge
        if p_is_fixed[var_idx]
            # Fixed edge constraint: p_ij = 0
            push!(I_idx, row)
            push!(J_idx, var_idx)
            push!(vals, 1.0)
            b[row] = 0.0
            continue
        end
        
        # Free edge: use energy derivative equation
        # Find faces incident to this edge
        edge_idx = edge_map[edge]
        incident_faces = edge_to_faces[edge_idx]
        
        if length(incident_faces) != 2
            # Boundary edge - add constraint p_ij = 0
            push!(I_idx, row)
            push!(J_idx, var_idx)
            push!(vals, 1.0)
            b[row] = 0.0
            continue
        end
        
        face_i, face_j = incident_faces
        # Kappa stored for edge in i→j direction
        kappa_ij = get(kappa, edge, 0.0)
        
        # ∂E/∂p_ij = π(θ_i + κ_ij + (π/2)p_ij - θ_j) = 0
        # Simplifies to: (π²/2)p_ij + πθ_i - πθ_j = -πκ_ij
        
        # Coefficient for p_ij itself: π² / 2
        push!(I_idx, row)
        push!(J_idx, var_idx)
        push!(vals, π^2 / 2)
        
        # Coefficient for θ_i: π (always include as variable)
        theta_var_i = theta_to_var[face_i]
        push!(I_idx, row)
        push!(J_idx, theta_var_i)
        push!(vals, π)
        
        # Coefficient for θ_j: -π (always include as variable)
        theta_var_j = theta_to_var[face_j]
        push!(I_idx, row)
        push!(J_idx, theta_var_j)
        push!(vals, -π)
        
        # Constant term from κ_ij
        b[row] -= π * kappa_ij
    end
    
    # Equations from ∂E/∂θ_k = 0 for each theta variable
    for (var_idx, face_k) in sort(collect(var_to_theta))
        row += 1
        
        # Check if this is a constrained face
        if theta_is_fixed[var_idx]
            # Fixed theta constraint: θ_k = constrained_value
            push!(I_idx, row)
            push!(J_idx, var_idx)
            push!(vals, 1.0)
            idx = findfirst(==(face_k), constrained_faces)
            b[row] = constrained_angles[idx]
            continue
        end
        
        # Free face: use energy derivative equation
        # Sum over all neighbors of face_k
        degree = length(dual_adj[face_k])
        
        # Diagonal term: 2 * degree
        push!(I_idx, row)
        push!(J_idx, var_idx)
        push!(vals, 2.0 * degree)
        
        for (face_j, edge) in dual_adj[face_k]
            kappa_kj = kappa[(face_k, face_j)]
            
            # Coefficient for p_kj: π (always include as variable)
            p_var = p_to_var[edge]
            push!(I_idx, row)
            push!(J_idx, p_var)
            push!(vals, π)
            
            # Coefficient for θ_j: -2 (always include as variable)
            theta_var_j = theta_to_var[face_j]
            push!(I_idx, row)
            push!(J_idx, theta_var_j)
            push!(vals, -2.0)
            
            # Constant term from κ_kj
            b[row] -= 2.0 * kappa_kj
        end
    end
    
    # Build sparse matrix
    A = sparse(I_idx, J_idx, vals, n_vars, n_vars)
    
    return A, b, var_to_p, var_to_theta
end

"""
    assemble_free_system(
        mesh::Mesh,
        potential_fixed_edges::Set{Tuple{Int, Int}},
        fixed_edges_per_face::Dict{Int, Tuple{Int, Int}},
        constrained_faces::Vector{Int}=Int[],
        constrained_angles::Vector{Float64}=Float64[],
        debug::Bool=false
    ) -> (SparseMatrixCSC, Vector{Float64}, Dict, Dict)

Assemble the linear system for ONLY the free variables (for use with greedy solver).

This implementation follows the reference MIQ approach:
- Accumulates energy derivatives edge-by-edge
- Each edge contributes to multiple equation rows
- Creates a properly conditioned SPD system

Fixed edges are eliminated by setting p_ij = 0 and substituting into equations.
Constrained faces are eliminated by substituting their fixed angles.

Returns a system where:
- First n_free_p variables are FREE period jumps (integers)
- Last n_free_theta variables are FREE face angles (reals)

This is the correct format for the greedy mixed-integer solver.
"""
function assemble_free_system(
    mesh::GeometryBasics.Mesh,
    potential_fixed_edges::Set{Tuple{Int, Int}},
    fixed_edges_per_face::Dict{Int, Tuple{Int, Int}},
    constrained_faces::Vector{Int}=Int[],
    constrained_angles::Vector{Float64}=Float64[],
    debug::Bool=false
)
    fs = faces(mesh)
    vs = coordinates(mesh)
    n_faces = length(fs)
    
    # Build topology
    edge_map = Dict{Tuple{Int, Int}, Int}()
    edge_to_faces = Dict{Int, Vector{Int}}()
    dual_adj = Dict{Int, Vector{Tuple{Int, Tuple{Int, Int}}}}()
    
    edge_idx = 0
    for i in 1:n_faces
        dual_adj[i] = []
    end
    
    # Build edge map and face adjacency
    for (i, face_i) in enumerate(fs)
        verts_i = [face_i[1], face_i[2], face_i[3]]
        
        for (j, face_j) in enumerate(fs)
            if i >= j
                continue
            end
            
            verts_j = [face_j[1], face_j[2], face_j[3]]
            shared = intersect(verts_i, verts_j)
            
            if length(shared) == 2
                v1, v2 = shared[1], shared[2]
                edge = v1 < v2 ? (v1, v2) : (v2, v1)
                
                if !haskey(edge_map, edge)
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
    
    # Compute transport angles κ indexed by edge (matches reference)
    kappa = Dict{Tuple{Int, Int}, Float64}()  # edge -> kappa value
    
    for (edge, edge_idx) in edge_map
        faces_list = edge_to_faces[edge_idx]
        if length(faces_list) == 2
            i, j = faces_list
            
            # Get vertices for Face i
            f_i = fs[i]
            pi1, pi2, pi3 = vs[f_i[1]], vs[f_i[2]], vs[f_i[3]]
            ei1 = pi2 - pi1
            xi = normalize(ei1)
            ni = normalize(cross(ei1, pi3 - pi1))
            yi = cross(ni, xi)
            
            # Get vertices for Face j
            f_j = fs[j]
            pj1, pj2, pj3 = vs[f_j[1]], vs[f_j[2]], vs[f_j[3]]
            ej1 = pj2 - pj1
            xj = normalize(ej1)
            nj = normalize(cross(ej1, pj3 - pj1))
            yj = cross(nj, xj)
            
            # Find shared edge vertices
            common = intersect([f_i[1], f_i[2], f_i[3]], [f_j[1], f_j[2], f_j[3]])
            if length(common) == 2
                p_u = vs[common[1]]
                p_v = vs[common[2]]
                edge_vec = normalize(p_v - p_u)
                
                # Angle of edge in frame i
                ang_i = atan(dot(edge_vec, yi), dot(edge_vec, xi))
                
                # Angle of edge in frame j  
                ang_j = atan(dot(edge_vec, yj), dot(edge_vec, xj))
                
                # Transport angle (stored for edge in i→j direction)
                k_val = ang_j - ang_i
                
                # Wrap to (-π, π]
                k_val = mod2pi(k_val + π) - π
                
                kappa[edge] = k_val
            else
                kappa[edge] = 0.0
            end
        end
    end
    
    # Debug kappa values
    if debug
        nonzero_kappas = filter(p -> abs(p.second) > 1e-10, kappa)
        println("  Kappa values computed: $(length(kappa))")
        println("  Non-zero kappas: $(length(nonzero_kappas))")
        if length(nonzero_kappas) > 0
            kappa_vals = [abs(v) for v in values(kappa)]
            println("  Kappa range: [$(minimum(kappa_vals)), $(maximum(kappa_vals))]")
            println("  First 5 non-zero kappas:")
            for (i, (edge, k)) in enumerate(collect(nonzero_kappas)[1:min(5, length(nonzero_kappas))])
                faces = edge_to_faces[edge_map[edge]]
                println("    κ[edge $edge, faces $(faces)] = $k")
            end
        end
    end
    
    # Determine free and fixed variables
    # Use the potential_fixed_edges (spanning forest) as the actual fixed edges
    # This matches the reference implementation which fixes all spanning forest edges
    fixed_edges = potential_fixed_edges
    all_edges = collect(keys(edge_map))
    free_edges = [e for e in all_edges if !(e in fixed_edges)]
    free_faces = setdiff(1:n_faces, constrained_faces)
    
    n_free_p = length(free_edges)
    n_free_theta = length(free_faces)
    n_vars = n_free_p + n_free_theta
    
    debug && println("Free system dimensions:")
    debug && println("  Free period jumps: $n_free_p")
    debug && println("  Free angles: $n_free_theta")
    debug && println("  Total free variables: $n_vars")
    
    # Create variable mappings (ONLY free variables)
    var_to_p = Dict{Int, Tuple{Int, Int}}()
    p_to_var = Dict{Tuple{Int, Int}, Int}()
    
    for (i, edge) in enumerate(free_edges)
        var_to_p[i] = edge
        p_to_var[edge] = i
    end
    
    var_to_theta = Dict{Int, Int}()
    theta_to_var = Dict{Int, Int}()
    
    for (i, face_idx) in enumerate(sort(collect(free_faces)))
        var_idx = n_free_p + i
        var_to_theta[var_idx] = face_idx
        theta_to_var[face_idx] = var_idx
    end
    
    # Build sparse system using edge-by-edge accumulation (like reference implementation)
    # This creates an SPD system by accumulating derivative contributions
    I_idx = Int[]
    J_idx = Int[]
    vals = Float64[]
    b = zeros(n_vars)
    
    # Helper to add matrix entries
    function add_entry!(row, col, val)
        push!(I_idx, row)
        push!(J_idx, col)
        push!(vals, val)
    end
    
    # Iterate over all edges and accumulate contributions
    # This follows the reference implementation pattern
    for (e_idx, edge) in enumerate(all_edges)
        incident_faces = edge_to_faces[edge_map[edge]]
        
        if length(incident_faces) != 2
            continue  # Skip boundary edges
        end
        
        face_i, face_j = incident_faces
        # Kappa is stored for edge in i→j direction (per edge_to_faces order)
        kappa_ij = get(kappa, edge, 0.0)
        
        # Determine which variables are free and get their indices
        row_i = haskey(theta_to_var, face_i) ? theta_to_var[face_i] : 0
        row_j = haskey(theta_to_var, face_j) ? theta_to_var[face_j] : 0
        row_p = haskey(p_to_var, edge) ? p_to_var[edge] : 0
        
        # Get fixed values for constrained variables
        theta_i_val = 0.0
        theta_j_val = 0.0
        p_val = 0.0
        
        if row_i == 0  # face_i is constrained
            idx = findfirst(==(face_i), constrained_faces)
            if !isnothing(idx)
                theta_i_val = constrained_angles[idx]
            end
        end
        
        if row_j == 0  # face_j is constrained
            idx = findfirst(==(face_j), constrained_faces)
            if !isnothing(idx)
                theta_j_val = constrained_angles[idx]
            end
        end
        
        if row_p == 0  # edge is fixed
            p_val = 0.0  # Fixed edges have p = 0
        end
        
        # Compute RHS constant term
        const_RHS = -kappa_ij
        if row_i == 0
            const_RHS -= theta_i_val
        end
        if row_j == 0
            const_RHS += theta_j_val
        end
        if row_p == 0
            const_RHS -= (π/2) * p_val
        end
        
        # Add contributions to equation for ∂E/∂θ_i (if face_i is free)
        if row_i > 0
            # Diagonal term
            add_entry!(row_i, row_i, 2.0)
            
            # Coupling with θ_j
            if row_j > 0
                add_entry!(row_i, row_j, -2.0)
            end
            
            # Coupling with p_ij
            if row_p > 0
                add_entry!(row_i, row_p, π)
            end
            
            # RHS contribution
            b[row_i] += 2.0 * const_RHS
        end
        
        # Add contributions to equation for ∂E/∂θ_j (if face_j is free)
        if row_j > 0
            # Coupling with θ_i
            if row_i > 0
                add_entry!(row_j, row_i, -2.0)
            end
            
            # Diagonal term
            add_entry!(row_j, row_j, 2.0)
            
            # Coupling with p_ij
            if row_p > 0
                add_entry!(row_j, row_p, -π)
            end
            
            # RHS contribution
            b[row_j] += -2.0 * const_RHS
        end
        
        # Add contributions to equation for ∂E/∂p_ij (if edge is free)
        if row_p > 0
            # Coupling with θ_i
            if row_i > 0
                add_entry!(row_p, row_i, π)
            end
            
            # Coupling with θ_j
            if row_j > 0
                add_entry!(row_p, row_j, -π)
            end
            
            # Diagonal term
            add_entry!(row_p, row_p, π^2 / 2.0)
            
            # RHS contribution
            b[row_p] += π * const_RHS
        end
    end
    
    # Build sparse matrix and make it symmetric
    A = sparse(I_idx, J_idx, vals, n_vars, n_vars)
    
    # The accumulation creates a symmetric matrix, ensure it's exactly symmetric
    A = (A + A') / 2.0
    
    return A, b, var_to_p, var_to_theta
end 