module Topology

using ..Types
using Printf
using LinearAlgebra # Added for cross/dot products

export build_topology, get_edge_vector, compute_transport_angles, compute_boundary_constraints

function build_topology(vertices, faces)
    n = length(faces)
    dual_adj = [Tuple{Int, Tuple{Int,Int}}[] for _ in 1:n]
    face_ref = Vector{Tuple{Int,Int}}(undef, n)
    edge_hist = Dict{Tuple{Int,Int}, Int}()
    for (i, tri) in enumerate(faces)
        face_ref[i] = (tri[1], tri[2]) # Reference edge is the first edge of the triangle
        vs = [tri[1], tri[2], tri[3]]
        for k in 1:3
            v1, v2 = vs[k], vs[mod1(k+1, 3)]
            key = minmax(v1, v2)
            if haskey(edge_hist, key)
                n_neighbor = edge_hist[key]
                push!(dual_adj[i], (n_neighbor, key))
                push!(dual_adj[n_neighbor], (i, key))
            else
                edge_hist[key] = i
            end
        end
    end
    return MeshTopology(vertices, faces, dual_adj, face_ref)
end

function get_edge_vector(topo, v1, v2)
    p1 = topo.vertices[v1]
    p2 = topo.vertices[v2]
    return (p2.x - p1.x, p2.y - p1.y, p2.z - p1.z)
end

# Helper to compute angle of vector 'v' relative to basis 'basis_x' in the plane with normal 'n'
function angle_in_plane(v_vec, basis_x, normal)
    # Project v onto basis vectors
    # x_comp = dot(v, basis_x)
    # y_vec = cross(normal, basis_x)
    # y_comp = dot(v, y_vec)
    
    x_comp = v_vec[1]*basis_x[1] + v_vec[2]*basis_x[2] + v_vec[3]*basis_x[3]
    
    # Cross product: normal x basis_x
    y_basis_1 = normal[2]*basis_x[3] - normal[3]*basis_x[2]
    y_basis_2 = normal[3]*basis_x[1] - normal[1]*basis_x[3]
    y_basis_3 = normal[1]*basis_x[2] - normal[2]*basis_x[1]
    
    y_comp = v_vec[1]*y_basis_1 + v_vec[2]*y_basis_2 + v_vec[3]*y_basis_3
    
    return atan(y_comp, x_comp)
end

function compute_transport_angles(topo)
    kappa = Dict{Tuple{Int,Int}, Float64}()
    println("\n--- Computing Transport Angles (kappa) ---")
    
    for i in 1:length(topo.faces)
        # 1. Construct Local Frame for Face i
        tri_i = topo.faces[i]
        
        # Get coordinates
        p1 = topo.vertices[tri_i[1]]
        p2 = topo.vertices[tri_i[2]]
        p3 = topo.vertices[tri_i[3]]
        
        # Edges
        u = (p2.x-p1.x, p2.y-p1.y, p2.z-p1.z)
        v = (p3.x-p1.x, p3.y-p1.y, p3.z-p1.z)
        
        # Normal
        nx = u[2]*v[3] - u[3]*v[2]
        ny = u[3]*v[1] - u[1]*v[3]
        nz = u[1]*v[2] - u[2]*v[1]
        n_len = sqrt(nx^2 + ny^2 + nz^2)
        normal_i = (nx/n_len, ny/n_len, nz/n_len)
        
        # Basis X (Reference Edge)
        re_i = topo.face_ref_edges[i]
        ref_vec_i = get_edge_vector(topo, re_i[1], re_i[2])
        rv_len_i = sqrt(ref_vec_i[1]^2 + ref_vec_i[2]^2 + ref_vec_i[3]^2)
        basis_x_i = (ref_vec_i[1]/rv_len_i, ref_vec_i[2]/rv_len_i, ref_vec_i[3]/rv_len_i)

        for (j, _) in topo.dual_adj[i]
            if i >= j continue end
            
            # 2. Construct Local Frame for Face j
            tri_j = topo.faces[j]
            jp1, jp2, jp3 = topo.vertices[tri_j[1]], topo.vertices[tri_j[2]], topo.vertices[tri_j[3]]
            ju = (jp2.x-jp1.x, jp2.y-jp1.y, jp2.z-jp1.z)
            jv = (jp3.x-jp1.x, jp3.y-jp1.y, jp3.z-jp1.z)
            jnx = ju[2]*jv[3] - ju[3]*jv[2]
            jny = ju[3]*jv[1] - ju[1]*jv[3]
            jnz = ju[1]*jv[2] - ju[2]*jv[1]
            jn_len = sqrt(jnx^2 + jny^2 + jnz^2)
            normal_j = (jnx/jn_len, jny/jn_len, jnz/jn_len)
            
            re_j = topo.face_ref_edges[j]
            ref_vec_j = get_edge_vector(topo, re_j[1], re_j[2])
            rv_len_j = sqrt(ref_vec_j[1]^2 + ref_vec_j[2]^2 + ref_vec_j[3]^2)
            basis_x_j = (ref_vec_j[1]/rv_len_j, ref_vec_j[2]/rv_len_j, ref_vec_j[3]/rv_len_j)
            
            # 3. Compute Angle of Shared Edge in both frames
            vi = topo.faces[i]
            vj = topo.faces[j]
            shared = intersect(vi, vj)
            
            # Edge vector from shared[1] to shared[2]
            ev = get_edge_vector(topo, shared[1], shared[2])
            
            # Angle in Face i
            ai = angle_in_plane(ev, basis_x_i, normal_i)
            
            # Angle in Face j
            aj = angle_in_plane(ev, basis_x_j, normal_j)
            
            # 4. Compute Transport Angle kappa_ij
            kappa[(i, j)] = mod2pi(aj - ai + π) - π
        end
    end
    println("Computed $(length(kappa)) transport angles.")
    return kappa
end

function compute_boundary_constraints(topo)
    constraints = Dict{Int, Float64}()
    
    for i in 1:length(topo.faces)
        tri = topo.faces[i]
        
        # Local Frame
        p1 = topo.vertices[tri[1]]
        p2 = topo.vertices[tri[2]]
        p3 = topo.vertices[tri[3]]
        
        u = (p2.x-p1.x, p2.y-p1.y, p2.z-p1.z)
        v = (p3.x-p1.x, p3.y-p1.y, p3.z-p1.z)
        
        nx = u[2]*v[3] - u[3]*v[2]
        ny = u[3]*v[1] - u[1]*v[3]
        nz = u[1]*v[2] - u[2]*v[1]
        n_len = sqrt(nx^2 + ny^2 + nz^2)
        normal = (nx/n_len, ny/n_len, nz/n_len)
        
        re = topo.face_ref_edges[i]
        rv = get_edge_vector(topo, re[1], re[2])
        rv_len = sqrt(rv[1]^2 + rv[2]^2 + rv[3]^2)
        basis_x = (rv[1]/rv_len, rv[2]/rv_len, rv[3]/rv_len)
        
        edges = [(tri[1], tri[2]), (tri[2], tri[3]), (tri[3], tri[1])]
        
        for (k, edge) in enumerate(edges)
            key = minmax(edge[1], edge[2])
            has_neighbor = false
            for (_, n_key) in topo.dual_adj[i]
                if n_key == key; has_neighbor = true; break; end
            end
            
            if !has_neighbor
                # Boundary Edge
                ev = get_edge_vector(topo, edge[1], edge[2])
                
                # Angle of this boundary edge in local frame
                edge_ang = angle_in_plane(ev, basis_x, normal)
                
                constraints[i] = edge_ang
            end
        end
    end
    return constraints
end

end