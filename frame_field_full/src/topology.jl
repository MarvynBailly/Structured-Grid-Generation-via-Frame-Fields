module Topology

using ..Types
using Printf

export build_topology, get_edge_vector, compute_transport_angles, compute_boundary_constraints

function build_topology(vertices, faces)
    n = length(faces)
    dual_adj = [Tuple{Int, Tuple{Int,Int}}[] for _ in 1:n]
    face_ref = Vector{Tuple{Int,Int}}(undef, n)
    edge_hist = Dict{Tuple{Int,Int}, Int}()
    for (i, tri) in enumerate(faces)
        face_ref[i] = (tri[1], tri[2])
        vs = [tri[1], tri[2], tri[3]]
        for k in 1:3
            v1, v2 = vs[k], vs[mod1(k+1, 3)]; key = minmax(v1, v2)
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
    p1 = topo.vertices[v1]; p2 = topo.vertices[v2]
    return (p2.x - p1.x, p2.y - p1.y, p2.z - p1.z)
end

function compute_transport_angles(topo)
    kappa = Dict{Tuple{Int,Int}, Float64}()
    println("\n--- Computing Transport Angles (kappa) ---")
    for i in 1:length(topo.faces)
        for (j, _) in topo.dual_adj[i]
            if i >= j continue end
            vi = topo.faces[i]; vj = topo.faces[j]; shared = intersect(vi, vj)
            ev = get_edge_vector(topo, shared[1], shared[2])
            
            ri = topo.face_ref_edges[i]; rvi = get_edge_vector(topo, ri[1], ri[2])
            ai = atan(ev[2], ev[1]) - atan(rvi[2], rvi[1])
            
            rj = topo.face_ref_edges[j]; rvj = get_edge_vector(topo, rj[1], rj[2])
            aj = atan(ev[2], ev[1]) - atan(rvj[2], rvj[1])
            
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
        edges = [(tri[1], tri[2]), (tri[2], tri[3]), (tri[3], tri[1])]
        
        for (k, edge) in enumerate(edges)
            key = minmax(edge[1], edge[2])
            has_neighbor = false
            for (_, n_key) in topo.dual_adj[i]
                if n_key == key; has_neighbor = true; break; end
            end
            
            if !has_neighbor
                third_idx = if k == 1; tri[3]; elseif k == 2; tri[1]; else; tri[2]; end
                p3 = topo.vertices[third_idx]
                
                p1 = topo.vertices[edge[1]]
                p2 = topo.vertices[edge[2]]
                
                dx = p2.x - p1.x
                dy = p2.y - p1.y
                global_tan = atan(dy, dx)
                
                if p3.y < 0.0
                    if dx < 0.0
                        global_tan += π
                    end
                elseif p3.y > 0.0
                    if dx > 0.0
                        global_tan += π
                    end
                end
                
                re = topo.face_ref_edges[i]
                rp1 = topo.vertices[re[1]]
                rp2 = topo.vertices[re[2]]
                ref_ang = atan(rp2.y - rp1.y, rp2.x - rp1.x)
                
                val = mod2pi(global_tan - ref_ang + π) - π
                constraints[i] = val
            end
        end
    end
    return constraints
end

end