module Extraction

using ..Types
using ..Topology
using ..Parameterization
using LinearAlgebra
using DataStructures
using Printf

export QuadMesh, extract_quad_mesh

# --- Data Structures ---

struct QuadMesh
    vertices::Vector{Point3D}
    quads::Vector{Tuple{Int,Int,Int,Int}}
end

# --- Helpers ---

"""
    get_tri_uvs(topo, sys, x, f_idx)

Returns the (u, v) coordinates of the 3 corners of face f_idx.
"""
function get_tri_uvs(topo, sys, x, f_idx)
    ids = sys.v_map[f_idx, :]
    N = sys.n_continuous_vars
    
    us = [x[ids[1]], x[ids[2]], x[ids[3]]]
    vs = [x[ids[1]+N], x[ids[2]+N], x[ids[3]+N]]
    
    return us, vs
end

"""
    solve_barycentric(p, a, b, c)

Solves for barycentric coords of p in triangle abc (2D).
"""
function solve_barycentric(p, a, b, c)
    v0 = b - a
    v1 = c - a
    v2 = p - a
    denom = v0[1]*v1[2] - v1[1]*v0[2]
    if abs(denom) < 1e-12; return [1/3, 1/3, 1/3]; end
    
    v = (v2[1]*v1[2] - v1[1]*v2[2]) / denom
    w = (v0[1]*v2[2] - v2[1]*v0[2]) / denom
    u = 1.0 - v - w
    return [u, v, w]
end

# --- Simple Grid Sampling ---

function extract_quad_mesh(topo::MeshTopology, 
                           u_phys::Vector{Float64}, 
                           v_phys::Vector{Float64}, 
                           x::Vector{Float64}, 
                           sys, 
                           rotations::Vector{Int})
    
    println("\n--- Extracting Quad Mesh (Simple Grid Sampling) ---")
    
    # Map (integer_u, integer_v) -> Vertex Index
    node_map = Dict{Tuple{Int,Int}, Int}()
    quad_verts = Point3D[]
    
    # Helper to get/create vertex index
    function get_node_idx(ui, vi, p3d)
        key = (ui, vi)
        if haskey(node_map, key)
            return node_map[key]
        else
            push!(quad_verts, p3d)
            idx = length(quad_verts)
            node_map[key] = idx
            return idx
        end
    end
    
    # Iterate all triangles and sample integer grid points
    println("  Sampling integer grid...")
    for f in 1:length(topo.faces)
        us, vs = get_tri_uvs(topo, sys, x, f)
        
        min_u, max_u = floor(Int, minimum(us)), ceil(Int, maximum(us))
        min_v, max_v = floor(Int, minimum(vs)), ceil(Int, maximum(vs))
        
        eps = 1e-5
        
        for i in min_u:max_u
            for j in min_v:max_v
                # Check if (i, j) is inside triangle
                bary = solve_barycentric([Float64(i), Float64(j)], 
                                         [us[1], vs[1]], [us[2], vs[2]], [us[3], vs[3]])
                
                if all(bary .>= -eps)
                    # Compute 3D position
                    p1, p2, p3 = topo.vertices[topo.faces[f][1]], 
                                 topo.vertices[topo.faces[f][2]], 
                                 topo.vertices[topo.faces[f][3]]
                    px = bary[1]*p1.x + bary[2]*p2.x + bary[3]*p3.x
                    py = bary[1]*p1.y + bary[2]*p2.y + bary[3]*p3.y
                    pz = bary[1]*p1.z + bary[2]*p2.z + bary[3]*p3.z
                    
                    get_node_idx(i, j, Point3D(px, py, pz))
                end
            end
        end
    end
    
    println("  Found $(length(quad_verts)) unique grid vertices.")
    
    # Generate Quads by connecting neighboring grid points
    quads = Tuple{Int,Int,Int,Int}[]
    
    for (key, idx_bl) in node_map
        i, j = key
        
        # Check for neighbors to form quad (BL -> BR -> TR -> TL)
        if haskey(node_map, (i+1, j)) && 
           haskey(node_map, (i+1, j+1)) && 
           haskey(node_map, (i, j+1))
           
            idx_br = node_map[(i+1, j)]
            idx_tr = node_map[(i+1, j+1)]
            idx_tl = node_map[(i, j+1)]
            
            # Form Quad (CCW)
            push!(quads, (idx_bl, idx_br, idx_tr, idx_tl))
        end
    end
    
    println("  Extracted $(length(quads)) quads.")
    return QuadMesh(quad_verts, quads)
end

end