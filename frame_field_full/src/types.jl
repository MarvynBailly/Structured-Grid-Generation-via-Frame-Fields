module Types

export Point3D, MeshTopology, CrossField

struct Point3D
    x::Float64; y::Float64; z::Float64
end

struct MeshTopology
    vertices::Vector{Point3D}
    faces::Vector{Tuple{Int,Int,Int}}
    dual_adj::Vector{Vector{Tuple{Int, Tuple{Int,Int}}}} 
    face_ref_edges::Vector{Tuple{Int,Int}}
end

mutable struct CrossField
    topology::MeshTopology
    theta::Vector{Float64}
    period_jumps::Dict{Tuple{Int,Int}, Int}
    transport_angles::Dict{Tuple{Int,Int}, Float64}
    constrained_faces::Dict{Int, Float64}
    fixed_edges::Set{Tuple{Int,Int}}
end

end