using GeometryBasics

include("meshIO.jl")

mesh = load_triangulation("SimpleAutoQuad/output/meshes/simple-square.msh")

V = length(coordinates(mesh))
F = length(faces(mesh))
E = sum(length(f) for f in faces(mesh)) / 2

chi = V - E + F

println("Vertices: $V")
println("Edges: $E")
println("Faces: $F")
println("Euler characteristic χ = V - E + F = $chi")
println()
println("For a 4-RoSy field, sum of singularity indices should equal χ = $chi")
