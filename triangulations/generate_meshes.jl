"""
Generate simple test triangulations for frame field testing.

Creates well-structured triangular meshes:
1. Regular square with diagonal splits
2. Regular hexagon
3. Disk with radial structure
"""

using GeometryBasics
using FileIO

"""
Generate a regular square domain split into triangles with a grid pattern.
n: number of subdivisions per side
"""
function generate_regular_square(n::Int=4)
    vertices = Point{2, Float64}[]
    faces = TriangleFace{Int}[]
    
    # Generate grid of vertices
    for j in 0:n
        for i in 0:n
            x = Float64(i) / n
            y = Float64(j) / n
            push!(vertices, Point2f(x, y))
        end
    end
    
    # Generate triangle faces
    for j in 0:(n-1)
        for i in 0:(n-1)
            # Bottom-left corner of this quad
            bl = j * (n+1) + i + 1
            br = bl + 1
            tl = bl + (n+1)
            tr = tl + 1
            
            # Split quad into two triangles
            push!(faces, TriangleFace(bl, br, tl))
            push!(faces, TriangleFace(br, tr, tl))
        end
    end
    
    return GeometryBasics.Mesh(vertices, faces)
end

"""
Generate a regular hexagon subdivided into triangles.
"""
function generate_hexagon(n_radial::Int=3)
    vertices = Point{2, Float64}[]
    faces = TriangleFace{Int}[]
    
    # Center vertex
    push!(vertices, Point2f(0.5, 0.5))
    center_idx = 1
    
    # Create concentric rings
    for ring in 1:n_radial
        n_sides = 6 * ring
        radius = 0.4 * ring / n_radial
        
        ring_start = length(vertices) + 1
        
        # Add vertices in this ring
        for i in 0:(n_sides-1)
            angle = 2π * i / n_sides
            x = 0.5 + radius * cos(angle)
            y = 0.5 + radius * sin(angle)
            push!(vertices, Point2f(x, y))
        end
        
        # Connect to previous ring
        if ring == 1
            # Connect first ring to center
            for i in 0:(n_sides-1)
                v1 = ring_start + i
                v2 = ring_start + (i + 1) % n_sides
                push!(faces, TriangleFace(center_idx, v1, v2))
            end
        else
            # Connect to previous ring
            prev_ring_start = ring_start - 6 * (ring - 1)
            prev_n = 6 * (ring - 1)
            
            # Simple connection (not perfect but works)
            for i in 0:(n_sides-1)
                v1 = ring_start + i
                v2 = ring_start + (i + 1) % n_sides
                
                # Find closest vertex in previous ring
                prev_idx = prev_ring_start + (i * prev_n ÷ n_sides) % prev_n
                
                push!(faces, TriangleFace(prev_idx, v1, v2))
            end
        end
    end
    
    return GeometryBasics.Mesh(vertices, faces)
end

"""
Generate a disk with radial triangulation.
"""
function generate_disk(n_radial::Int=4, n_angular::Int=12)
    vertices = Point{2, Float64}[]
    faces = TriangleFace{Int}[]
    
    # Center vertex
    push!(vertices, Point2f(0.5, 0.5))
    center_idx = 1
    
    # Create concentric rings
    for ring in 1:n_radial
        radius = 0.4 * ring / n_radial
        ring_start = length(vertices) + 1
        
        # Add vertices in this ring
        for i in 0:(n_angular-1)
            angle = 2π * i / n_angular
            x = 0.5 + radius * cos(angle)
            y = 0.5 + radius * sin(angle)
            push!(vertices, Point2f(x, y))
        end
        
        if ring == 1
            # Connect first ring to center
            for i in 0:(n_angular-1)
                v1 = ring_start + i
                v2 = ring_start + (i + 1) % n_angular
                push!(faces, TriangleFace(center_idx, v1, v2))
            end
        else
            # Connect to previous ring
            prev_ring_start = ring_start - n_angular
            
            for i in 0:(n_angular-1)
                v1_inner = prev_ring_start + i
                v2_inner = prev_ring_start + (i + 1) % n_angular
                v1_outer = ring_start + i
                v2_outer = ring_start + (i + 1) % n_angular
                
                # Create two triangles for each quad
                push!(faces, TriangleFace(v1_inner, v1_outer, v2_inner))
                push!(faces, TriangleFace(v1_outer, v2_outer, v2_inner))
            end
        end
    end
    
    return GeometryBasics.Mesh(vertices, faces)
end

"""
Write mesh to Gmsh .msh format (version 4.1).
"""
function write_gmsh(filename::String, mesh::GeometryBasics.Mesh)
    vs = coordinates(mesh)
    fs = faces(mesh)
    
    open(filename, "w") do io
        # Header
        println(io, "\$MeshFormat")
        println(io, "4.1 0 8")
        println(io, "\$EndMeshFormat")
        
        # Physical names
        println(io, "\$PhysicalNames")
        println(io, "1")
        println(io, "2 1 \"domain\"")
        println(io, "\$EndPhysicalNames")
        
        # Entities (simplified - just one surface)
        println(io, "\$Entities")
        println(io, "0 0 1 0")
        
        # Find bounding box
        xs = [v[1] for v in vs]
        ys = [v[2] for v in vs]
        minx, maxx = minimum(xs), maximum(xs)
        miny, maxy = minimum(ys), maximum(ys)
        
        println(io, "1 $minx $miny 0 $maxx $maxy 0 1 1 0 0")
        println(io, "\$EndEntities")
        
        # Nodes
        println(io, "\$Nodes")
        n_verts = length(vs)
        println(io, "1 $n_verts 1 $n_verts")
        println(io, "2 1 0 $n_verts")
        for i in 1:n_verts
            println(io, i)
        end
        for v in vs
            # Gmsh uses 3D coordinates
            println(io, "$(v[1]) $(v[2]) 0")
        end
        println(io, "\$EndNodes")
        
        # Elements
        println(io, "\$Elements")
        n_faces = length(fs)
        println(io, "1 $n_faces 1 $n_faces")
        println(io, "2 1 2 $n_faces")
        for (i, f) in enumerate(fs)
            println(io, "$i $(f[1]) $(f[2]) $(f[3])")
        end
        println(io, "\$EndElements")
    end
end

# Generate meshes
println("Generating test triangulations...")

println("  1. Regular square (4×4)")
mesh1 = generate_regular_square(4)
write_gmsh(joinpath(@__DIR__, "regular-square-4x4.msh"), mesh1)

println("  2. Regular square (8×8)")
mesh2 = generate_regular_square(8)
write_gmsh(joinpath(@__DIR__, "regular-square-8x8.msh"), mesh2)

println("  3. Disk (radial)")
mesh3 = generate_disk(4, 12)
write_gmsh(joinpath(@__DIR__, "disk-radial.msh"), mesh3)

println("  4. Disk (fine)")
mesh4 = generate_disk(6, 16)
write_gmsh(joinpath(@__DIR__, "disk-radial-fine.msh"), mesh4)

println("\nGenerated meshes:")
println("  - regular-square-4x4.msh: $(length(coordinates(mesh1))) vertices, $(length(faces(mesh1))) triangles")
println("  - regular-square-8x8.msh: $(length(coordinates(mesh2))) vertices, $(length(faces(mesh2))) triangles")
println("  - disk-radial.msh: $(length(coordinates(mesh3))) vertices, $(length(faces(mesh3))) triangles")
println("  - disk-radial-fine.msh: $(length(coordinates(mesh4))) vertices, $(length(faces(mesh4))) triangles")
println("\nDone!")
