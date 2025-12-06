using FileIO
using GeometryBasics

"""
    load_triangulation(filepath::String)

Loads a mesh and standardizes it into a GeometryBasics Mesh 
containing only 1-based integer TriangleFaces.
"""
function load_triangulation(filepath::String)
    if !isfile(filepath)
        error("File not found: $filepath")
    end

    raw_mesh = load(filepath)

    # Standardize Positions and Faces
    verts = decompose(Point3f, raw_mesh)
    faces_idx = decompose(TriangleFace{Int}, raw_mesh)

    return GeometryBasics.Mesh(verts, faces_idx)
end
