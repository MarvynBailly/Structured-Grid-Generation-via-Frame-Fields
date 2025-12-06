# Frame Fields
In this folder, we have functions to:
- read in a triangulation and user defined constraints
- Compute the feasable set of period jumps that may be set to zero using the Dijkstra tree method
- For each face, pick a period junmp to fix to zero out of the feasable set
- The remaining period jumps and unconstrained angles are set up in a matrix to be passed to the greedy solver



## Loading the Triangulation
We load a `.msh` file into Julia's GeometryBasics format:

```julia
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

    return Mesh(verts, faces_idx)
end
```