# meshIO.jl
using FileIO
using MeshIO
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



"""
    save_cross_field_to_vtu(filename, mesh, thetas, constraint_vecs, topo, p_map, x_sol)

Saves the cross field results to a VTU file for visualization in ParaView.
Includes:
- Mesh geometry (vertices and triangles)
- Cross field directions (4 vectors per face)
- Constraint markers and desired directions
- Period jumps (singularities)
- Face angles (theta values)
"""
function save_cross_field_to_vtu(filename::String, mesh::Mesh, thetas::Vector{Float64}, 
                                  constraint_vecs::Dict{Int, Vec3f}, 
                                  topo::MeshTopology, p_map::Dict{Int, Int}, x_sol::Vector{Float64})
    
    vs = coordinates(mesh)
    fs = faces(mesh)
    n_faces = length(fs)
    
    # Convert mesh to WriteVTK format
    # Extract vertices as matrix (3 x n_vertices)
    points = zeros(3, length(vs))
    for (i, v) in enumerate(vs)
        points[1, i] = v[1]
        points[2, i] = v[2]
        points[3, i] = v[3]
    end
    
    # Extract connectivity (triangles)
    cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, [f[1], f[2], f[3]]) for f in fs]
    
    # Compute cross field directions at each face
    dir1_x, dir1_y, dir1_z = zeros(n_faces), zeros(n_faces), zeros(n_faces)
    dir2_x, dir2_y, dir2_z = zeros(n_faces), zeros(n_faces), zeros(n_faces)
    dir3_x, dir3_y, dir3_z = zeros(n_faces), zeros(n_faces), zeros(n_faces)
    dir4_x, dir4_y, dir4_z = zeros(n_faces), zeros(n_faces), zeros(n_faces)
    
    constraint_marker = zeros(Int, n_faces)
    constraint_dir_x, constraint_dir_y, constraint_dir_z = zeros(n_faces), zeros(n_faces), zeros(n_faces)
    
    for (i, face) in enumerate(fs)
        p1, p2, p3 = vs[face[1]], vs[face[2]], vs[face[3]]
        _, xaxis, yaxis = get_local_frame(p1, p2, p3)
        
        angle = thetas[i]
        
        # Four cross directions
        for k in 0:3
            rot_ang = angle + k * (pi/2)
            dir = xaxis * cos(rot_ang) + yaxis * sin(rot_ang)
            
            if k == 0
                dir1_x[i], dir1_y[i], dir1_z[i] = dir[1], dir[2], dir[3]
            elseif k == 1
                dir2_x[i], dir2_y[i], dir2_z[i] = dir[1], dir[2], dir[3]
            elseif k == 2
                dir3_x[i], dir3_y[i], dir3_z[i] = dir[1], dir[2], dir[3]
            else
                dir4_x[i], dir4_y[i], dir4_z[i] = dir[1], dir[2], dir[3]
            end
        end
        
        # Mark constraints
        if haskey(constraint_vecs, i)
            constraint_marker[i] = 1
            dir_vec = normalize(constraint_vecs[i])
            constraint_dir_x[i] = dir_vec[1]
            constraint_dir_y[i] = dir_vec[2]
            constraint_dir_z[i] = dir_vec[3]
        end
    end
    
    # Extract period jumps (singularities)
    period_jumps = zeros(n_faces)
    for (e_idx, mat_idx) in p_map
        p_val = round(x_sol[mat_idx])
        if abs(p_val) > 0.5  # Non-zero period jump
            faces_on_edge = topo.edge_to_faces[e_idx]
            if length(faces_on_edge) == 2
                # Mark both adjacent faces with the period jump
                period_jumps[faces_on_edge[1]] += p_val / 2
                period_jumps[faces_on_edge[2]] += p_val / 2
            end
        end
    end
    
    # Write VTU file
    vtk_grid(filename, points, cells) do vtk
        # Cell data (per-triangle)
        vtk["theta"] = thetas
        vtk["direction1"] = (dir1_x, dir1_y, dir1_z)
        vtk["direction2"] = (dir2_x, dir2_y, dir2_z)
        vtk["direction3"] = (dir3_x, dir3_y, dir3_z)
        vtk["direction4"] = (dir4_x, dir4_y, dir4_z)
        vtk["constraint_marker"] = constraint_marker
        vtk["user_direction"] = (constraint_dir_x, constraint_dir_y, constraint_dir_z)
        vtk["period_jump"] = period_jumps
    end
    
    println("Saved cross field to $filename.vtu")
end