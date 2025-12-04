using GeometryBasics
using WGLMakie


# A helper struct to handle mesh topology explicitly
struct MeshTopology
    # Maps (v1, v2) sorted tuple to a unique edge index
    edge_map::Dict{Tuple{Int, Int}, Int} 
    # Maps face index to list of (neighbor_face_index, shared_edge_index)
    dual_adj::Vector{Vector{Tuple{Int, Int}}} 
    # We need this to iterate edges in the solver later
    edge_to_faces::Dict{Int, Vector{Int}}
    n_edges::Int
end

"""
    compute_face_centroids(mesh::Mesh)

Returns a Vector{Point3f} containing the centroid of every triangle in the mesh.
"""
function compute_face_centroids(mesh::Mesh)
    # Get all faces (triangles) and vertices
    fs = faces(mesh)
    vs = coordinates(mesh)
    
    centroids = Vector{Point3f}(undef, length(fs))
    
    for (i, face) in enumerate(fs)
        # GeometryBasics faces are 1-based indices into vertices
        p1 = vs[face[1]]
        p2 = vs[face[2]]
        p3 = vs[face[3]]
        
        # Centroid is just the average of the 3 vertices
        centroids[i] = (p1 + p2 + p3) / 3f0
    end
    
    return centroids
end


"""
    get_dual_graph_segments(topo::MeshTopology, centroids::Vector{Point3f})

Returns a vector of point pairs (p1, p2) representing the edges of the dual graph.
Intended for use with `linesegments!`.
"""
function get_dual_graph_segments(topo::MeshTopology, centroids::Vector{Point3f})
    segments = Point3f[]
    
    # We iterate through the dual adjacency list
    # topo.dual_adj[i] contains a list of neighbors for face i
    for f_idx in 1:length(topo.dual_adj)
        neighbors = topo.dual_adj[f_idx]
        
        p_start = centroids[f_idx]
        
        for (neighbor_idx, edge_idx) in neighbors
            # To avoid drawing every edge twice (once from A->B and B->A),
            # we only add the segment if f_idx < neighbor_idx
            if f_idx < neighbor_idx
                p_end = centroids[neighbor_idx]
                push!(segments, p_start)
                push!(segments, p_end)
            end
        end
    end
    
    return segments
end


# function visualize_mesh(mesh::MeshTopology)
#     fig = Figure(size=(800, 600))
#     ax = Axis3(fig[1, 1], title="Mesh Visualization", aspect=:data)
    
#     # loop over the edge of the mesh, give a rnadom color and plot
#     mesh_edges = 


# end











function visualize_forest(mesh::Mesh, topo::MeshTopology, fixed_edges::BitVector, constrained_faces::Vector{Int})
    # 1. Compute Geometry
    centroids = compute_face_centroids(mesh)
    dual_segments = get_dual_graph_segments(topo, centroids)
    
    # 2. Setup Figure
    fig = Figure(size=(1000, 800))
    ax = Axis3(fig[1, 1], aspect = :data, title = "Dual Graph & Constrained Roots")#, azimuth = 0.0, elevation = pi/2)
    
    # 3. First, find connected components (trees) in the fixed edges
    parent = collect(1:length(centroids))
    
    function find_root(x)
        if parent[x] != x
            parent[x] = find_root(parent[x])
        end
        return parent[x]
    end
    
    function union!(x, y)
        root_x = find_root(x)
        root_y = find_root(y)
        if root_x != root_y
            parent[root_y] = root_x
        end
    end
    
    # Build union-find structure for fixed edges
    for f_idx in 1:length(topo.dual_adj)
        for (neighbor_idx, edge_idx) in topo.dual_adj[f_idx]
            if fixed_edges[edge_idx]
                union!(f_idx, neighbor_idx)
            end
        end
    end
    
    # Assign unique color to each tree
    tree_ids = Dict{Int, Int}()
    color_idx = 0
    for i in 1:length(centroids)
        root = find_root(i)
        if !haskey(tree_ids, root)
            color_idx += 1
            tree_ids[root] = color_idx
        end
    end
    
    # Generate distinct colors
    colors = Makie.wong_colors()
    default_color = RGBf(0.8, 0.8, 0.8)
    
    # Get mesh faces and vertices
    fs = faces(mesh)
    vs = coordinates(mesh)
    
    # Group faces by tree
    tree_faces = Dict{Int, Vector{Int}}()
    non_tree_faces = Int[]
    
    for i in 1:length(centroids)
        root = find_root(i)
        if haskey(tree_ids, root)
            tree_id = tree_ids[root]
            if !haskey(tree_faces, tree_id)
                tree_faces[tree_id] = Int[]
            end
            push!(tree_faces[tree_id], i)
        else
            push!(non_tree_faces, i)
        end
    end
    
    # 4. Plot non-tree faces in gray
    if !isempty(non_tree_faces)
        non_tree_mesh_faces = [fs[i] for i in non_tree_faces]
        non_tree_mesh = Mesh(vs, non_tree_mesh_faces)
        mesh!(ax, non_tree_mesh, color = default_color, shading = NoShading)
    end
    
    # 5. Plot each tree with solid color
    for (tree_id, face_indices) in tree_faces
        tree_mesh_faces = [fs[i] for i in face_indices]
        tree_mesh = Mesh(vs, tree_mesh_faces)
        color = colors[mod1(tree_id, length(colors))]
        mesh!(ax, tree_mesh, color = RGBf(color.r, color.g, color.b), shading = NoShading)
    end
    
    # 6. Plot the mesh edges as wireframe
    wireframe!(ax, mesh, color = (:black, 0.3), linewidth = 0.5)
    
    
    # 7. Plot the Full Dual Graph (light background)
    linesegments!(ax, dual_segments, color = (:blue, 0.2), linewidth = 1)
    
    # 8. Plot Standard Dual Nodes (Black)
    scatter!(ax, centroids, color = :black, markersize = 5)

    # 9. Highlight Constrained Roots (Green)
    # These are the starting points of the trees
    if !isempty(constrained_faces)
        root_points = centroids[constrained_faces]
        scatter!(ax, root_points, color = :green, markersize = 20, label="User Constrained Faces")
    end

    axislegend(ax)
    return fig
end







function plot_cross_field(mesh::Mesh, thetas::Vector{Float64}, constrained_faces; 
                          constraint_vecs::Union{Nothing, Dict{Int, Vec3f}}=nothing, 
                          sings=nothing,
                          topo::Union{Nothing, MeshTopology}=nothing,
                          p_values::Union{Nothing, Dict{Int, Int}}=nothing)
    # Compute centers and frames
    centers = compute_face_centroids(mesh)
    vs = coordinates(mesh)
    fs = faces(mesh)
    
    cross_segments = Point3f[]
    # scale = 0.5 # Adjust based on mesh size
    scale = 0.03 # Adjust based on mesh size
    scale = 0.01 # Adjust based on mesh size
    
    for (i, face) in enumerate(fs)
        p1, p2, p3 = vs[face[1]], vs[face[2]], vs[face[3]]
        center, xaxis, yaxis = get_local_frame(p1, p2, p3)
        
        # The angle theta_i
        angle = thetas[i]
        
        # Generate 4 vectors for the cross
        for k in 0:3
            rot_ang = angle + k * (pi/2)
            dir = xaxis * cos(rot_ang) + yaxis * sin(rot_ang)
            
            push!(cross_segments, center)
            push!(cross_segments, center + dir * scale)
        end
    end

    fig = Figure(size=(1000, 800))
    ax = Axis3(fig[1, 1], title="Generated Cross Field", aspect=:data)
    
    # Mesh
    # Create color array where each edge gets a unique color based on its index
    n_edges = length(faces(mesh)) * 3  # Approximate number of edges
    edge_colors = [RGBf(rand(), rand(), rand()) for _ in 1:n_edges]
    wireframe!(ax, mesh, color = :grey)#edge_colors)
    
    # Cross Field
    linesegments!(ax, cross_segments, color=:black, linewidth=1)
    
    # Dual Edges Visualization (Period Jumps)
    if !isnothing(topo) && !isnothing(p_values)
        dual_segments_colored = Point3f[]
        dual_colors = RGBAf[]
        
        # Iterate over all edges in topology
        for (e_idx, faces_list) in topo.edge_to_faces
            if length(faces_list) == 2
                f1, f2 = faces_list
                p_val = get(p_values, e_idx, 0)
                
                c1 = centers[f1]
                c2 = centers[f2]
                
                push!(dual_segments_colored, c1)
                push!(dual_segments_colored, c2)
                
                # Color based on p value
                color = if p_val == 0
                    RGBAf(0.5, 0.5, 0.5, 0.1) # Faint gray for 0
                elseif p_val > 0
                    RGBAf(1.0, 0.0, 0.0, 1.0) # Red for positive
                else
                    RGBAf(0.0, 0.0, 1.0, 1.0) # Blue for negative
                end
                push!(dual_colors, color)
                push!(dual_colors, color)
            end
        end
        
        if !isempty(dual_segments_colored)
            linesegments!(ax, dual_segments_colored, color=dual_colors, linewidth=3, label="Dual Edges (P)")
        end
    end

    # User Desired Directions (if provided)
    if !isnothing(constraint_vecs) && !isempty(constraint_vecs)
        constraint_segments = Point3f[]
        for (f_idx, dir_vec) in constraint_vecs
            center = centers[f_idx]
            normalized_dir = normalize(dir_vec)
            
            # Draw the constraint direction as a thicker, colored arrow
            push!(constraint_segments, center)
            push!(constraint_segments, center + normalized_dir * scale * 1.5)
        end
        
        linesegments!(ax, constraint_segments, color=:red, linewidth=4, label="User Direction")
    end
    
    # Highlight Constraints
    if !isempty(constrained_faces)
        scatter!(ax, centers[constrained_faces], color=:red, markersize=15, label="Constraint Faces")
    end

    if !isnothing(sings) && !isempty(sings)
        # Highlight Singularities
        sing_points = Point3f[]
        sing_colors = RGBf[]
        for (v_idx, index) in sings
            # Get vertex position
            v_pos = vs[v_idx]
            push!(sing_points, v_pos)
            if index > 0
                push!(sing_colors, RGBf(0.0, 1.0, 0.0))  # Red for positive index
            else
                push!(sing_colors, RGBf(0.0, 0.0, 1.0))  # Blue for negative index
            end
        end
        scatter!(ax, sing_points, color=sing_colors, markersize=20, label="Singularities")
    end
    
    axislegend(ax)
    return fig
end

function plot_field_smoothness(mesh::Mesh, thetas::Vector{Float64}, topo::MeshTopology, kappas::Dict{Int, Float64}, p_values::Dict{Int, Int})
    # Compute smoothness score for each face
    n_faces = length(faces(mesh))
    smoothness_scores = zeros(Float64, n_faces)
    
    for f_idx in 1:n_faces
        max_err = 0.0
        
        # Iterate over neighbors
        for (neighbor_idx, edge_idx) in topo.dual_adj[f_idx]
            # Determine direction
            f1, f2 = topo.edge_to_faces[edge_idx]
            
            # Retrieve kappa and p
            k_val = get(kappas, edge_idx, 0.0)
            p_val = get(p_values, edge_idx, 0)
            
            # Calculate mismatch based on direction
            # Equation: theta_j - theta_i + kappa + p*pi/2 = 0
            # If f_idx is f1 (i), neighbor is f2 (j) -> err = theta_j - theta_i + k + p*pi/2
            # If f_idx is f2 (j), neighbor is f1 (i) -> err = theta_i - theta_j - k - p*pi/2
            
            err = 0.0
            if f_idx == f1
                # i -> j
                # Solver enforces: theta_j - theta_i = k + p*pi/2
                # Residual: theta_j - theta_i - k - p*pi/2
                err = thetas[neighbor_idx] - thetas[f_idx] - k_val - p_val * (pi/2)
            else
                # j -> i
                # Residual: theta_i - theta_j + k + p*pi/2
                err = thetas[neighbor_idx] - thetas[f_idx] + k_val + p_val * (pi/2)
            end
            
            max_err = max(max_err, abs(err))
        end
        smoothness_scores[f_idx] = max_err
    end
    
    # Normalize scores for coloring
    max_score = maximum(smoothness_scores)
    if max_score == 0
        max_score = 1.0 # Avoid division by zero
    end
    
    # Generate colors
    # White (0) -> Red (Max)
    # face_colors = [RGBAf(1.0, 1.0 - (s/max_score), 1.0 - (s/max_score), 1.0) for s in smoothness_scores]
    
    # Create a new mesh with duplicated vertices for per-face coloring
    vs = coordinates(mesh)
    fs = faces(mesh)
    
    new_vs = Point3f[]
    new_fs = TriangleFace{Int}[]
    new_values = Float64[]
    
    for (i, face) in enumerate(fs)
        # Get original vertices
        p1 = vs[face[1]]
        p2 = vs[face[2]]
        p3 = vs[face[3]]
        
        # Add to new list
        push!(new_vs, p1)
        push!(new_vs, p2)
        push!(new_vs, p3)
        
        # Add new face (indices are 1-based, current length is end)
        base_idx = length(new_vs) - 3
        push!(new_fs, TriangleFace(base_idx+1, base_idx+2, base_idx+3))
        
        # Add value for each vertex of this face
        s = smoothness_scores[i]
        push!(new_values, s)
        push!(new_values, s)
        push!(new_values, s)
    end
    
    new_mesh = Mesh(new_vs, new_fs)
    
    fig = Figure(size=(1000, 800))
    ax = Axis3(fig[1, 1], title="Field Smoothness (White=Smooth, Red=Error)", aspect=:data)
    
    # Define colormap: White -> Red
    cmap = cgrad([:white, :red])
    
    # Plot mesh with values mapped to colormap
    m = mesh!(ax, new_mesh, color=new_values, colormap=cmap, colorrange=(0, max_score), shading=NoShading)
    wireframe!(ax, mesh, color=(:black, 0.1), linewidth=0.5)
    
    # Add Colorbar
    Colorbar(fig[1, 2], m, label="Smoothness Error (Residual)")
    
    return fig
end