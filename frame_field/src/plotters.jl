using CairoMakie
using GeometryBasics
using Statistics  


############################################################
############################################################
############################################################
############################################################
############################################################
############################################################


"""
    plot_triangulation(mesh::GeometryBasics.Mesh; savepath=nothing, show_edges=true)

Plot a 2D projection of a triangular `GeometryBasics.Mesh` using CairoMakie.
If `savepath` is provided the figure will be saved to that path.
"""
function plot_triangulation(mesh::GeometryBasics.Mesh; savepath=nothing, show_edges=true,
                             show_dual=true, dual_color=:orange, dual_width=1.0,
                             edge_color=:black, vertex_color=:red, face_color=(0.7,0.85,1.0,0.9))
    fig = Figure(Scene=(900, 700))
    ax = Axis(fig[1, 1], aspect=DataAspect())

    # Extract vertex coordinates as matrix Nx3 using GeometryBasics helpers
    vs = coordinates(mesh)
    V = zeros(Float64, length(vs), 3)
    for (i, p) in enumerate(vs)
        V[i, 1] = p[1]
        V[i, 2] = p[2]
        V[i, 3] = p[3]
    end

    # Iterate faces and draw polygons (project to XY)
    fs = faces(mesh)

    # Build edge->faces map (used for dual graph and Euler characteristic)
    edge_map = Dict{Tuple{Int,Int}, Vector{Int}}()
    for (fi, f) in enumerate(fs)
        idxs = Tuple(f)
        edges = (
            (min(idxs[1], idxs[2]), max(idxs[1], idxs[2])),
            (min(idxs[2], idxs[3]), max(idxs[2], idxs[3])),
            (min(idxs[3], idxs[1]), max(idxs[3], idxs[1]))
        )
        for e in edges
            push!(get!(edge_map, e, Int[]), fi)
        end
    end

    for f in fs
        # f is a TriangleFace with indices (1-based)
        idxs = Tuple(f)
        pts = Point2f[]
        for i in idxs
            push!(pts, Point2f(V[i, 1], V[i, 2]))
        end
        poly!(ax, pts, color=face_color, strokecolor=edge_color, strokewidth=0.5)
        if show_edges
            lines!(ax, [p[1] for p in pts], [p[2] for p in pts], color=edge_color, linewidth=0.3)
        end
    end

    # Compute face centroids (for dual graph) and draw triangles
    fs = faces(mesh)

    # Scatter vertices on top
    scatter!(ax, V[:, 1], V[:, 2], color=:red, markersize=6)

    # If requested, compute and draw dual graph (centroids connected across edges)
    if show_dual
        nfaces = length(fs)
        centroids = zeros(Float64, nfaces, 2)
        for (fi, f) in enumerate(fs)
            idxs = Tuple(f)
            centroids[fi, 1] = mean(V[[idxs...], 1])
            centroids[fi, 2] = mean(V[[idxs...], 2])
        end

        # Draw dual edges for interior edges that have two adjacent faces
        for (_, faces_on_edge) in edge_map
            if length(faces_on_edge) == 2
                f1, f2 = faces_on_edge[1], faces_on_edge[2]
                xcoords = [centroids[f1, 1], centroids[f2, 1]]
                ycoords = [centroids[f1, 2], centroids[f2, 2]]
                lines!(ax, xcoords, ycoords, color=dual_color, linewidth=dual_width, linestyle=:dash)
            end
        end

        # Optional: draw centroid markers
        scatter!(ax, centroids[:, 1], centroids[:, 2], color=dual_color, markersize=4, alpha=0.9)
    end

        # --- Legend and Euler characteristic ---
        # Create small invisible plot objects to represent legend entries
        dummy_edge = lines!(ax, [0.0, 0.0], [0.0, 0.0], color=edge_color, linewidth=1.0, visible=false)
        dummy_dual = lines!(ax, [0.0, 0.0], [0.0, 0.0], color=dual_color, linewidth=dual_width, linestyle=:dash, visible=false)
        dummy_vertex_edge = scatter!(ax, [0.0], [0.0], color=dual_color, markersize=4, visible=false)
        dummy_vertex = scatter!(ax, [0.0], [0.0], color=vertex_color, markersize=6, visible=false)



        # Compute Euler characteristic: V - E + F
        # E is number of unique edges in edge_map if it exists, otherwise approximate
        # Count unique edges from the previously-built edge_map
        Ecount = length(edge_map)

        Vcount = size(V, 1)
        Fcount = length(fs)
        euler_char = Vcount - Ecount + Fcount
        euler_char_vertex = scatter!(ax, [0.0], [0.0], color=:black, markersize=0, visible=false)  # dummy for legend

        # Place legend to the right of the axis
        legend = Legend(fig[1, 2], [dummy_edge, dummy_dual, dummy_vertex_edge, dummy_vertex, euler_char_vertex], ["Edge", "Dual edge", "Dual vertex", "Vertex", "χ = $(euler_char)"])
        fig[1, 2] = legend

    if savepath !== nothing
        save(savepath, fig)
    end
    return fig
end



############################################################
############################################################
############################################################
############################################################
############################################################
############################################################


"""
    plot_forest(mesh::GeometryBasics.Mesh; constrained_faces=Int[], 
                potential_fixed_edges=nothing,
                savepath=nothing, edge_color=:black, constrained_color=:red, 
                tree_color=:green,
                face_color=(0.95,0.97,1.0,0.95))

Plot the triangular mesh and highlight constrained faces and the spanning forest.
- Triangulation edges: `edge_color` (default black)
- Constrained face edges: `constrained_color` (default red)
- Tree (dual) edges corresponding to `potential_fixed_edges`: `tree_color` (default green)
- Faces colored by their parent tree component

If `potential_fixed_edges` is `nothing`, the function will call `compute_spanning_forest(mesh; constrained_faces=...)`.
"""
function plot_forest(mesh::GeometryBasics.Mesh; constrained_faces=Int[], 
                     potential_fixed_edges=nothing,
                     savepath=nothing, edge_color=:black, constrained_color=:red, 
                     tree_color=:green,
                     face_color=(0.95,0.97,1.0,0.95),
                     display_background=true)
    fig = Figure(Scene=(900,700))
    ax = Axis(fig[1,1], aspect=DataAspect())

    # Extract vertex coordinates
    vs = coordinates(mesh)
    V = zeros(Float64, length(vs), 3)
    for (i, p) in enumerate(vs)
        V[i,1] = p[1]
        V[i,2] = p[2]
        V[i,3] = p[3]
    end

    fs = faces(mesh)

    # Build edge->faces map
    edge_map = Dict{Tuple{Int,Int}, Vector{Int}}()
    for (fi, f) in enumerate(fs)
        idxs = Tuple(f)
        edges = (
            (min(idxs[1], idxs[2]), max(idxs[1], idxs[2])),
            (min(idxs[2], idxs[3]), max(idxs[2], idxs[3])),
            (min(idxs[3], idxs[1]), max(idxs[3], idxs[1]))
        )
        for e in edges
            push!(get!(edge_map, e, Int[]), fi)
        end
    end

    # Compute potential fixed edges if not provided
    if potential_fixed_edges === nothing
        # Import the function from dijkstra_forest module
        include("dijkstra_forest.jl")
        potential_fixed_edges = compute_spanning_forest(mesh; constrained_faces=constrained_faces)
    end

    # Build dual adjacency from edge_map and potential_fixed_edges
    nfaces = length(fs)
    dual_adj = [Vector{Int}() for _ in 1:nfaces]
    for (key, faces_on_edge) in edge_map
        if length(faces_on_edge) == 2 && (key in potential_fixed_edges)
            f1, f2 = faces_on_edge[1], faces_on_edge[2]
            push!(dual_adj[f1], f2)
            push!(dual_adj[f2], f1)
        end
    end

    # Assign each face to a tree component (root node)
    face_to_root = zeros(Int, nfaces)
    for cf in constrained_faces
        if 1 <= cf <= nfaces
            # BFS from this constrained face
            queue = [cf]
            visited_local = Set{Int}()
            while !isempty(queue)
                cur = popfirst!(queue)
                if cur in visited_local || face_to_root[cur] != 0
                    continue
                end
                push!(visited_local, cur)
                face_to_root[cur] = cf
                for nbr in dual_adj[cur]
                    if !(nbr in visited_local) && face_to_root[nbr] == 0
                        push!(queue, nbr)
                    end
                end
            end
        end
    end

    # Assign colors to each root (tree component)
    unique_roots = unique(filter(x -> x != 0, face_to_root))
    n_roots = length(unique_roots)
    
    # Generate distinct colors using a colormap
    if n_roots > 0
        colormap = cgrad(:tab20, n_roots, categorical=true)
        root_colors = Dict(root => colormap[i] for (i, root) in enumerate(unique_roots))
    else
        root_colors = Dict{Int, Any}()
    end

    # Compute centroids for dual edges
    centroids = zeros(Float64, nfaces, 2)
    for (fi, f) in enumerate(fs)
        idxs = Tuple(f)
        centroids[fi,1] = mean(V[[idxs...], 1])
        centroids[fi,2] = mean(V[[idxs...], 2])
    end

    # Draw faces with colors based on parent node
    for (fi, f) in enumerate(fs)
        idxs = Tuple(f)
        pts = Point2f[]
        for i in idxs
            push!(pts, Point2f(V[i,1], V[i,2]))
        end
        
        # Determine face color based on root
        if face_to_root[fi] != 0 && display_background
            fcolor = root_colors[face_to_root[fi]]
        else
            fcolor = face_color  # default color for unconnected faces
        end
        
        poly!(ax, pts, color=fcolor, strokecolor=edge_color, strokewidth=0.5)
        lines!(ax, [p[1] for p in pts], [p[2] for p in pts], color=edge_color, linewidth=0.6)
    end

    # Highlight constrained face edges in red
    for cf in constrained_faces
        if 1 <= cf <= nfaces
            f = fs[cf]
            idxs = Tuple(f)
            # Draw each edge of the constrained triangle
            for k in 1:3
                a = idxs[k]
                b = idxs[mod1(k+1, 3)]
                xcoords = [V[a,1], V[b,1]]
                ycoords = [V[a,2], V[b,2]]
                lines!(ax, xcoords, ycoords, color=constrained_color, linewidth=2.5)
            end
        end
    end

    # Draw tree (dual) edges in green for potential_fixed_edges
    for (key, faces_on_edge) in edge_map
        if length(faces_on_edge) == 2 && (key in potential_fixed_edges)
            f1, f2 = faces_on_edge[1], faces_on_edge[2]
            xcoords = [centroids[f1, 1], centroids[f2, 1]]
            ycoords = [centroids[f1, 2], centroids[f2, 2]]
            lines!(ax, xcoords, ycoords, color=tree_color, linewidth=2.0)
        end
    end

    # Draw centroid markers for faces participating in tree
    participates = falses(nfaces)
    for (key, faces_on_edge) in edge_map
        if length(faces_on_edge) == 2 && (key in potential_fixed_edges)
            participates[faces_on_edge[1]] = true
            participates[faces_on_edge[2]] = true
        end
    end
    scatter!(ax, centroids[participates,1], centroids[participates,2], color=tree_color, markersize=4)

    # Legend
    dummy_edge = lines!(ax, [0.0,0.0], [0.0,0.0], color=edge_color, linewidth=1.0, visible=false)
    dummy_constrained = lines!(ax, [0.0,0.0], [0.0,0.0], color=constrained_color, linewidth=2.5, visible=false)
    dummy_tree = lines!(ax, [0.0,0.0], [0.0,0.0], color=tree_color, linewidth=2.0, visible=false)
    legend = Legend(fig[1,2], [dummy_edge, dummy_constrained, dummy_tree], 
                    ["Triangulation edge", "Constrained face edge", "Tree (dual) edge"]) 
    fig[1,2] = legend

    if savepath !== nothing
        save(savepath, fig)
    end
    return fig
end

"""
    plot_cross_field(
        mesh::GeometryBasics.Mesh,
        theta::Vector{Float64};
        scale::Float64=0.3,
        savepath=nothing,
        edge_color=:black,
        field_color=:blue,
        field_width=2.0,
        constrained_faces::Vector{Int}=Int[],
        constrained_angles::Vector{Float64}=Float64[],
        constraint_color=:red,
        constraint_width=3.0,
        boundary_faces::Vector{Int}=Int[],
        boundary_color=:orange,
        boundary_width=2.5,
        subsample::Int=1
    )

Plot the computed cross field (4-RoSy field) on a triangular mesh.

Each face displays a cross (4 directions at 90° intervals) based on the computed
angle θ for that face. The cross represents the preferred directions for the
quadrilateral mesh.

Constrained faces are highlighted with their user-specified directions drawn
in a different color. Boundary faces (auto-aligned) shown in a third color.

# Arguments
- `mesh` - Triangular mesh
- `theta` - Vector of face angles (length = number of faces)
- `scale` - Scale factor for cross arm length relative to face size
- `savepath` - Optional path to save the figure
- `edge_color` - Color for mesh edges
- `field_color` - Color for cross field directions
- `field_width` - Line width for cross field vectors
- `constrained_faces` - Indices of faces with user constraints
- `constrained_angles` - Corresponding constraint angles
- `constraint_color` - Color for user constrained face directions
- `constraint_width` - Line width for constraint vectors
- `boundary_faces` - Indices of boundary faces (auto-aligned)
- `boundary_color` - Color for boundary-aligned directions
- `boundary_width` - Line width for boundary vectors
- `subsample` - Plot every nth face (1 = all faces, 2 = every other face, etc.)
"""
function plot_cross_field(
    mesh::GeometryBasics.Mesh,
    theta::Vector{Float64};
    scale::Float64=0.3,
    savepath=nothing,
    edge_color=:black,
    field_color=:blue,
    field_width=2.0,
    constrained_faces::Vector{Int}=Int[],
    constrained_angles::Vector{Float64}=Float64[],
    constraint_color=:red,
    constraint_width=3.0,
    boundary_faces::Vector{Int}=Int[],
    boundary_color=:orange,
    boundary_width=2.5,
    subsample::Int=1
)
    fig = Figure(size=(900, 700))
    ax = Axis(fig[1, 1], aspect=DataAspect(), title="Cross Field Visualization")

    # Extract vertex coordinates
    vs = coordinates(mesh)
    V = zeros(Float64, length(vs), 3)
    for (i, p) in enumerate(vs)
        V[i, 1] = p[1]
        V[i, 2] = p[2]
        V[i, 3] = p[3]
    end

    fs = faces(mesh)
    nfaces = length(fs)

    # Draw mesh faces and edges
    for f in fs
        idxs = Tuple(f)
        pts = Point2f[]
        for i in idxs
            push!(pts, Point2f(V[i, 1], V[i, 2]))
        end
        poly!(ax, pts, color=(0.9, 0.9, 0.9, 0.3), strokecolor=edge_color, strokewidth=0.5)
    end

    # Compute face centroids and characteristic sizes
    centroids = zeros(Float64, nfaces, 2)
    face_sizes = zeros(Float64, nfaces)
    
    for (fi, f) in enumerate(fs)
        idxs = Tuple(f)
        # Centroid
        centroids[fi, 1] = mean(V[[idxs...], 1])
        centroids[fi, 2] = mean(V[[idxs...], 2])
        
        # Characteristic size (average edge length)
        pts = [V[idxs[1], 1:2], V[idxs[2], 1:2], V[idxs[3], 1:2]]
        edge_lengths = [
            norm(pts[2] - pts[1]),
            norm(pts[3] - pts[2]),
            norm(pts[1] - pts[3])
        ]
        face_sizes[fi] = mean(edge_lengths)
    end

    # Create sets for quick lookup
    constrained_set = Set(constrained_faces)
    boundary_set = Set(boundary_faces)
    
    # Draw cross field for each face (with subsampling)
    for fi in 1:nfaces
        if fi > length(theta)
            continue  # Skip if theta not provided for this face
        end
        
        # Determine face type
        is_boundary = fi in boundary_set
        is_constrained = fi in constrained_set
        is_special = is_boundary || is_constrained
        
        # Subsample: skip faces unless it's special or matches subsample rate
        if !is_special && (fi % subsample != 0)
            continue  # Skip this face
        end
        
        cx, cy = centroids[fi, 1], centroids[fi, 2]
        θ = theta[fi]
        arm_length = scale * face_sizes[fi]
        
        # Determine colors and widths based on face type
        if is_constrained
            color = constraint_color
            width = constraint_width
        elseif is_boundary
            color = boundary_color
            width = boundary_width
        else
            color = field_color
            width = field_width
        end
        
        # Draw 4 directions at 90° intervals (cross field / 4-RoSy)
        for k in 0:3
            angle = θ + k * π/2
            dx = arm_length * cos(angle)
            dy = arm_length * sin(angle)
            
            # Draw line from center in both directions
            lines!(ax, 
                [cx - dx, cx + dx], 
                [cy - dy, cy + dy], 
                color=color, 
                linewidth=width)
        end
        
        # Draw a small circle at the centroid
        marker_size = is_special ? 6 : 3
        scatter!(ax, [cx], [cy], color=color, markersize=marker_size, alpha=0.7)
    end

    # Legend
    dummy_edge = lines!(ax, [0.0, 0.0], [0.0, 0.0], color=edge_color, linewidth=0.5, visible=false)
    dummy_field = lines!(ax, [0.0, 0.0], [0.0, 0.0], color=field_color, linewidth=field_width, visible=false)
    
    legend_items = [dummy_edge, dummy_field]
    legend_labels = ["Mesh edges", "Cross field directions"]
    
    # Add boundary entry to legend if there are boundary faces
    if !isempty(boundary_faces)
        dummy_boundary = lines!(ax, [0.0, 0.0], [0.0, 0.0], color=boundary_color, linewidth=boundary_width, visible=false)
        push!(legend_items, dummy_boundary)
        push!(legend_labels, "Boundary aligned ($(length(boundary_faces)) faces)")
    end
    
    # Add constraint entry to legend if there are any user constraints
    if !isempty(constrained_faces)
        dummy_constraint = lines!(ax, [0.0, 0.0], [0.0, 0.0], color=constraint_color, linewidth=constraint_width, visible=false)
        push!(legend_items, dummy_constraint)
        push!(legend_labels, "User constraints ($(length(constrained_faces)) faces)")
    end
    
    legend = Legend(fig[1, 2], legend_items, legend_labels)
    
    if savepath !== nothing
        save(savepath, fig)
        println("Saved cross field plot to: $savepath")
    end
    
    return fig
end


"""
    plot_cross_field_with_singularities(
        mesh, theta, singularities;
        scale=0.25, savepath=nothing, 
        field_color=:darkblue, field_width=1.5,
        zero_rotation_color=:red, zero_rotation_width=2.5,
        boundary_faces=Int[], boundary_color=:orange, boundary_width=2.0,
        show_singularity_indices=true, singularity_marker_size=15,
        subsample=1
    )

Enhanced cross field visualization with singularities highlighted.

# Features
- Highlights arrows corresponding to zero rotations (θ direction) in a different color
- Plots singularities with their index values
- Shows boundary faces in a distinct color

# Arguments
- `mesh`: GeometryBasics.Mesh triangulation
- `theta`: Vector of face angles (length = n_faces)
- `singularities`: Dict mapping vertex index → singularity index
- `scale`: Arrow length scale factor
- `field_color`: Color for regular cross field arms
- `zero_rotation_color`: Color to highlight the zero-rotation (θ) direction
- `boundary_faces`: Indices of boundary faces
- `show_singularity_indices`: Whether to show index numbers on singularities
- `singularity_marker_size`: Size of singularity markers
- `subsample`: Plot every Nth face (1 = all faces)
"""
function plot_cross_field_with_singularities(
    mesh::GeometryBasics.Mesh,
    theta::Vector{Float64},
    singularities::Dict{Int, Float64};
    scale::Float64=0.25,
    savepath=nothing,
    field_color=:darkblue,
    field_width::Float64=1.5,
    zero_rotation_color=:red,
    zero_rotation_width::Float64=2.5,
    boundary_faces::Vector{Int}=Int[],
    boundary_color=:orange,
    boundary_width::Float64=2.0,
    show_singularity_indices::Bool=true,
    singularity_marker_size::Int=15,
    subsample::Int=1,
    edge_color=:gray,
    edge_width::Float64=0.5
)
    fig = Figure(size=(1200, 1000))
    ax = Axis(fig[1, 1], aspect=DataAspect(), title="Cross Field with Singularities")
    
    vs = coordinates(mesh)
    fs = faces(mesh)
    n_faces = length(fs)
    
    # Extract coordinates
    V = zeros(Float64, length(vs), 3)
    for (i, p) in enumerate(vs)
        V[i, 1] = p[1]
        V[i, 2] = p[2]
        V[i, 3] = p[3]
    end
    
    # Draw mesh edges
    for f in fs
        idxs = Tuple(f)
        for i in 1:3
            i1 = idxs[i]
            i2 = idxs[mod1(i+1, 3)]
            lines!(ax, [V[i1,1], V[i2,1]], [V[i1,2], V[i2,2]], 
                   color=edge_color, linewidth=edge_width, alpha=0.3)
        end
    end
    
    # Determine mesh scale for arrow lengths
    x_range = maximum(V[:,1]) - minimum(V[:,1])
    y_range = maximum(V[:,2]) - minimum(V[:,2])
    mesh_scale = min(x_range, y_range)
    arm_length = mesh_scale * scale / sqrt(n_faces)
    
    # Compute centroids
    centroids = zeros(Float64, n_faces, 2)
    for (fi, f) in enumerate(fs)
        idxs = Tuple(f)
        centroids[fi,1] = mean(V[[idxs...], 1])
        centroids[fi,2] = mean(V[[idxs...], 2])
    end
    
    # Plot cross field
    boundary_set = Set(boundary_faces)
    
    for fi in 1:subsample:n_faces
        θ = theta[fi]
        cx, cy = centroids[fi, :]
        
        is_boundary = fi in boundary_set
        
        # Draw 4 directions at 90° intervals
        for k in 0:1
            angle = θ + k * π/2
            dx = arm_length * cos(angle)
            dy = arm_length * sin(angle)
            
            # Highlight the zero-rotation direction (k=0, corresponding to θ)
            if k == 0
                # color = is_boundary ? boundary_color : zero_rotation_color
                color = zero_rotation_color
                width = is_boundary ? boundary_width : zero_rotation_width
            else
                color = is_boundary ? boundary_color : field_color
                width = is_boundary ? boundary_width : field_width
            end
            
            # Draw line from center in both directions
            lines!(ax, 
                [cx, cx + dx], 
                [cy, cy + dy], 
                color=color, 
                linewidth=width,
                alpha=0.8)
        end
        
        # Small marker at centroid
        scatter!(ax, [cx], [cy], 
                color=is_boundary ? boundary_color : field_color, 
                markersize=3, alpha=0.5)
    end
    
    # Plot singularities
    if !isempty(singularities)
        sing_positions = Point2f[]
        sing_indices = Float64[]
        sing_colors = []
        
        for (v_idx, index) in singularities
            if 1 <= v_idx <= length(vs)
                push!(sing_positions, Point2f(V[v_idx, 1], V[v_idx, 2]))
                push!(sing_indices, index)
                
                # Color code by index: positive = red, negative = blue
                if index > 0
                    push!(sing_colors, :red)
                elseif index < 0
                    push!(sing_colors, :blue)
                else
                    push!(sing_colors, :gray)
                end
            end
        end
        
        # Draw singularity markers
        scatter!(ax, sing_positions, 
                color=sing_colors, 
                markersize=singularity_marker_size,
                marker=:circle,
                strokewidth=2,
                strokecolor=:black,
                alpha=0.8)
        
        # Add index labels if requested
        if show_singularity_indices
            for (i, (v_idx, index)) in enumerate(singularities)
                if 1 <= v_idx <= length(vs)
                    x, y = V[v_idx, 1], V[v_idx, 2]
                    # Offset label slightly
                    offset = mesh_scale * 0.02
                    index_str = string(round(index, digits=3))
                    text!(ax, x + offset, y + offset, 
                          text=index_str,
                          fontsize=10,
                          color=:black,
                          font=:bold)
                end
            end
        end
    end
    
    # Legend
    legend_items = []
    legend_labels = String[]
    
    # Zero rotation direction
    dummy_zero = lines!(ax, [0.0, 0.0], [0.0, 0.0], 
                       color=zero_rotation_color, linewidth=zero_rotation_width, visible=false)
    push!(legend_items, dummy_zero)
    push!(legend_labels, "Zero rotation (θ direction)")
    
    # Regular field
    dummy_field = lines!(ax, [0.0, 0.0], [0.0, 0.0], 
                        color=field_color, linewidth=field_width, visible=false)
    push!(legend_items, dummy_field)
    push!(legend_labels, "Other directions (θ + π/2, π, 3π/2)")
    
    # Boundary
    if !isempty(boundary_faces)
        dummy_boundary = lines!(ax, [0.0, 0.0], [0.0, 0.0], 
                               color=boundary_color, linewidth=boundary_width, visible=false)
        push!(legend_items, dummy_boundary)
        push!(legend_labels, "Boundary aligned ($(length(boundary_faces)))")
    end
    
    # Singularities
    if !isempty(singularities)
        n_pos = count(>(0), values(singularities))
        n_neg = count(<(0), values(singularities))
        
        if n_pos > 0
            dummy_pos = scatter!(ax, [0.0], [0.0], color=:red, markersize=singularity_marker_size, 
                               marker=:circle, strokewidth=2, strokecolor=:black, visible=false)
            push!(legend_items, dummy_pos)
            push!(legend_labels, "Positive singularities ($n_pos)")
        end
        
        if n_neg > 0
            dummy_neg = scatter!(ax, [0.0], [0.0], color=:blue, markersize=singularity_marker_size,
                               marker=:circle, strokewidth=2, strokecolor=:black, visible=false)
            push!(legend_items, dummy_neg)
            push!(legend_labels, "Negative singularities ($n_neg)")
        end
    end
    
    # Add total index sum as a label in the legend (create a dummy invisible item for it)
    if !isempty(singularities)
        total_index = sum(values(singularities))
        # Create dummy item for total index label
        dummy_total = lines!(ax, [0.0, 0.0], [0.0, 0.0], color=:transparent, visible=false)
        push!(legend_items, dummy_total)
        push!(legend_labels, "Total index: $(round(total_index, digits=3))")
    end
    
    Legend(fig[1, 2], legend_items, legend_labels)
    
    if savepath !== nothing
        save(savepath, fig)
        println("Saved enhanced cross field plot to: $savepath")
    end
    
    return fig
end


"""
    plot_cut_graph(
        mesh, cut_graph, singularities=Dict{Int,Float64}();
        savepath=nothing,
        cut_color=:red, cut_width=3.0,
        singularity_size=20,
        show_vertex_labels=false
    )

Visualize the cut graph on the mesh.

Shows:
- Mesh edges in gray
- Cut/seam edges in red (thick)
- Singularities as colored circles
- Optional vertex labels
"""
function plot_cut_graph(
    mesh::GeometryBasics.Mesh,
    cut_graph,  # CutGraph struct
    singularities::Dict{Int,Float64}=Dict{Int,Float64}();
    savepath=nothing,
    cut_color=:red,
    cut_width::Float64=3.0,
    singularity_size::Int=20,
    show_vertex_labels::Bool=false,
    mesh_edge_color=:gray,
    mesh_edge_width::Float64=0.5
)
    fig = Figure(size=(1200, 1000))
    ax = Axis(fig[1, 1], aspect=DataAspect(), title="Cut Graph Visualization")
    
    vs = coordinates(mesh)
    fs = faces(mesh)
    
    # Extract coordinates
    V = zeros(Float64, length(vs), 3)
    for (i, p) in enumerate(vs)
        V[i, 1] = p[1]
        V[i, 2] = p[2]
        V[i, 3] = p[3]
    end
    
    # Draw all mesh edges in light gray
    for f in fs
        idxs = Tuple(f)
        for i in 1:3
            i1 = idxs[i]
            i2 = idxs[mod1(i+1, 3)]
            lines!(ax, [V[i1,1], V[i2,1]], [V[i1,2], V[i2,2]], 
                   color=mesh_edge_color, linewidth=mesh_edge_width, alpha=0.3)
        end
    end
    
    # Draw cut/seam edges in red (thick)
    seam_edges = cut_graph.seam_edges
    if !isempty(seam_edges)
        for (v1, v2) in seam_edges
            lines!(ax, [V[v1,1], V[v2,1]], [V[v1,2], V[v2,2]], 
                   color=cut_color, linewidth=cut_width, alpha=0.9)
        end
    end
    
    # Draw singularities
    if !isempty(singularities)
        sing_positions = Point2f[]
        sing_colors = []
        
        for (v_idx, index) in singularities
            if 1 <= v_idx <= length(vs)
                push!(sing_positions, Point2f(V[v_idx, 1], V[v_idx, 2]))
                
                # Color code by index
                if index > 0
                    push!(sing_colors, :blue)
                elseif index < 0
                    push!(sing_colors, :green)
                else
                    push!(sing_colors, :yellow)
                end
            end
        end
        
        scatter!(ax, sing_positions, 
                color=sing_colors, 
                markersize=singularity_size,
                marker=:circle,
                strokewidth=2,
                strokecolor=:black,
                alpha=0.8)
        
        # Add labels
        for (v_idx, index) in singularities
            if 1 <= v_idx <= length(vs)
                x, y = V[v_idx, 1], V[v_idx, 2]
                text!(ax, x, y, 
                      text="v$v_idx\n$(round(index, digits=2))",
                      fontsize=9,
                      color=:black,
                      align=(:center, :center),
                      font=:bold)
            end
        end
    end
    
    # Optionally show vertex labels
    if show_vertex_labels
        for i in 1:length(vs)
            x, y = V[i, 1], V[i, 2]
            text!(ax, x, y, 
                  text="$i",
                  fontsize=7,
                  color=:darkgray,
                  align=(:center, :center))
        end
    end
    
    # Legend
    legend_items = []
    legend_labels = String[]
    
    # Mesh edges
    dummy_mesh = lines!(ax, [0.0, 0.0], [0.0, 0.0], 
                       color=mesh_edge_color, linewidth=mesh_edge_width, visible=false)
    push!(legend_items, dummy_mesh)
    push!(legend_labels, "Mesh edges")
    
    # Cut edges
    if !isempty(seam_edges)
        dummy_cut = lines!(ax, [0.0, 0.0], [0.0, 0.0], 
                          color=cut_color, linewidth=cut_width, visible=false)
        push!(legend_items, dummy_cut)
        push!(legend_labels, "Cut/Seam edges ($(length(seam_edges)))")
    end
    
    # Singularities
    if !isempty(singularities)
        n_pos = count(>(0), values(singularities))
        n_neg = count(<(0), values(singularities))
        
        if n_pos > 0
            dummy_pos = scatter!(ax, [0.0], [0.0], color=:blue, markersize=singularity_size,
                               marker=:circle, strokewidth=2, strokecolor=:black, visible=false)
            push!(legend_items, dummy_pos)
            push!(legend_labels, "Positive singularities ($n_pos)")
        end
        
        if n_neg > 0
            dummy_neg = scatter!(ax, [0.0], [0.0], color=:green, markersize=singularity_size,
                               marker=:circle, strokewidth=2, strokecolor=:black, visible=false)
            push!(legend_items, dummy_neg)
            push!(legend_labels, "Negative singularities ($n_neg)")
        end
        
        total_index = sum(values(singularities))
        dummy_total = lines!(ax, [0.0, 0.0], [0.0, 0.0], color=:transparent, visible=false)
        push!(legend_items, dummy_total)
        push!(legend_labels, "Total index: $(round(total_index, digits=3))")
    end
    
    Legend(fig[1, 2], legend_items, legend_labels)
    
    if savepath !== nothing
        save(savepath, fig)
        println("Saved cut graph plot to: $savepath")
    end
    
    return fig
end
