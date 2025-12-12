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


############################################################
############################################################
############################################################
############################################################
############################################################
############################################################


"""
    plot_frame_field(mesh::GeometryBasics.Mesh, angles::Vector{Float64};
                     savepath=nothing, arrow_scale=0.3, 
                     show_triangulation=true, constrained_faces=Int[])

Plot the frame field (cross field) on a triangular mesh.

Each triangle face shows a cross (4-way directional field) based on the angle value.
The cross represents the principal directions at π/2 intervals.

# Arguments
- `mesh::GeometryBasics.Mesh` - The triangular mesh
- `angles::Vector{Float64}` - Angle for each face (in radians)
- `savepath` - Path to save the figure (optional)
- `arrow_scale` - Scale factor for the cross arms (default: 0.3)
- `show_triangulation` - Whether to show the mesh edges (default: true)
- `constrained_faces` - Faces with constrained directions (highlighted in red)
- `singularities` - Vector of (vertex_idx, index_value, position, is_boundary) tuples to display
- `show_boundary_singularities` - Whether to show boundary singularities (default: false)
"""
function plot_frame_field(mesh::GeometryBasics.Mesh, angles::Vector{Float64};
                          savepath=nothing, arrow_scale=0.3,
                          show_triangulation=true, constrained_faces=Int[],
                          constrained_color=:red,
                          singularities=[],
                          show_boundary_singularities=false,
                          title="Frame Field",
                          period_jumps=nothing,
                          topology=nothing)
    fig = Figure(size=(900, 700))
    ax = Axis(fig[1, 1], aspect=DataAspect(), title=title)

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

    # Draw triangulation if requested
    if show_triangulation
        for f in fs
            idxs = Tuple(f)
            pts = Point2f[]
            for i in idxs
                push!(pts, Point2f(V[i, 1], V[i, 2]))
            end
            poly!(ax, pts, color=(:lightblue, 0.3), strokecolor=:gray, strokewidth=0.5)
        end
    end

    # Compute face centroids
    centroids = zeros(Float64, nfaces, 2)
    face_radii = zeros(Float64, nfaces)
    for (fi, f) in enumerate(fs)
        idxs = Tuple(f)
        centroids[fi, 1] = mean(V[[idxs...], 1])
        centroids[fi, 2] = mean(V[[idxs...], 2])
        
        # Compute approximate face radius (for scaling arrows)
        dists = [norm([V[idxs[i], 1] - centroids[fi, 1], 
                      V[idxs[i], 2] - centroids[fi, 2]]) for i in 1:3]
        face_radii[fi] = mean(dists)
    end

    # Draw frame field (crosses)
    for (fi, angle) in enumerate(angles)
        cx, cy = centroids[fi, 1], centroids[fi, 2]
        r = face_radii[fi] * arrow_scale
        
        # Draw 4 directions at π/2 intervals
        for k in 0:3
            θ = angle + k * π/2
            dx = r * cos(θ)
            dy = r * sin(θ)
            
            # Draw line from center outward
            color = (fi in constrained_faces) ? constrained_color : :black
            linewidth = (fi in constrained_faces) ? 2.5 : 1.5
            
            lines!(ax, [cx - dx, cx + dx], [cy - dy, cy + dy], 
                   color=color, linewidth=linewidth)
        end
    end

    # Draw centroids
    scatter!(ax, centroids[:, 1], centroids[:, 2], color=:blue, markersize=3, alpha=0.5)

    # Draw dual graph with period jump coloring
    if period_jumps !== nothing && topology !== nothing
        # Build edge map for dual graph
        edge_map = Dict{Tuple{Int,Int}, Vector{Int}}()
        for (fi, f) in enumerate(fs)
            idxs = Tuple(f)
            for j in 1:3
                v1 = idxs[j]
                v2 = idxs[mod1(j+1, 3)]
                edge = v1 < v2 ? (v1, v2) : (v2, v1)
                if !haskey(edge_map, edge)
                    edge_map[edge] = Int[]
                end
                push!(edge_map[edge], fi)
            end
        end
        
        # Color map for period jumps (0->gray, 1->green, 2->yellow, 3->red)
        p_colors = Dict(0 => :gray, 1 => :green, 2 => :yellow, 3 => :red)
        
        # Draw dual edges colored by period jump
        for (edge, adj_faces) in edge_map
            if length(adj_faces) == 2
                face_i, face_j = adj_faces[1], adj_faces[2]
                if face_i > face_j
                    face_i, face_j = face_j, face_i
                end
                
                # Get period jump value
                p_val = get(period_jumps, (face_i, face_j), 0)
                color = get(p_colors, Int(round(p_val)), :black)
                
                # Draw edge between centroids
                x_coords = [centroids[face_i, 1], centroids[face_j, 1]]
                y_coords = [centroids[face_i, 2], centroids[face_j, 2]]
                lines!(ax, x_coords, y_coords, color=color, linewidth=2.5, alpha=0.7)
            end
        end
        
        # Add legend for period jumps
        legend_elements = [
            LineElement(color=:gray, linewidth=2.5),
            LineElement(color=:green, linewidth=2.5),
            LineElement(color=:yellow, linewidth=2.5),
            LineElement(color=:red, linewidth=2.5)
        ]
        legend_labels = ["p=0", "p=1", "p=2", "p=3"]
        Legend(fig[1, 2], legend_elements, legend_labels, "Period Jumps")
    end

    # Highlight constrained faces
    if !isempty(constrained_faces)
        for cf in constrained_faces
            if 1 <= cf <= nfaces
                f = fs[cf]
                idxs = Tuple(f)
                pts = Point2f[]
                for i in idxs
                    push!(pts, Point2f(V[i, 1], V[i, 2]))
                end
                poly!(ax, pts, color=(:red, 0.2), strokecolor=constrained_color, strokewidth=2.5)
            end
        end
    end

    # Draw singularities
    if !isempty(singularities)
        for sing in singularities
            # Handle both 3-tuple (v_idx, index, pos) and 4-tuple (v_idx, index, pos, is_boundary)
            v_idx = sing[1]
            index = sing[2]
            pos = sing[3]
            is_boundary = length(sing) >= 4 ? sing[4] : false
            
            # Skip boundary singularities unless requested
            if is_boundary && !show_boundary_singularities
                continue
            end
            
            x, y = pos[1], pos[2]
            
            # Color based on index sign: positive (green) = lower valence, negative (purple) = higher valence
            if is_boundary
                sing_color = :orange
                marker = :rect
            elseif index > 0
                sing_color = :green
                marker = :star5
            else
                sing_color = :purple
                marker = :diamond
            end
            
            # Size based on absolute index
            markersize = 12 + abs(index) * 20
            
            scatter!(ax, [x], [y], color=sing_color, markersize=markersize, 
                     marker=marker, strokecolor=:black, strokewidth=1.5)
            
            # Add text label with index as fraction
            index_quarters = round(Int, index * 4)
            if abs(index - index_quarters/4) < 1e-6
                index_str = index_quarters > 0 ? "+$(index_quarters)/4" : "$(index_quarters)/4"
            else
                index_str = "$(round(index, digits=2))"
            end
            text!(ax, x, y + 0.03 * (maximum(V[:, 2]) - minimum(V[:, 2])), 
                  text=index_str, fontsize=10, align=(:center, :bottom), color=:black)
        end
    end

    if savepath !== nothing
        save(savepath, fig)
    end
    return fig
end
