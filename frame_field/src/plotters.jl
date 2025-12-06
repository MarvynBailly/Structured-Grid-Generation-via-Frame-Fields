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
