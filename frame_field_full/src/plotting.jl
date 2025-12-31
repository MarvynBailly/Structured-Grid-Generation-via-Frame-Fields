module Plotting

using ..Types
using ..Analysis
using CairoMakie
using GLMakie
using GeometryBasics
using Printf
using LinearAlgebra

export visualize, plot_global_rotations, plot_panel_topology!, plot_panel_field!, plot_panel_quads!

# --- Helper: Convert Local Theta to Global Vector ---
function get_global_vector(topo, f_idx, theta)
    # 1. Get Reference Frame
    tri = topo.faces[f_idx]
    ref_edge = topo.face_ref_edges[f_idx] # (v_start, v_end)
    
    p1 = topo.vertices[ref_edge[1]]
    p2 = topo.vertices[ref_edge[2]]
    
    # Basis X (e1) = Normalized Reference Edge
    e1 = [p2.x - p1.x, p2.y - p1.y, p2.z - p1.z]
    normalize!(e1)
    
    # Normal (n)
    v1 = topo.vertices[tri[1]]
    v2 = topo.vertices[tri[2]]
    v3 = topo.vertices[tri[3]]
    u = [v2.x - v1.x, v2.y - v1.y, v2.z - v1.z]
    v = [v3.x - v1.x, v3.y - v1.y, v3.z - v1.z]
    n = cross(u, v)
    normalize!(n)
    
    # Basis Y (e2) = n x e1
    e2 = cross(n, e1)
    
    # 2. Rotate in Tangent Plane
    v_global = cos(theta) .* e1 .+ sin(theta) .* e2
    return v_global
end

function get_face_center(topo, f_idx)
    tri = topo.faces[f_idx]
    p1 = topo.vertices[tri[1]]
    p2 = topo.vertices[tri[2]]
    p3 = topo.vertices[tri[3]]
    return Point3f((p1.x+p2.x+p3.x)/3, (p1.y+p2.y+p3.y)/3, (p1.z+p2.z+p3.z)/3)
end


# --- Unified Visualization Function ---

"""
    visualize(filename; kw...)
    visualize(filename, topo; kw...)
    visualize(filename, field; kw...)

Unified plotting function. Automatically determines the layout (1, 2, or 3 panels)
based on the provided data.

# Arguments
- `filename::String`: Output path (e.g., "out.png")
- `topo::MeshTopology`: (Optional) The base mesh.
- `field::CrossField`: (Optional) The cross field to visualize.
- `quad_mesh::QuadMesh`: (Optional) The extracted quad mesh.

# Options
- `show_singularities::Bool`: Plot singularity dots (+ red, - cyan). Default `true`.
- `show_constraints::Bool`: Highlight constrained faces and directions. Default `true`.
- `show_period_jumps::Bool`: Color edges based on period jumps. Default `false`.
- `cut_edges`: Collection of edges (tuples) to draw as cuts (red lines).
- `verbose::Bool`: Print status messages.
"""
function visualize(filename::String; 
                   topo::Union{MeshTopology, Nothing}=nothing,
                   field::Union{CrossField, Nothing}=nothing,
                   quad_mesh=nothing,
                   show_singularities::Bool=false,
                   show_period_jumps::Bool=false,
                   show_constraints::Bool=true,
                   cut_edges::Union{Set{Tuple{Int,Int}}, Vector{Tuple{Int,Int}}, Nothing}=nothing,
                   title_prefix::String="",
                   verbose::Bool=false)

    # 1. Resolve Topology
    if !isnothing(field) && isnothing(topo)
        topo = field.topology
    end
    if isnothing(topo)
        error("Topology must be provided explicitly or via a CrossField object.")
    end

    # 2. Determine Layout (1, 2, or 3 panels)
    has_field = !isnothing(field)
    has_quads = !isnothing(quad_mesh)
    
    n_panels = 1
    if has_field; n_panels += 1; end
    if has_quads; n_panels += 1; end
    
    # Width scales with panels
    fig_size = (600 * n_panels, 600)
    fig = Figure(size=fig_size)
    
    current_col = 1
    
    # --- Panel 1: Initial Mesh ---
    ax1 = Axis(fig[1, current_col], aspect=DataAspect(), title="$(title_prefix)Initial Mesh")
    plot_panel_topology!(ax1, topo)
    current_col += 1
    
    # --- Panel 2: Cross Field (if exists) ---
    if has_field
        ax2 = Axis(fig[1, current_col], aspect=DataAspect(), title="Cross Field")
        plot_panel_field!(ax2, field; 
                          show_singularities=show_singularities,
                          show_period_jumps=show_period_jumps,
                          show_constraints=show_constraints,
                          cut_edges=cut_edges)
        
        # Legend for Field
        if show_singularities || show_period_jumps || show_constraints
            axislegend(ax2, position=:rt, merge=true, unique=true)
        end
        current_col += 1
    end
    
    # --- Panel 3: Quad Mesh (if exists) ---
    if has_quads
        ax3 = Axis(fig[1, current_col], aspect=DataAspect(), title="Extracted Quads")
        plot_panel_quads!(ax3, quad_mesh; topo=topo)
        current_col += 1
    end
    
    save(filename, fig)
    if verbose
        println("Saved visualization to $filename")
        if has_field && show_singularities
             sings = compute_singularities(field)
             println("  Singularities plotted: $(length(sings))")
        end
    end
end

# Convenience dispatches
visualize(filename::String, topo::MeshTopology; kwargs...) = visualize(filename; topo=topo, kwargs...)
visualize(filename::String, field::CrossField; kwargs...) = visualize(filename; field=field, kwargs...)

# Backwards compatibility alias
const plot_results = visualize


# --- Panel Implementations ---

function plot_panel_topology!(ax, topo::MeshTopology)
    # Draw wireframe of input mesh
    for tri in topo.faces
        pts = [topo.vertices[i] for i in tri]; push!(pts, pts[1])
        lines!(ax, [p.x for p in pts], [p.y for p in pts], color=(:gray, 0.4), linewidth=0.5)
    end
end

function plot_panel_field!(ax, field::CrossField; 
                           show_singularities=true, 
                           show_period_jumps=false, 
                           show_constraints=true,
                           cut_edges=nothing)
    
    topo = field.topology
    
    # 1. Background Mesh (Faint)
    for tri in topo.faces
        pts = [topo.vertices[i] for i in tri]; push!(pts, pts[1])
        lines!(ax, [p.x for p in pts], [p.y for p in pts], color=(:gray, 0.1), linewidth=0.3)
    end
    
    # 2. Constraints (Faces & Directions)
    if show_constraints && !isempty(field.constrained_faces)
        c_xs, c_ys, c_us, c_vs = Float64[], Float64[], Float64[], Float64[]
        scale = 0.065
        
        for (f_idx, angle) in field.constrained_faces
            tri = topo.faces[f_idx]
            pts = [topo.vertices[i] for i in tri]
            poly!(ax, [Point2f(p.x, p.y) for p in pts], color=(:orange, 0.4))
            
            # Constraint Arrow
            v_global = get_global_vector(topo, f_idx, angle)
            c = get_face_center(topo, f_idx)
            push!(c_xs, c[1]); push!(c_ys, c[2])
            push!(c_us, v_global[1]*scale); push!(c_vs, v_global[2]*scale)
        end
        
        if !isempty(c_xs)
            arrows2d!(ax, c_xs, c_ys, c_us, c_vs, color=(:red, 0.6), label="Constraint")
        end
    end
    
    # 3. Field Vectors (Sampled)
    n_faces = length(topo.faces)
    sample_rate = if n_faces < 1000; 1
        elseif n_faces < 10000; 2
        elseif n_faces < 50000; 5
        else; 10; end
        
    xs, ys, us, vs = Float64[], Float64[], Float64[], Float64[]
    xs_sec, ys_sec, us_sec, vs_sec = Float64[], Float64[], Float64[], Float64[]
    scale = 0.03
    
    for i in 1:n_faces
        if i % sample_rate != 0; continue; end
        
        # Primary vector
        v_global = get_global_vector(topo, i, field.theta[i])
        c = get_face_center(topo, i)
        
        push!(xs, c[1]); push!(ys, c[2])
        push!(us, v_global[1]*scale); push!(vs, v_global[2]*scale)
        
        # Cross vectors (90, 180, 270)
        for k in 1:3
            v_ortho = get_global_vector(topo, i, field.theta[i] + k*π/2)
            push!(xs_sec, c[1]); push!(ys_sec, c[2])
            push!(us_sec, v_ortho[1]*scale); push!(vs_sec, v_ortho[2]*scale)
        end
    end
    
    arrows2d!(ax, xs_sec, ys_sec, us_sec, vs_sec, color=(:blue))
    arrows2d!(ax, xs, ys, us, vs, color=:blue, label="Frame")

    # 4. Period Jumps
    if show_period_jumps && !isempty(field.period_jumps)
         max_p = isempty(values(field.period_jumps)) ? 1 : maximum(abs.(values(field.period_jumps)))
         
         for ((i, j), p) in field.period_jumps
             if p == 0; continue; end
             
             tri_i = topo.faces[i]; tri_j = topo.faces[j]
             ci = get_face_center(topo, i); cj = get_face_center(topo, j)
             
             color_val = abs(p) / max(max_p, 1)
             edge_color = abs(p) == 1 ? (:yellow, 0.8) : (:red, min(0.9, 0.5 + 0.4 * color_val))
             lw = abs(p) == 1 ? 2.0 : 3.0
             
             lines!(ax, [ci[1], cj[1]], [ci[2], cj[2]], color=edge_color, linewidth=lw)
         end
         # Dummy lines for legend
         lines!(ax, [NaN], [NaN], color=(:yellow, 0.8), linewidth=2, label="|p|=1")
         lines!(ax, [NaN], [NaN], color=(:red, 0.9), linewidth=3, label="|p|≥2")
    end
    
    # 5. Cut Edges
    if !isnothing(cut_edges) && !isempty(cut_edges)
        cx, cy = Float64[], Float64[]
        for (u, v) in cut_edges
            p1, p2 = topo.vertices[u], topo.vertices[v]
            push!(cx, p1.x, p2.x, NaN)
            push!(cy, p1.y, p2.y, NaN)
        end
        lines!(ax, cx, cy, color=:red, linewidth=3.0, label="Cut Graph")
    end
    
    # 6. Singularities
    if show_singularities
        sings = compute_singularities(field)
        pos, neg = Point2f[], Point2f[]
        for (v, I) in sings
            p = topo.vertices[v]
            pt = Point2f(p.x, p.y)
            if I > 0; push!(pos, pt); else; push!(neg, pt); end
        end
        if !isempty(pos); scatter!(ax, pos, color=:red, markersize=12, strokecolor=:black, strokewidth=1, label="Sing. (+)"); end
        if !isempty(neg); scatter!(ax, neg, color=:cyan, markersize=12, strokecolor=:black, strokewidth=1, label="Sing. (-)"); end
    end
end

function plot_panel_quads!(ax, quad_mesh; topo=nothing)
    # Draw Background (if provided)
    if !isnothing(topo)
        for tri in topo.faces
            pts = [topo.vertices[i] for i in tri]; push!(pts, pts[1])
            lines!(ax, [p.x for p in pts], [p.y for p in pts], color=(:gray, 0.1), linewidth=0.3)
        end
    end
    
    # Draw Quads (2D Projection)
    q_edges_x, q_edges_y = Float64[], Float64[]
    
    for (v1, v2, v3, v4) in quad_mesh.quads
        pts = [quad_mesh.vertices[v1], quad_mesh.vertices[v2], 
               quad_mesh.vertices[v3], quad_mesh.vertices[v4], quad_mesh.vertices[v1]]
        
        for p in pts
            push!(q_edges_x, p.x)
            push!(q_edges_y, p.y)
        end
        push!(q_edges_x, NaN) # Separator
        push!(q_edges_y, NaN)
    end
    
    lines!(ax, q_edges_x, q_edges_y, color=:black, linewidth=1.5)
end


# --- Supplementary Visualization (Smooth Field) ---

function plot_global_rotations(field::CrossField, rotations::Vector{Int}; 
                               save_path::String=nothing, verbose=false, 
                               title="Global Smoothed Field", final_cuts=nothing)
    topo = field.topology
    fig = Figure(size=(1000, 800))
    ax = Axis3(fig[1, 1], aspect=:data, title=title) # Using 3D axis for this specific debug view

    # Wireframe
    pts = [Point3f(v.x, v.y, v.z) for v in topo.vertices]
    tris = [GLMakie.TriangleFace(f[1], f[2], f[3]) for f in topo.faces]
    geo_mesh = GeometryBasics.normal_mesh(pts, tris)
    mesh!(ax, geo_mesh, color=(:gray, 0.1), transparency=true, shading=NoShading)
    wireframe!(ax, geo_mesh, color=(:black), linewidth=0.2)

    # Smoothed Global Vectors
    us, vs = Vec3f[], Vec3f[]
    centers = Point3f[]
    
    # Scale calculation
    avg_len = sum(norm(pts[f[1]] - pts[f[2]]) for f in topo.faces) / length(topo.faces)
    vec_scale = avg_len * 0.4

    for i in 1:length(topo.faces)
        r = rotations[i]
        c = get_face_center(topo, i)
        u = get_global_vector(topo, i, field.theta[i] + r * π/2)
        v = get_global_vector(topo, i, field.theta[i] + r * π/2 + π/2)
        push!(centers, c)
        push!(us, Vec3f(u...)); push!(vs, Vec3f(v...))
    end

    arrows2d!(ax, centers, us, lengthscale=vec_scale, color=:red, label="Global u")
    arrows2d!(ax, centers, vs, lengthscale=vec_scale, color=:blue, label="Global v")
    
    if !isnothing(final_cuts)
        cut_segs = Point3f[]
        for (u, v) in final_cuts
            push!(cut_segs, pts[u]); push!(cut_segs, pts[v])
        end
        linesegments!(ax, cut_segs, color=:red, linewidth=4.0, label="Cuts")
    end

    Legend(fig[1, 2], ax)
    
    if !isnothing(save_path)
        save(save_path, fig)
    end
    display(fig)
end

end