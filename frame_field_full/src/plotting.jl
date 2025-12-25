module Plotting

using ..Types
using ..Analysis
using CairoMakie
using Printf

export plot_frame, plot_results, plot_smooth_global_field, plot_quad_mesh

function plot_frame(field, filename, title_text; cut_edges=nothing, show_singularities=false)
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], aspect=DataAspect(), title=title_text)
    
    for tri in field.topology.faces
        pts = [field.topology.vertices[i] for i in tri]; push!(pts, pts[1])
        lines!(ax, [p.x for p in pts], [p.y for p in pts], color=(:gray, 0.2), linewidth=0.3)
    end
    
    for (face_idx, _) in field.constrained_faces
        tri = field.topology.faces[face_idx]
        pts = [field.topology.vertices[i] for i in tri]; push!(pts, pts[1])
        lines!(ax, [p.x for p in pts], [p.y for p in pts], color=:orange, linewidth=1.5)
    end
    
    n_faces = length(field.topology.faces)
    sample_rate = max(1, n_faces ÷ 500)
    
    xs, ys, us, vs = Float64[], Float64[], Float64[], Float64[]
    for i in 1:n_faces
        if i % sample_rate != 0; continue; end
        
        re = field.topology.face_ref_edges[i]
        p1 = field.topology.vertices[re[1]]; p2 = field.topology.vertices[re[2]]
        ref_ang = atan(p2.y - p1.y, p2.x - p1.x)
        ang = ref_ang + field.theta[i]
        
        tri = field.topology.faces[i]
        c = (field.topology.vertices[tri[1]].x + field.topology.vertices[tri[2]].x + field.topology.vertices[tri[3]].x)/3
        cy = (field.topology.vertices[tri[1]].y + field.topology.vertices[tri[2]].y + field.topology.vertices[tri[3]].y)/3
        
        push!(xs, c); push!(ys, cy)
        push!(us, cos(ang)*0.03); push!(vs, sin(ang)*0.03)
    end
    
    arrows!(ax, xs, ys, us, vs, color=:blue)
    
    save(filename, fig)
    return nothing
end

function plot_results(field, filename; verbose=false, cut_edges=nothing, show_period_jumps=false)
    fig = Figure(size=(1000, 800))
    ax = Axis(fig[1, 1], aspect=DataAspect(), title="MIQ Solution (Constrained Faces & Singularities)")
    
    for tri in field.topology.faces
        pts = [field.topology.vertices[i] for i in tri]; push!(pts, pts[1])
        lines!(ax, [p.x for p in pts], [p.y for p in pts], color=:gray90, linewidth=0.5)
    end
    
    n_faces = length(field.topology.faces)
    sample_rate = if n_faces < 1000; 1
    elseif n_faces < 10000; 2
    elseif n_faces < 50000; 5
    elseif n_faces < 200000; 10
    else; 500; end
    
    if verbose && sample_rate > 1
        println("Plotting every $(sample_rate)th cross field (mesh has $n_faces faces)")
    end
    
    if show_period_jumps
        p_values = [abs(p) for p in values(field.period_jumps)]
        if !isempty(p_values)
            max_p = maximum(p_values)
            
            for ((i, j), p) in field.period_jumps
                tri_i = field.topology.faces[i]
                cx_i = (field.topology.vertices[tri_i[1]].x + field.topology.vertices[tri_i[2]].x + field.topology.vertices[tri_i[3]].x)/3
                cy_i = (field.topology.vertices[tri_i[1]].y + field.topology.vertices[tri_i[2]].y + field.topology.vertices[tri_i[3]].y)/3
                
                tri_j = field.topology.faces[j]
                cx_j = (field.topology.vertices[tri_j[1]].x + field.topology.vertices[tri_j[2]].x + field.topology.vertices[tri_j[3]].x)/3
                cy_j = (field.topology.vertices[tri_j[1]].y + field.topology.vertices[tri_j[2]].y + field.topology.vertices[tri_j[3]].y)/3
                
                color_val = abs(p) / max(max_p, 1)
                edge_color = if p == 0; (:green, 0.6)
                elseif abs(p) == 1; (:yellow, 0.8)
                else; (:red, min(0.9, 0.5 + 0.4 * color_val)); end
                
                linewidth = abs(p) == 0 ? 1.5 : 2.5
                lines!(ax, [cx_i, cx_j], [cy_i, cy_j], color=edge_color, linewidth=linewidth)
            end
            
            if verbose
                println("Period jumps range: [$(minimum(values(field.period_jumps))), $(maximum(values(field.period_jumps)))]")
            end
        end
    else
        for edge in field.fixed_edges
            i, j = edge
            tri_i = field.topology.faces[i]
            cx_i = (field.topology.vertices[tri_i[1]].x + field.topology.vertices[tri_i[2]].x + field.topology.vertices[tri_i[3]].x)/3
            cy_i = (field.topology.vertices[tri_i[1]].y + field.topology.vertices[tri_i[2]].y + field.topology.vertices[tri_i[3]].y)/3
            
            tri_j = field.topology.faces[j]
            cx_j = (field.topology.vertices[tri_j[1]].x + field.topology.vertices[tri_j[2]].x + field.topology.vertices[tri_j[3]].x)/3
            cy_j = (field.topology.vertices[tri_j[1]].y + field.topology.vertices[tri_j[2]].y + field.topology.vertices[tri_j[3]].y)/3
            
            lines!(ax, [cx_i, cx_j], [cy_i, cy_j], color=:green, linewidth=2, alpha=0.6)
        end
    end
    
    for (face_idx, _) in field.constrained_faces
        tri = field.topology.faces[face_idx]
        pts = [field.topology.vertices[i] for i in tri]; push!(pts, pts[1])
        lines!(ax, [p.x for p in pts], [p.y for p in pts], color=:orange, linewidth=2)
        poly!(ax, [Point2f(p.x, p.y) for p in pts[1:3]], color=(:orange, 0.2))
    end
    
    xs_main, ys_main, us_main, vs_main = Float64[], Float64[], Float64[], Float64[]
    xs_sec, ys_sec, us_sec, vs_sec = Float64[], Float64[], Float64[], Float64[]
    
    for i in 1:length(field.topology.faces)
        if i % sample_rate != 0
            continue
        end
        
        re = field.topology.face_ref_edges[i]
        p1 = field.topology.vertices[re[1]]; p2 = field.topology.vertices[re[2]]
        ref_ang = atan(p2.y - p1.y, p2.x - p1.x)
        ang = ref_ang + field.theta[i]
        
        tri = field.topology.faces[i]
        c = (field.topology.vertices[tri[1]].x + field.topology.vertices[tri[2]].x + field.topology.vertices[tri[3]].x)/3
        cy = (field.topology.vertices[tri[1]].y + field.topology.vertices[tri[2]].y + field.topology.vertices[tri[3]].y)/3
        
        push!(xs_main, c); push!(ys_main, cy); push!(us_main, cos(ang)*0.03); push!(vs_main, sin(ang)*0.03)
        
        for k in 1:3
            a = ang + k*π/2
            push!(xs_sec, c); push!(ys_sec, cy); push!(us_sec, cos(a)*0.03); push!(vs_sec, sin(a)*0.03)
        end
    end
    arrows!(ax, xs_sec, ys_sec, us_sec, vs_sec, color=:blue)
    arrows!(ax, xs_main, ys_main, us_main, vs_main, color=:red)
    
    sings = compute_singularities(field; verbose=verbose)
    pos = Point2f[]; neg = Point2f[]
    for (v, I) in sings
        p = field.topology.vertices[v]
        pt = Point2f(p.x, p.y)
        if I > 0; push!(pos, pt); else; push!(neg, pt); end
    end
    scatter!(ax, pos, color=:red, markersize=15, strokecolor=:black, strokewidth=1, label="+ Index")
    scatter!(ax, neg, color=:cyan, markersize=15, strokecolor=:black, strokewidth=1, label="- Index")
    
    if cut_edges !== nothing
        cx, cy = Float64[], Float64[]
        for (u, v) in cut_edges
            p1, p2 = field.topology.vertices[u], field.topology.vertices[v]
            push!(cx, p1.x, p2.x, NaN)
            push!(cy, p1.y, p2.y, NaN)
        end
        lines!(ax, cx, cy, color=:red, linewidth=5.0)
        lines!(ax, [NaN], [NaN], color=:red, linewidth=5.0, label="Mesh Cut Edges")
    end

    if !isempty(field.constrained_faces)
        lines!(ax, [NaN], [NaN], color=:orange, linewidth=2, label="Constrained Faces")
    end
    
    if show_period_jumps
        lines!(ax, [NaN], [NaN], color=(:green, 0.6), linewidth=1.5, label="p = 0")
        lines!(ax, [NaN], [NaN], color=(:yellow, 0.8), linewidth=2.5, label="|p| = 1")
        lines!(ax, [NaN], [NaN], color=(:red, 0.9), linewidth=2.5, label="|p| ≥ 2")
    else
        if !isempty(field.fixed_edges)
            lines!(ax, [NaN], [NaN], color=:green, linewidth=2, alpha=0.6, label="Spanning Tree")
        end
    end
    
    if !isempty(pos) || !isempty(neg) || !isempty(field.constrained_faces) || !isempty(field.fixed_edges) || show_period_jumps
        axislegend(ax)
    end
    
    save(filename, fig)
    if verbose
        println("Saved visualization to $filename")
        println("Found $(length(sings)) singularities.")
    end
end

function plot_smooth_global_field(field::CrossField, rotations::Vector{Int}, filename::String; cut_edges=nothing, verbose=false)
    fig = Figure(size=(1000, 800))
    ax = Axis(fig[1, 1], aspect=DataAspect(), title="Globally Smoothed 'u' Direction (Seams at Cuts)")
    
    for tri in field.topology.faces
        pts = [field.topology.vertices[i] for i in tri]; push!(pts, pts[1])
        lines!(ax, [p.x for p in pts], [p.y for p in pts], color=(:gray, 0.3), linewidth=0.5)
    end

    if cut_edges !== nothing
        cx, cy = Float64[], Float64[]
        for (u, v) in cut_edges
            p1, p2 = field.topology.vertices[u], field.topology.vertices[v]
            push!(cx, p1.x, p2.x, NaN)
            push!(cy, p1.y, p2.y, NaN)
        end
        lines!(ax, cx, cy, color=:red, linewidth=4.0, label="Cut Graph")
    end

    xs, ys, us, vs = Float64[], Float64[], Float64[], Float64[]
    
    n_faces = length(field.topology.faces)
    sample_rate = n_faces > 10000 ? (n_faces ÷ 2000) : 1
    
    for i in 1:n_faces
        if i % sample_rate != 0; continue; end

        re = field.topology.face_ref_edges[i]
        p1 = field.topology.vertices[re[1]]
        p2 = field.topology.vertices[re[2]]
        ref_ang = atan(p2.y - p1.y, p2.x - p1.x)
        
        local_theta = field.theta[i]
        
        r = rotations[i]
        
        global_angle = ref_ang + local_theta + (r * π / 2.0)
        
        tri = field.topology.faces[i]
        c = (field.topology.vertices[tri[1]].x + field.topology.vertices[tri[2]].x + field.topology.vertices[tri[3]].x)/3
        cy = (field.topology.vertices[tri[1]].y + field.topology.vertices[tri[2]].y + field.topology.vertices[tri[3]].y)/3
        
        push!(xs, c); push!(ys, cy)
        push!(us, cos(global_angle) * 0.04)
        push!(vs, sin(global_angle) * 0.04)
    end

    arrows!(ax, xs, ys, us, vs, color=:blue, label="Global 'u' Field")
    
    sings = compute_singularities(field)
    pos = Point2f[]; neg = Point2f[]
    for (v, I) in sings
        p = field.topology.vertices[v]
        if I > 0; push!(pos, Point2f(p.x, p.y)); else; push!(neg, Point2f(p.x, p.y)); end
    end
    if !isempty(pos); scatter!(ax, pos, color=:red, markersize=12, label="Singularity (+)"); end
    if !isempty(neg); scatter!(ax, neg, color=:cyan, markersize=12, label="Singularity (-)"); end

    axislegend(ax)
    save(filename, fig)
    if verbose; println("Saved smooth field plot to $filename"); end
end

function plot_quad_mesh(topo::MeshTopology, u_coords::Vector{Float64}, v_coords::Vector{Float64}, 
                        quads::Vector{Tuple{Int,Int,Int,Int}}, filename::String; verbose=false)
    println("\n--- Plotting Quad Mesh ---")
    
    fig = Figure(size=(1600, 500))
    
    ax1 = Axis(fig[1, 1], aspect=DataAspect(), title="UV Parameter Space")
    
    for quad in quads
        v1, v2, v3, v4 = quad
        us = [u_coords[v1], u_coords[v2], u_coords[v3], u_coords[v4], u_coords[v1]]
        vs = [v_coords[v1], v_coords[v2], v_coords[v3], v_coords[v4], v_coords[v1]]
        lines!(ax1, us, vs, color=:blue, linewidth=1.5)
    end
    
    scatter!(ax1, u_coords, v_coords, color=:red, markersize=4)
    
    ax2 = Axis(fig[1, 2], aspect=DataAspect(), title="Parameterization on Mesh")
    
    for tri in topo.faces
        pts = [topo.vertices[i] for i in tri]
        xs = [pts[1].x, pts[2].x, pts[3].x, pts[1].x]
        ys = [pts[1].y, pts[2].y, pts[3].y, pts[1].y]
        lines!(ax2, xs, ys, color=(:gray, 0.15), linewidth=0.3)
    end
    
    xs = [v.x for v in topo.vertices]
    ys = [v.y for v in topo.vertices]
    scatter!(ax2, xs, ys, color=u_coords, colormap=:viridis, markersize=4)
    
    ax3 = Axis(fig[1, 3], aspect=DataAspect(), title="Quad Mesh ($(length(quads)) quads)")
    
    for tri in topo.faces
        pts = [topo.vertices[i] for i in tri]
        xs = [pts[1].x, pts[2].x, pts[3].x, pts[1].x]
        ys = [pts[1].y, pts[2].y, pts[3].y, pts[1].y]
        lines!(ax3, xs, ys, color=(:gray, 0.08), linewidth=0.2)
    end
    
    quad_qualities = Float64[]
    
    for quad in quads
        v1, v2, v3, v4 = quad
        p1, p2, p3, p4 = topo.vertices[v1], topo.vertices[v2], topo.vertices[v3], topo.vertices[v4]
        
        e1 = sqrt((p2.x-p1.x)^2 + (p2.y-p1.y)^2)
        e2 = sqrt((p3.x-p2.x)^2 + (p3.y-p2.y)^2)
        e3 = sqrt((p4.x-p3.x)^2 + (p4.y-p3.y)^2)
        e4 = sqrt((p1.x-p4.x)^2 + (p1.y-p4.y)^2)
        
        min_edge = min(e1, e2, e3, e4)
        max_edge = max(e1, e2, e3, e4)
        quality = max_edge > 1e-10 ? min_edge / max_edge : 0.0
        
        push!(quad_qualities, quality)
    end
    
    for (q_idx, quad) in enumerate(quads)
        v1, v2, v3, v4 = quad
        pts = [topo.vertices[v1], topo.vertices[v2], topo.vertices[v3], topo.vertices[v4]]
        xs = [pts[1].x, pts[2].x, pts[3].x, pts[4].x, pts[1].x]
        ys = [pts[1].y, pts[2].y, pts[3].y, pts[4].y, pts[1].y]
        
        quality = quad_qualities[q_idx]
        color = quality > 0.5 ? :green : (quality > 0.3 ? :orange : :red)
        
        lines!(ax3, xs, ys, color=color, linewidth=2)
    end
    
    Colorbar(fig[1, 4], limits=(minimum(u_coords), maximum(u_coords)), 
             colormap=:viridis, label="u-coordinate")
    
    save(filename, fig)
    
    if verbose
        println("Saved quad mesh visualization to $filename")
        if !isempty(quad_qualities)
            @printf("Quad quality - min: %.3f, avg: %.3f, max: %.3f\n", 
                    minimum(quad_qualities), sum(quad_qualities)/length(quad_qualities), maximum(quad_qualities))
        end
    end
end

end