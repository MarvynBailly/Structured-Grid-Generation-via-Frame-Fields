module Plotting

using ..Types
using ..Analysis
using CairoMakie
using GLMakie
using GeometryBasics
using Printf
using LinearAlgebra

export plot_frame, plot_results, plot_smooth_global_field, plot_quad_mesh, plot_global_rotations, plot_extracted_quads


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



function plot_global_rotations(field::CrossField, rotations::Vector{Int}; save_path::String=nothing, verbose=false, title="Global Rotations", two_d=false, final_cuts=nothing)
    topo = field.topology

    fig = Figure(size=(1200, 1000))

    # Configure Axis based on 2D/3D preference
    ax_settings = Dict{Symbol, Any}(
        :title => title, 
        :aspect => :data
    )
    
    if two_d
        # Top-down view settings
        ax_settings[:elevation] = π/2       # Look from Top (Z-axis)
        ax_settings[:azimuth] = -π/2        # Rotate so X is Right, Y is Up
        ax_settings[:perspectiveness] = 0.0 # Orthographic (No depth distortion)
    end
    
    ax = Axis3(fig[1, 1]; ax_settings...)

    # 2. Plot Mesh (Wireframe)
    pts = [Point3f(v.x, v.y, v.z) for v in topo.vertices]
    tris = [GLMakie.TriangleFace(f[1], f[2], f[3]) for f in topo.faces]
    
    geo_mesh = GeometryBasics.normal_mesh(pts, tris)
    
    mesh!(ax, geo_mesh, color=(:gray, 0.3), transparency=true, shading=NoShading)
    wireframe!(ax, geo_mesh, color=(:black), linewidth=0.5)

    # plot Global U and V Vectors
    us, vs = Vec3f[], Vec3f[]
    centers = Point3f[]

    # Compute avg edge length for scaling
    total_len = 0.0
    count = 0
    for i in 1:length(topo.faces)
        tri = topo.faces[i]
        p1, p2 = pts[tri[1]], pts[tri[2]]
        total_len += norm(p1 - p2)
        count += 1
    end
    avg_len = count > 0 ? total_len / count : 1.0
    vec_scale = avg_len * 0.2

    for i in 1:length(topo.faces)
        r = rotations[i]
        c = get_face_center(topo, i)
        u = get_global_vector(topo, i, field.theta[i] + r * π/2)
        v = get_global_vector(topo, i, field.theta[i] + r * π/2 + π/2)
        
        push!(centers, c)
        push!(us, Vec3f(u...))
        push!(vs, Vec3f(v...))
    end
    
    if verbose
        # Print rotation statistics
        rotation_counts = Dict(0=>0, 1=>0, 2=>0, 3=>0, -1=>0)
        for r in rotations
            rotation_counts[r] = get(rotation_counts, r, 0) + 1
        end
        println("Rotation distribution:")
        for (r, count) in sort(collect(rotation_counts))
            if count > 0
                println("  r=$r: $count faces")
            end
        end
    end

    arrows2d!(ax, centers, us, lengthscale=vec_scale, color=:red, label="Global u")
    arrows2d!(ax, centers, vs, lengthscale=vec_scale, color=:blue, label="Global v")

    jump_segments = Point3f[]
    for (edge_key, p_val) in field.period_jumps
        if p_val != 0
            f_i, f_j = edge_key
            for (neighbor, primal_key) in topo.dual_adj[f_i]
                if neighbor == f_j
                    v1, v2 = primal_key
                    push!(jump_segments, pts[v1])
                    push!(jump_segments, pts[v2])
                    break
                end
            end
        end
    end
    


    if !isempty(jump_segments)
        linesegments!(ax, jump_segments, color=:orange, linewidth=4.0, label="Period Jumps")
    end


    if final_cuts !== nothing && !isempty(final_cuts)
        cut_segments = Point3f[]
        for (u, v) in final_cuts
            p1, p2 = topo.vertices[u], topo.vertices[v]
            push!(cut_segments, Point3f(p1.x, p1.y, p1.z))
            push!(cut_segments, Point3f(p2.x, p2.y, p2.z))
        end
        linesegments!(ax, cut_segments, color=:red, linewidth=5.0, label="Final Cut Edges")
    end



    Legend(fig[1, 2], ax)

    if !isnothing(save_path)
        save(save_path, fig)
        println("Figure saved to: $save_path")
    end
    
    display(fig)
end




function plot_smooth_global_field(field::CrossField, rotations::Vector{Int}, filename::String; 
                                  cut_edges=Set{Tuple{Int,Int}}(), root_face::Int=0, verbose=false,)
    fig = Figure(size=(1200, 1000))
    ax = Axis(fig[1, 1], aspect=DataAspect(), title="Global Parametrization Field (Smoothed)")
    
    topo = field.topology
    
    # --- 1. Plot Mesh Wireframe ---
    for tri in topo.faces
        pts = [topo.vertices[i] for i in tri]; push!(pts, pts[1])
        lines!(ax, [p.x for p in pts], [p.y for p in pts], color=(:gray, 0.2), linewidth=0.5)
    end

    # --- 2. Highlight Root Face ---
    if root_face > 0 && root_face <= length(topo.faces)
        tri = topo.faces[root_face]
        pts = [topo.vertices[i] for i in tri]
        # Draw filled polygon
        poly!(ax, [Point2f(p.x, p.y) for p in pts], color=(:yellow, 0.5), strokecolor=:black, strokewidth=2, label="Root Face")
    end

    # --- 3. Plot Global U and V Vectors ---
    xs, ys = Float64[], Float64[]
    us_u, vs_u = Float64[], Float64[] # U vector components
    us_v, vs_v = Float64[], Float64[] # V vector components
    
    n_faces = length(topo.faces)
    # Adjust sample rate for density
    sample_rate = n_faces > 2000 ? (n_faces ÷ 1000) : 1
    scale = 0.04 # Adjust based on mesh scale

    for i in 1:n_faces
        if i % sample_rate != 0; continue; end
        
        # Determine Global Rotation Angle
        # Global = Reference + Local_Theta + (r * pi/2)
        re = topo.face_ref_edges[i]
        p1 = topo.vertices[re[1]]
        p2 = topo.vertices[re[2]]
        ref_ang = atan(p2.y - p1.y, p2.x - p1.x)
        
        local_theta = field.theta[i]
        r = rotations[i]
        
        # The Smooth Global Angle
        angle_u = ref_ang + local_theta + (r * π / 2.0)
        angle_v = angle_u + (π / 2.0)
        
        # Face Center
        tri = topo.faces[i]
        c = (topo.vertices[tri[1]].x + topo.vertices[tri[2]].x + topo.vertices[tri[3]].x)/3
        cy = (topo.vertices[tri[1]].y + topo.vertices[tri[2]].y + topo.vertices[tri[3]].y)/3
        
        push!(xs, c); push!(ys, cy)
        
        # U Vector (Red)
        push!(us_u, cos(angle_u) * scale)
        push!(vs_u, sin(angle_u) * scale)
        
        # V Vector (Green)
        push!(us_v, cos(angle_v) * scale)
        push!(vs_v, sin(angle_v) * scale)
    end

    arrows2d!(ax, xs, ys, us_u, vs_u, color=:red, label="Global u")
    arrows2d!(ax, xs, ys, us_v, vs_v, color=:blue, label="Global v")

    # --- 4. Highlight Edges with Non-Zero Transition Rotation ---
    # We check every edge. If the global orientation implies a rotation mismatch, we plot it.
    # Ideally, this only happens on Cut Edges.
    
    mismatch_segments = Point2f[]
    mismatch_colors = Symbol[]
    
    for i in 1:n_faces
        for (j, edge_key) in topo.dual_adj[i]
            if i < j # Process each edge once
                # Get Period Jump p_ij (defined for i->j direction if i<j)
                p = get(field.period_jumps, edge_key, 0)
                
                # Global Mismatch = (r_j - r_i) + p
                # If smooth, this sum should be 0 (mod 4).
                # The "Transition Rotation" is exactly this residue.
                
                # Note: Rotations are defined mod 4.
                # transition = (r_j - r_i + p) % 4
                trans = mod(rotations[j] - rotations[i] + p, 4)
                
                # If trans == 0, the chart merge is identity (Smooth).
                # If trans != 0, there is a rotation (Seam).
                if trans != 0
                    v1 = topo.vertices[edge_key[1]]
                    v2 = topo.vertices[edge_key[2]]
                    push!(mismatch_segments, Point2f(v1.x, v1.y))
                    push!(mismatch_segments, Point2f(v2.x, v2.y))
                    push!(mismatch_segments, Point2f(NaN, NaN))
                    
                    # Color code by rotation amount?
                    # 1 or 3 (-1) are 90 degree turns. 2 is 180 degree turn.
                    if trans == 2
                        push!(mismatch_colors, :cyan) # 180 flip
                    else
                        push!(mismatch_colors, :magenta) # 90 turn
                    end
                end
            end
        end
    end

    if !isempty(mismatch_segments)
        lines!(ax, mismatch_segments, color=:magenta, linewidth=4.0, label="Rotation Jump (r ≠ 0)")
    end
    
    # --- 5. Plot Cut Graph (for Reference) ---
    if !isempty(cut_edges)
        cut_segs = Point2f[]
        for (u, v) in cut_edges
            p1 = topo.vertices[u]
            p2 = topo.vertices[v]
            push!(cut_segs, Point2f(p1.x, p1.y))
            push!(cut_segs, Point2f(p2.x, p2.y))
            push!(cut_segs, Point2f(NaN, NaN))
        end
        lines!(ax, cut_segs, color=:black, linestyle=:dash, linewidth=2.0, label="Cut Graph")
    end

    # --- 6. Plot Singularities ---
    sings = compute_singularities(field)
    pos = Point2f[]
    neg = Point2f[]
    for (v, I) in sings
        p = topo.vertices[v]
        if I > 0
            push!(pos, Point2f(p.x, p.y))
        else
            push!(neg, Point2f(p.x, p.y))
        end
    end
    if !isempty(pos); scatter!(ax, pos, color=:red, markersize=12, label="Singularity (+)"); end
    if !isempty(neg); scatter!(ax, neg, color=:blue, markersize=12, label="Singularity (-)"); end

    axislegend(ax, position=:rt)
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


"""
    plot_extracted_quads(vertices, quads, filename; topo=nothing)

Visualizes the final extracted quad mesh in 3D.
- vertices: Vector of Point3D (the nodes of the quad mesh)
- quads: Vector of 4-tuples (indices into vertices)
- topo: Optional background triangle mesh for reference
"""
function plot_extracted_quads(vertices::Vector{Point3D}, quads::Vector{Tuple{Int,Int,Int,Int}}, filename::String; 
                              topo::Union{MeshTopology, Nothing}=nothing, verbose=false)
    if verbose; println("\n--- Plotting Extracted Quad Mesh ---"); end
    
    fig = Figure(size=(1000, 800))
    ax = Axis3(fig[1, 1], aspect=:data, title="Extracted Quad Mesh")
    
    # 1. Plot Background (Original Triangle Mesh) - Ghosted
    if !isnothing(topo)
        t_pts = [Point3f(v.x, v.y, v.z) for v in topo.vertices]
        t_faces = [GLMakie.TriangleFace(f[1], f[2], f[3]) for f in topo.faces]
        t_mesh = GeometryBasics.normal_mesh(t_pts, t_faces)
        # Very transparent gray
        mesh!(ax, t_mesh, color=(:gray, 0.1), transparency=true, shading=NoShading)
    end
    
    # 2. Prepare Quad Mesh Geometry
    # Convert vertices to Float32 points for Makie
    q_pts = [Point3f(v.x, v.y, v.z) for v in vertices]
    
    # Makie requires triangles for the 'mesh' command.
    # We split each quad (v1,v2,v3,v4) into two triangles: (v1,v2,v3) and (v1,v3,v4)
    tri_faces = GLMakie.TriangleFace{Int}[]
    wire_segments = Point3f[]
    
    for (v1, v2, v3, v4) in quads
        # Surface Triangles
        push!(tri_faces, GLMakie.TriangleFace(v1, v2, v3))
        push!(tri_faces, GLMakie.TriangleFace(v1, v3, v4))
        
        # Wireframe Edges
        p1, p2, p3, p4 = q_pts[v1], q_pts[v2], q_pts[v3], q_pts[v4]
        push!(wire_segments, p1, p2)
        push!(wire_segments, p2, p3)
        push!(wire_segments, p3, p4)
        push!(wire_segments, p4, p1)
    end
    
    # 3. Plot Surface (Cyan)
    if !isempty(tri_faces)
        q_mesh = GeometryBasics.normal_mesh(q_pts, tri_faces)
        mesh!(ax, q_mesh, color=(:cyan, 0.8), shading=NoShading)
    end
    
    # 4. Plot Wireframe (Black)
    if !isempty(wire_segments)
        linesegments!(ax, wire_segments, color=:black, linewidth=1.5)
    end
    
    save(filename, fig)
    if verbose; println("Saved to $filename"); end
end

end