using GLMakie
using LinearAlgebra
using FrameFieldFull
using FrameFieldFull.Types
using FrameFieldFull.Topology
using FrameFieldFull.Analysis
import GeometryBasics



# --- Main Visualization Function ---
function visualize_field(field::CrossField; title="Frame Field Debug", final_cuts=nothing, save_path=nothing, two_d=false)
    topo = field.topology
    
    # 1. Setup Scene
    fig = Figure(size=(1000, 800))
    
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
    wireframe!(ax, geo_mesh, color=(:black, 0.1), linewidth=0.5)
    
    # 3. Plot Field Vectors
    main_direction = Vec3f[]
    main_centers = Point3f[]
    other_directions = Vec3f[]
    other_centers = Point3f[]

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
        c = get_face_center(topo, i)
        v = get_global_vector(topo, i, field.theta[i])
        
        push!(main_centers, c)
        push!(main_direction, Vec3f(v...))
        
        # Compute and store other directions (90, 180, 270 degrees rotations)
        for k in 1:3
            theta_rot = field.theta[i] + k * (π/2)
            v_rot = get_global_vector(topo, i, theta_rot)
            push!(other_directions, Vec3f(v_rot...))
            push!(other_centers, c)
        end
    end
    
    arrows2d!(ax, main_centers, main_direction, lengthscale=vec_scale, 
            shaftcolor=:red, tipcolor=:red, label="Major Direction")

    arrows2d!(ax, other_centers, other_directions, lengthscale=vec_scale, 
            shaftcolor=:black, tipcolor=:black, label="Other Directions")
    
            
    # 4. Plot Singularities
    sings = compute_singularities(field; verbose=false)
    pos_sings = Point3f[]
    neg_sings = Point3f[]
    
    for (v_idx, idx_val) in sings
        p = pts[v_idx]
        if idx_val > 0.1
            push!(pos_sings, p)
        elseif idx_val < -0.1
            push!(neg_sings, p)
        end
    end
    
    if !isempty(pos_sings)
        meshscatter!(ax, pos_sings, color=:red, markersize=avg_len*0.3, label="Index > 0")
    end
    if !isempty(neg_sings)
        meshscatter!(ax, neg_sings, color=:blue, markersize=avg_len*0.3, label="Index < 0")
    end
    
    # 5. Plot Period Jumps (Cuts)
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

    # 6. Plot Constraints
    constraint_centers = Point3f[]
    constraint_vecs = Vec3f[]
    constraint_edges = Point3f[]

    for (f_idx, angle) in field.constrained_faces
        # Find which edge is the boundary
        tri = topo.faces[f_idx]
        local_edges = [(tri[1], tri[2]), (tri[2], tri[3]), (tri[3], tri[1])]
        
        found_boundary = false
        for edge in local_edges
            key = minmax(edge[1], edge[2])
            has_neighbor = false
            for (_, n_key) in topo.dual_adj[f_idx]
                if n_key == key; has_neighbor = true; break; end
            end
            
            if !has_neighbor
                p1 = pts[edge[1]]
                p2 = pts[edge[2]]
                push!(constraint_edges, p1)
                push!(constraint_edges, p2)
                
                mid = (p1 + p2) / 2.0
                push!(constraint_centers, mid)
                
                v_global = get_global_vector(topo, f_idx, angle)
                push!(constraint_vecs, Vec3f(v_global...))
                
                found_boundary = true
            end
        end
        
        if !found_boundary
            push!(constraint_centers, get_face_center(topo, f_idx))
            v_global = get_global_vector(topo, f_idx, angle)
            push!(constraint_vecs, Vec3f(v_global...))
        end
    end

    if !isempty(constraint_edges)
        linesegments!(ax, constraint_edges, color=:cyan, linewidth=5.0, label="Constrained Edges")
        arrows2d!(ax, constraint_centers, constraint_vecs, lengthscale=vec_scale*1.2, 
                shaftcolor=:cyan, tipcolor=:cyan, label="Target Dir")
    end
    
    # 7. Plot Final Cuts
    if !isnothing(final_cuts)
        cut_segments = Point3f[]
        for (v1, v2) in final_cuts
             push!(cut_segments, pts[v1])
             push!(cut_segments, pts[v2])
        end
        
        if length(cut_segments) != 2 * length(final_cuts)
             @warn "Mismatch in cut segments generation"
        end
        
        if !isempty(cut_segments)
            linesegments!(ax, cut_segments, color=:magenta, linewidth=6.0, label="Mesh Cuts")
        end
    end
    
    Legend(fig[1, 2], ax)
    
    if !isnothing(save_path)
        save(save_path, fig)
        println("Figure saved to: $save_path")
    end
    
    display(fig)
end



if abspath(PROGRAM_FILE) == @__FILE__
    println("Load this file and call visualize_field(field) after running your pipeline.")
end