"""
    visualization.jl

Visualization utilities for frame fields.

Note: Requires GLMakie or WGLMakie for plotting.
These are optional dependencies - install separately if needed:
    using Pkg; Pkg.add("GLMakie")
"""

using GeometryBasics

# Visualization functions are stubs - implement when plotting library is available

"""
    plot_frame_field(mesh::Mesh, result::FrameFieldResult; kwargs...)

Plot the frame field on the mesh.

Displays:
- Triangle mesh edges
- Frame directions at each face center
- Color-coded by angle value

# Optional Arguments
- `arrow_scale::Float64=0.1` - Scale factor for frame direction arrows
- `show_all_directions::Bool=false` - Show all 4 directions or just first
"""
function plot_frame_field(mesh::Mesh, result::FrameFieldResult; 
                         arrow_scale::Float64=0.1,
                         show_all_directions::Bool=false)
    @info "plot_frame_field: Visualization requires GLMakie or WGLMakie"
    @info "Install with: using Pkg; Pkg.add(\"GLMakie\")"
    
    # Stub implementation
    println("Frame Field Visualization:")
    println("  Faces: $(length(result.angles))")
    println("  Energy: $(result.energy)")
    println("  Singularities: $(length(result.singularities))")
    
    # TODO: Implement actual plotting when Makie is available
end

"""
    plot_singularities(mesh::Mesh, result::FrameFieldResult; kwargs...)

Plot singular vertices on the mesh.

Displays:
- Mesh with singularities highlighted
- Color-coded by singularity index
- Labels showing index values

# Optional Arguments
- `marker_size::Float64=20.0` - Size of singularity markers
- `show_labels::Bool=true` - Show index values as labels
"""
function plot_singularities(mesh::Mesh, result::FrameFieldResult;
                           marker_size::Float64=20.0,
                           show_labels::Bool=true)
    @info "plot_singularities: Visualization requires GLMakie or WGLMakie"
    
    println("\nSingularity Analysis:")
    println("  Total singular vertices: $(length(result.singularities))")
    
    if !isempty(result.singularities)
        println("\n  Vertex Index | Singularity Index | Quad Valence")
        println("  " * "-"^52)
        for (v_idx, index) in result.singularities
            # Convert index to expected quad valence
            valence = 4 + round(Int, 4 * index)
            println("  $(lpad(v_idx, 12)) | $(lpad(round(index, digits=4), 17)) | $(lpad(valence, 12))")
        end
    end
    
    # TODO: Implement actual plotting when Makie is available
end

"""
    plot_period_jumps(mesh::Mesh, topology::MeshTopology, result::FrameFieldResult)

Plot period jumps on mesh edges.

Shows edges colored by their period jump values.
"""
function plot_period_jumps(mesh::Mesh, topology::MeshTopology, result::FrameFieldResult)
    @info "plot_period_jumps: Visualization requires GLMakie or WGLMakie"
    
    # Count period jump distribution
    p_values = collect(values(result.period_jumps))
    unique_vals = unique(p_values)
    
    println("\nPeriod Jump Distribution:")
    for val in sort(unique_vals)
        count = sum(p_values .== val)
        println("  p = $val: $count edges")
    end
    
    # TODO: Implement actual plotting when Makie is available
end

"""
    export_frame_field_vtk(filename::String, mesh::Mesh, result::FrameFieldResult)

Export frame field to VTK format for visualization in ParaView.

# Arguments
- `filename::String` - Output filename (without extension)
- `mesh::Mesh` - Triangular mesh
- `result::FrameFieldResult` - Frame field solution
"""
function export_frame_field_vtk(filename::String, mesh::Mesh, result::FrameFieldResult)
    @info "export_frame_field_vtk: Requires WriteVTK package"
    @info "Install with: using Pkg; Pkg.add(\"WriteVTK\")"
    
    # TODO: Implement VTK export when WriteVTK is available
    println("Would export frame field to: $filename.vtu")
end

"""
    print_frame_field_summary(result::FrameFieldResult)

Print a text summary of the frame field solution.
"""
function print_frame_field_summary(result::FrameFieldResult)
    println("\n" * "="^60)
    println("Frame Field Summary")
    println("="^60)
    
    println("\nOptimization:")
    println("  Converged: $(result.converged)")
    println("  Iterations: $(result.iterations)")
    println("  Final Energy: $(result.energy)")
    
    println("\nFrame Field:")
    println("  Faces: $(length(result.angles))")
    println("  Angle range: [$(minimum(result.angles)), $(maximum(result.angles))]")
    
    println("\nPeriod Jumps:")
    println("  Edges: $(length(result.period_jumps))")
    p_vals = collect(values(result.period_jumps))
    if !isempty(p_vals)
        println("  Range: [$(minimum(p_vals)), $(maximum(p_vals))]")
        println("  Non-zero: $(sum(p_vals .!= 0))")
    end
    
    println("\nSingularities:")
    println("  Singular vertices: $(length(result.singularities))")
    if !isempty(result.singularities)
        indices = [idx for (_, idx) in result.singularities]
        println("  Index range: [$(minimum(indices)), $(maximum(indices))]")
    end
    
    println("\n" * "="^60 * "\n")
end
