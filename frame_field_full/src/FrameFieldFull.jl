module FrameFieldFull

# Dependencies
using LinearAlgebra
using SparseArrays
using Printf
using CairoMakie
using DataStructures

# Include source files
include("types.jl")
include("mesh_io.jl")
include("topology.jl")
include("analysis.jl")
include("plotting.jl")
include("solver.jl")
include("cutting.jl")
include("parameterization.jl")
include("extraction.jl")

# Import from submodules
using .Types
using .MeshIO
using .Topology
using .Solver
using .Plotting
using .Analysis
using .Cutting
using .Parameterization
using .Extraction

# Export user-facing functions
export read_mesh, 
       build_topology, 
       compute_boundary_constraints,
       initialize_field,
       solve_greedy!,
       optimize_singularities!,
       compute_singularities,
       compute_cut_graph,
       propagate_orientations,
       compute_parameterization_least_squares,
       extract_quad_mesh,
       plot_results,
       plot_smooth_global_field,
       plot_quad_mesh,
       plot_frame,
       QuadMesh, 
       extract_quad_mesh

end # module