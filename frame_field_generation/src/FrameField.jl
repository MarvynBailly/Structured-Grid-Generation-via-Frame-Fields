"""
    FrameField

Module for generating smooth frame fields over triangular meshes for 
quadrilateral mesh generation.

Based on the methodology from:
Bommes et al. (2009) "Mixed-integer quadrangulation"

# Exports
- `generate_frame_field` - Main function to generate frame field
- `FrameFieldResult` - Result container
- `MeshTopology` - Mesh topology structure
- `compute_singularities` - Find singular vertices
- `plot_frame_field` - Visualization function
"""
module FrameField

using LinearAlgebra
using SparseArrays
using DataStructures
using Statistics
using GeometryBasics

# Include submodules
include("mesh_topology.jl")
include("dijkstra_forest.jl")
include("energy.jl")
include("frame_field_solver.jl")
include("visualization.jl")

# Export main types and functions
export FrameFieldResult, MeshTopology, DijkstraForest
export generate_frame_field, compute_singularities
export build_mesh_topology, build_dijkstra_forest
export assemble_system_matrix, compute_smoothness_energy
export get_frame_directions
export plot_frame_field, plot_singularities, print_frame_field_summary

end # module
