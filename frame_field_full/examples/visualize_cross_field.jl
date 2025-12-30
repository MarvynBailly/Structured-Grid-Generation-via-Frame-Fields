using GLMakie
using LinearAlgebra
using FrameFieldFull
using FrameFieldFull.Types
using FrameFieldFull.Topology
using FrameFieldFull.Analysis
import GeometryBasics



# --- Main Visualization Function ---




if abspath(PROGRAM_FILE) == @__FILE__
    println("Load this file and call visualize_field(field) after running your pipeline.")
end