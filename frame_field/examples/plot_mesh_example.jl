# Runner: load mesh and save a PNG visualization
include("../src/meshio.jl")
include("../src/plotters.jl")

case_folder_name = "triangulations"

case_names = [
    "simple-square",
    "300_polygon_sphere_100mm",
]

for case_name in case_names
    mesh_path = joinpath(@__DIR__, "..", "..", "triangulations", "$case_name.msh")

    println("Loading mesh: ", mesh_path)
    mesh = load_triangulation(mesh_path)

    outpath = joinpath(@__DIR__, "..", "output", "$(case_folder_name)")
    mkpath(outpath)
    outpath = joinpath(outpath, "triangulation_$(case_name).png")
    println("Saving visualization to: ", outpath)
    fig = plot_triangulation(mesh; savepath=outpath, show_edges=true)
    println("Done.")
end