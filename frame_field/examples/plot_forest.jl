include("../src/meshio.jl")
include("../src/plotters.jl")
include("../src/dijkstra_forest.jl")

case_folder_name = "forests"

case_name = "simple-square"
# case_name = "300_polygon_sphere_100mm"
# case_names = [
#     "simple-square",
#     "300_polygon_sphere_100mm",
# ]


# for case_name in case_names
    mesh_path = joinpath(@__DIR__, "..", "..", "triangulations", "$case_name.msh")
    println("Loading mesh: ", mesh_path)
    mesh = load_triangulation(mesh_path)

    # randomly pick n constrained faces
    n = 2
    total_faces = length(faces(mesh))
    constrained_faces = rand(1:total_faces, n)

    # compute the set of edges that can be set to zero
    potential_fixed_edges = compute_spanning_forest(mesh; constrained_faces=constrained_faces)

    # select a subset of these edges to be fixed
    fixed_edges_per_face = fix_suitable_edges(mesh, potential_fixed_edges)

    out_dir = joinpath(@__DIR__, "..", "output", "$(case_folder_name)")
    mkpath(out_dir)
    out_path = joinpath(out_dir, "triangulation_forest_$(case_name).png")

    fig = plot_forest(mesh; constrained_faces=constrained_faces, 
                      potential_fixed_edges=potential_fixed_edges, 
                      savepath=out_path)

    # plot the results showing only the fixed edges
    out_path = joinpath(out_dir, "triangulation_forest_fixed_$(case_name).png")

    fig = plot_forest(mesh; constrained_faces=constrained_faces, 
                      potential_fixed_edges=fixed_edges_per_face, 
                      savepath=out_path, display_background=false)


    println("Saved forest plot to ", out_path)
# end