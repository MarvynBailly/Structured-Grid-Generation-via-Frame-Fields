# ---------------------------------------------------------
# Julia Script to Generate a 2-Element Airfoil Gmsh File
# ---------------------------------------------------------

function naca4_points(t_pct, chord, n_points)
    beta = range(0, π, length=n_points)
    x = (1.0 .- cos.(beta)) ./ 2.0 .* chord
    yt = 5 * t_pct * chord * (
        0.2969 .* sqrt.(x/chord) .- 0.1260 .* (x/chord) .-
        0.3516 .* (x/chord).^2 .+ 0.2843 .* (x/chord).^3 .-
        0.1015 .* (x/chord).^4
    )
    # Upper surface
    x_up, y_up = x, yt
    # Lower surface (reverse to create loop)
    x_lo, y_lo = x[end-1:-1:1], -yt[end-1:-1:1]
    
    return vcat(x_up, x_lo), vcat(y_up, y_lo)
end

function transform(x, y, dx, dy, ang_deg)
    rad = deg2rad(ang_deg)
    xn = @. x * cos(rad) - y * sin(rad) + dx
    yn = @. x * sin(rad) + y * cos(rad) + dy
    return xn, yn
end

# --- Configuration ---
filename = "two_element.geo"
farfield_radius = 2.0
mesh_size_airfoil = 0.02  # Fine resolution on wall
mesh_size_farfield = 1.0  # Coarse resolution far away

# Generate Geometry Arrays
x_m, y_m = naca4_points(0.12, 1.0, 60) # Main
x_f, y_f = naca4_points(0.12, 0.4, 40) # Flap
x_f, y_f = transform(x_f, y_f, 1.05, -0.05, -20.0) # Move Flap

open(filename, "w") do io
    println(io, "// Gmsh Geometry for SU2 - Two Element Airfoil")
    
    # 1. Write Points for Main Element
    main_pt_ids = Int[]
    pid = 1
    for i in 1:length(x_m)
        println(io, "Point($pid) = {$(x_m[i]), $(y_m[i]), 0, $mesh_size_airfoil};")
        push!(main_pt_ids, pid)
        pid += 1
    end
    # Close the loop explicitly in geometry if needed, or rely on Spline
    # Note: NACA trailing edge isn't perfectly sharp in coordinates, but we close it here
    println(io, "Spline(1) = {$(join(main_pt_ids, ", ")), $(main_pt_ids[1])};")
    println(io, "Curve Loop(1) = {1};")

    # 2. Write Points for Flap Element
    flap_start_pid = pid
    flap_pt_ids = Int[]
    for i in 1:length(x_f)
        println(io, "Point($pid) = {$(x_f[i]), $(y_f[i]), 0, $mesh_size_airfoil};")
        push!(flap_pt_ids, pid)
        pid += 1
    end
    println(io, "Spline(2) = {$(join(flap_pt_ids, ", ")), $(flap_pt_ids[1])};")
    println(io, "Curve Loop(2) = {2};")

    # 3. Create Circular Farfield
    # We define 4 points to make a circle
    c_pt = pid 
    R = farfield_radius
    println(io, "Point($pid) = {0, 0, 0, $mesh_size_farfield};") # Center
    println(io, "Point($(pid+1)) = {$R, 0, 0, $mesh_size_farfield};") 
    println(io, "Point($(pid+2)) = {0, $R, 0, $mesh_size_farfield};")
    println(io, "Point($(pid+3)) = {-$R, 0, 0, $mesh_size_farfield};")
    println(io, "Point($(pid+4)) = {0, -$R, 0, $mesh_size_farfield};")
    
    # Arcs
    println(io, "Circle(10) = {$(pid+1), $pid, $(pid+2)};")
    println(io, "Circle(11) = {$(pid+2), $pid, $(pid+3)};")
    println(io, "Circle(12) = {$(pid+3), $pid, $(pid+4)};")
    println(io, "Circle(13) = {$(pid+4), $pid, $(pid+1)};")
    println(io, "Curve Loop(3) = {10, 11, 12, 13};")

    # 4. Define Surface (Farfield - Main - Flap)
    println(io, "Plane Surface(1) = {3, 1, 2};")

    # 5. Define Physical Groups (Crucial for SU2 Markers)
    # "airfoil" will be the wall boundary
    # "farfield" will be the outer boundary
    # "fluid" will be the internal mesh area
    println(io, "Physical Curve(\"airfoil\") = {1, 2};")
    println(io, "Physical Curve(\"farfield\") = {10, 11, 12, 13};")
    println(io, "Physical Surface(\"fluid\") = {1};")
end

println("Successfully wrote $filename")
println("Next Step: Run 'gmsh $filename -2 -format su2 -o mesh.su2'")