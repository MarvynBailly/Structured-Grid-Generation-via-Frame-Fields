using LinearAlgebra
using Printf

# --- Configuration ---
INPUT_FILE = "triangulations/2DMEA_RevM_wake_L1.b8.ugrid"  # Your large mesh
OUTPUT_FILE = "crmhl_test.su2"
RADIUS = 1.0                      # Distance from (0,0) to keep
CENTER = (0.0, 0.0)               # Approximate center of the airfoil

# --- Data Structures ---
struct Point3D
    x::Float64; y::Float64; z::Float64
end

# --- UGRID Binary Reader (Fortran Unformatted / C3D Binary Format) ---
function read_ugrid_binary(filename::String)
    println("Reading binary UGRID file: $filename...")
    
    f = open(filename, "r")
    
    # Read header (7 integers: nnodes, ntria, nquad, ntetra, npyra, nprism, nhex)
    # Try big-endian first (common for UGRID binary files)
    counts_raw = read!(f, Vector{Int32}(undef, 7))
    counts = ntoh.(counts_raw)  # Network to host byte order (big-endian to native)
    
    n_nodes = Int(counts[1])
    n_tri = Int(counts[2])
    n_quad = Int(counts[3])
    
    # Sanity check - if values are crazy, file might be little-endian or different format
    if n_nodes < 0 || n_nodes > 1_000_000_000 || n_tri < 0 || n_tri > 1_000_000_000
        println("  Warning: Values seem wrong with big-endian, trying little-endian...")
        seekstart(f)
        counts_raw = read!(f, Vector{Int32}(undef, 7))
        counts = counts_raw  # Use as-is (little-endian)
        n_nodes = Int(counts[1])
        n_tri = Int(counts[2])
        n_quad = Int(counts[3])
    end
    
    println("  Original Nodes: $n_nodes")
    println("  Original Triangles: $n_tri")
    println("  Original Quads: $n_quad")
    println("  Original Faces: $(n_tri + n_quad)")
    
    # Read vertices (3 * n_nodes doubles)
    vertices = Point3D[]
    for i in 1:n_nodes
        x = ntoh(read(f, Float64))
        y = ntoh(read(f, Float64))
        z = ntoh(read(f, Float64))
        push!(vertices, Point3D(x, y, z))
    end
    
    # Read triangle connectivity (3 * n_tri integers)
    faces = Tuple{Int,Int,Int}[]
    if n_tri > 0
        for i in 1:n_tri
            v1 = Int(ntoh(read(f, Int32)))
            v2 = Int(ntoh(read(f, Int32)))
            v3 = Int(ntoh(read(f, Int32)))
            push!(faces, (v1, v2, v3))
        end
    end
    
    # Read quad connectivity (4 * n_quad integers) and split into triangles
    if n_quad > 0
        for i in 1:n_quad
            v1 = Int(ntoh(read(f, Int32)))
            v2 = Int(ntoh(read(f, Int32)))
            v3 = Int(ntoh(read(f, Int32)))
            v4 = Int(ntoh(read(f, Int32)))
            # Split quad into two triangles
            push!(faces, (v1, v2, v3))
            push!(faces, (v1, v3, v4))
        end
    end
    
    close(f)
    println("  Successfully read $(length(vertices)) vertices and $(length(faces)) faces")
    return vertices, faces
end

# --- The O-Grid Extractor ---
function extract_ogrid(vertices, faces, center, radius)
    println("\n--- Extracting O-Grid (Radius: $radius) ---")
    
    valid_faces = Tuple{Int,Int,Int}[]
    
    # 1. Filter Faces by Distance
    for tri in faces
        p1, p2, p3 = vertices[tri[1]], vertices[tri[2]], vertices[tri[3]]
        
        # Calculate Centroid
        cx = (p1.x + p2.x + p3.x) / 3.0
        cy = (p1.y + p2.y + p3.y) / 3.0
        
        dist = sqrt((cx - center[1])^2 + (cy - center[2])^2)
        
        # Keep if inside radius
        if dist < radius
            push!(valid_faces, tri)
        end
    end
    
    if isempty(valid_faces)
        error("Radius too small! No faces found within $radius of $center.")
    end

    # 2. Collect Used Vertices (Compact)
    used_mask = fill(false, length(vertices))
    for f in valid_faces
        used_mask[f[1]] = true
        used_mask[f[2]] = true
        used_mask[f[3]] = true
    end
    
    old_to_new = zeros(Int, length(vertices))
    new_vertices = Point3D[]
    
    for (i, used) in enumerate(used_mask)
        if used
            push!(new_vertices, vertices[i])
            old_to_new[i] = length(new_vertices)
        end
    end
    
    # 3. Remap Face Indices
    new_faces = Tuple{Int,Int,Int}[]
    for f in valid_faces
        push!(new_faces, (old_to_new[f[1]], old_to_new[f[2]], old_to_new[f[3]]))
    end
    
    println("  New Nodes: $(length(new_vertices))")
    println("  New Faces: $(length(new_faces))")
    
    return new_vertices, new_faces
end

# --- Save to SU2 Format ---
function save_su2(filename, vertices, faces)
    println("\nSaving to $filename...")
    open(filename, "w") do io
        println(io, "NDIME= 2")
        println(io, "NELEM= $(length(faces))")
        for (i, f) in enumerate(faces)
            # SU2 Element Type 5 = Triangle
            # Indices are 0-based in SU2
            println(io, "5 $(f[1]-1) $(f[2]-1) $(f[3]-1) $(i-1)")
        end
        println(io, "NPOIN= $(length(vertices))")
        for (i, v) in enumerate(vertices)
            # Output x y point_index (0-based)
            println(io, "$(v.x) $(v.y) $(i-1)") 
        end
        
        # Create a dummy boundary marker for the outer edge
        println(io, "NMARK= 1")
        println(io, "MARKER_TAG= farfield")
        println(io, "MARKER_ELEMS= 1")
        println(io, "3 0 1") # Dummy line element
    end
    println("Done.")
end

# --- Run ---
if !isfile(INPUT_FILE)
    println("Error: Input file '$INPUT_FILE' not found.")
else
    verts, faces = read_ugrid_binary(INPUT_FILE)
    new_verts, new_faces = extract_ogrid(verts, faces, CENTER, RADIUS)
    save_su2(OUTPUT_FILE, new_verts, new_faces)
end