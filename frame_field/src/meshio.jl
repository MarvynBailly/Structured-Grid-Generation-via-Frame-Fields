using FileIO
using GeometryBasics

"""
    load_msh_4(filepath::String)

Custom parser for MSH format 4.x files.
Handles the MSH 4.1 format that MeshIO has trouble with.
"""
function load_msh_4(filepath::String)
    vertices = Point3f[]
    triangles = TriangleFace{Int}[]
    
    open(filepath, "r") do io
        # Read through file
        in_nodes = false
        in_elements = false
        node_blocks = 0
        nodes_to_read = 0
        current_node_id = 0
        node_map = Dict{Int, Int}()  # Maps node ID to index
        
        elem_blocks = 0
        elems_to_read = 0
        
        while !eof(io)
            line = readline(io)
            
            # Start of Nodes section
            if startswith(line, "\$Nodes")
                in_nodes = true
                # Next line has: numEntityBlocks numNodes minNodeTag maxNodeTag
                line = readline(io)
                parts = split(line)
                node_blocks = parse(Int, parts[1])
                total_nodes = parse(Int, parts[2])
                continue
            end
            
            # End of Nodes section
            if startswith(line, "\$EndNodes")
                in_nodes = false
                continue
            end
            
            # Process node blocks
            if in_nodes
                if nodes_to_read == 0
                    # Block header: entityDim entityTag parametric numNodesInBlock
                    parts = split(line)
                    if length(parts) >= 4
                        nodes_to_read = parse(Int, parts[4])
                        # Read node IDs first
                        node_ids = Int[]
                        for _ in 1:nodes_to_read
                            nid = parse(Int, readline(io))
                            push!(node_ids, nid)
                        end
                        # Then read coordinates
                        for nid in node_ids
                            coord_line = readline(io)
                            coords = split(coord_line)
                            x = parse(Float32, coords[1])
                            y = parse(Float32, coords[2])
                            z = parse(Float32, coords[3])
                            push!(vertices, Point3f(x, y, z))
                            node_map[nid] = length(vertices)
                        end
                        nodes_to_read = 0
                    end
                end
                continue
            end
            
            # Start of Elements section
            if startswith(line, "\$Elements")
                in_elements = true
                # Next line has: numEntityBlocks numElements minElementTag maxElementTag
                line = readline(io)
                parts = split(line)
                elem_blocks = parse(Int, parts[1])
                continue
            end
            
            # End of Elements section
            if startswith(line, "\$EndElements")
                in_elements = false
                continue
            end
            
            # Process element blocks
            if in_elements
                if elems_to_read == 0
                    # Block header: entityDim entityTag elementType numElementsInBlock
                    parts = split(line)
                    if length(parts) >= 4
                        elem_type = parse(Int, parts[3])
                        elems_to_read = parse(Int, parts[4])
                        
                        # Element type 2 = 3-node triangle
                        if elem_type == 2
                            for _ in 1:elems_to_read
                                elem_line = readline(io)
                                parts = split(elem_line)
                                # Format: elemTag node1 node2 node3
                                n1 = parse(Int, parts[2])
                                n2 = parse(Int, parts[3])
                                n3 = parse(Int, parts[4])
                                
                                # Convert to 1-based indices
                                i1 = node_map[n1]
                                i2 = node_map[n2]
                                i3 = node_map[n3]
                                
                                push!(triangles, TriangleFace{Int}(i1, i2, i3))
                            end
                        else
                            # Skip non-triangle elements
                            for _ in 1:elems_to_read
                                readline(io)
                            end
                        end
                        elems_to_read = 0
                    end
                end
                continue
            end
        end
    end
    
    return GeometryBasics.Mesh(vertices, triangles)
end

"""
    load_triangulation(filepath::String)

Loads a mesh and standardizes it into a GeometryBasics Mesh 
containing only 1-based integer TriangleFaces.

Falls back to custom MSH 4.x parser if MeshIO fails.
"""
function load_triangulation(filepath::String)
    if !isfile(filepath)
        error("File not found: $filepath")
    end

    # Try custom MSH parser first for .msh files
    if endswith(lowercase(filepath), ".msh")
        try
            return load_msh_4(filepath)
        catch e
            println("Custom MSH parser failed: $e")
            println("Trying MeshIO...")
        end
    end
    
    # Fall back to MeshIO
    raw_mesh = load(filepath)

    # Standardize Positions and Faces
    verts = decompose(Point3f, raw_mesh)
    faces_idx = decompose(TriangleFace{Int}, raw_mesh)

    return GeometryBasics.Mesh(verts, faces_idx)
end
