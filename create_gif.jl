using FileIO, ImageIO

frame_dir = "output/animation_frames"
output_file = "greedy_solver_animation.gif"
fps = 5

println("Creating GIF from frames in $frame_dir...")

# Get all PNG files
frame_files = filter(f -> endswith(f, ".png"), readdir(frame_dir))
sort!(frame_files)

println("Found $(length(frame_files)) frames")

# Load all frames
images = []
for (i, fname) in enumerate(frame_files)
    fpath = joinpath(frame_dir, fname)
    img = load(fpath)
    push!(images, img)
    if i % 10 == 0
        println("  Loaded $i/$(length(frame_files)) frames")
    end
end

println("Creating GIF with $fps fps...")
save(output_file, images, fps=fps)

println("✓ Successfully created: $output_file")
file_size = filesize(output_file) / (1024 * 1024)
println("  File size: $(round(file_size, digits=2)) MB")
