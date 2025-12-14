"""
Create animated GIF from frames using Python and PIL/Pillow
This is a fallback method when Julia's native GIF creation doesn't work
"""

from PIL import Image
import os
from pathlib import Path

def create_gif_from_frames(frame_dir, output_file, fps=5, optimize=True):
    """
    Create an animated GIF from PNG frames
    
    Args:
        frame_dir: Directory containing frame_####.png files
        output_file: Output GIF filename
        fps: Frames per second (default: 5)
        optimize: Whether to optimize GIF size (default: True)
    """
    print(f"\n--- Creating Animated GIF with Python ---")
    print(f"Frame directory: {frame_dir}")
    
    # Get all PNG files
    frame_files = sorted([f for f in os.listdir(frame_dir) if f.endswith('.png')])
    
    if not frame_files:
        print(f"ERROR: No PNG frames found in {frame_dir}")
        return False
    
    print(f"Found {len(frame_files)} frames")
    
    # Load all frames
    images = []
    print("Loading frames...")
    for i, fname in enumerate(frame_files):
        fpath = os.path.join(frame_dir, fname)
        img = Image.open(fpath)
        # Convert to RGB if necessary (GIF doesn't support RGBA)
        if img.mode == 'RGBA':
            # Create white background
            bg = Image.new('RGB', img.size, (255, 255, 255))
            bg.paste(img, mask=img.split()[3])  # Use alpha channel as mask
            img = bg
        elif img.mode != 'RGB':
            img = img.convert('RGB')
        images.append(img)
        
        if (i + 1) % 10 == 0:
            print(f"  Loaded {i+1}/{len(frame_files)} frames")
    
    # Calculate duration per frame in milliseconds
    duration = int(1000 / fps)
    
    print(f"Creating GIF with {fps} fps (duration: {duration}ms per frame)...")
    
    # Save as animated GIF
    images[0].save(
        output_file,
        save_all=True,
        append_images=images[1:],
        duration=duration,
        loop=0,  # 0 = infinite loop
        optimize=optimize
    )
    
    print(f"✓ Successfully created: {output_file}")
    
    # Report file size
    file_size = os.path.getsize(output_file) / (1024 * 1024)  # MB
    print(f"  File size: {file_size:.2f} MB")
    
    return True

if __name__ == "__main__":
    frame_dir = "output/animation_frames"
    output_file = "greedy_solver_animation.gif"
    fps = 5
    
    # Check if frame directory exists
    if not os.path.exists(frame_dir):
        print(f"ERROR: Frame directory not found: {frame_dir}")
        print("Run main.jl first to generate frames")
        exit(1)
    
    success = create_gif_from_frames(frame_dir, output_file, fps=fps)
    
    if success:
        print("\n✓ Animation complete!")
        print(f"  View your animation: {output_file}")
    else:
        print("\n✗ Failed to create animation")
        exit(1)
