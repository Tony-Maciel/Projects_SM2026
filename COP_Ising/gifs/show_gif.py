import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

def create_ising_timelapse(input_file, output_gif, fps=10):
    if not os.path.exists(input_file):
        print(f"Error: File '{input_file}' not found.")
        return

    print(f"Reading data from {input_file}...")
    times = []
    configs = []
    
    # Parse file
    with open(input_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            
            times.append(float(parts[0]))
            # Spins start from the second column
            configs.append(np.array([int(s) for s in parts[1:]], dtype=np.int8))
            
    times = np.array(times)
    configs = np.array(configs)
    
    # Determine system size 
    num_spins = configs.shape[1]
    L = int(np.round(np.sqrt(num_spins))) # assumes 2D
    
    print(f"Detected system size: L = {L} (Total spins = {num_spins})")
    print(f"Number of frames to render: {len(times)}")

    # Setup the 2D plot using imshow
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # cmap='gray_r' maps the minimum value (-1) to white and maximum (1) to black
    # .T transposes the matrix, and origin='lower' ensures the coordinate orientation 
    initial_grid = configs[0].reshape((L, L)).T
    im = ax.imshow(initial_grid, cmap='gray_r', vmin=-1, vmax=1, origin='lower')
    
    ax.set_title("2D COP Ising Model Timelapse")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    
    # Initialize time label
    time_text = ax.text(0.05, 0.95, "", transform=ax.transAxes, fontsize=12, 
                        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
                        va='top')

    # Animation update function
    def update(frame):
        # Update the image data directly (much faster than scatter)
        grid = configs[frame].reshape((L, L)).T
        im.set_array(grid)
        
        # Update the time display
        time_text.set_text(f"Time: {times[frame]:.2e}")
        
        # Print progress for long renderings
        if frame % max(1, len(times) // 10) == 0:
            print(f"Rendering frame {frame}/{len(times)}...")
            
        return [im, time_text]

    print(f"Generating animation ({fps} FPS)...")
    # blit=True optimizes rendering by only redrawing changing parts
    ani = animation.FuncAnimation(fig, update, frames=len(times), interval=1000/fps, blit=True)
    
    print(f"Saving GIF to {output_gif}...")
    ani.save(output_gif, writer='pillow', fps=fps)

if __name__ == "__main__":
    L = "512"
    input_file  = f"configs_{L}_0.600_57.dat"
    output_file = f"timelapse{L}.gif"
    fps         = 10
    create_ising_timelapse(input_file, output_file, fps)
