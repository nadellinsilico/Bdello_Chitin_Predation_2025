# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 09:46:50 2026

@author: holtj
"""
import numpy as np
import matplotlib.pyplot as plt

def create_3d_spherical_grid(n, p):
    """Create a 3D grid (2n x 2n x 2n) with a digital spherical mask."""
    size = 2 * n
    center = (size - 1) / 2

    X, Y, Z = np.indices((size, size, size))
    dist = np.sqrt((X - center)**2 + (Y - center)**2 + (Z - center)**2)

    # digital sphere: max 1-voxel steps
    mask = np.floor(dist + 0.5) <= n

    grid = (np.random.rand(size, size, size) < p) & mask
    return grid, mask

def dfs_iterative(grid, visited, start_x, start_y, start_z, mask):
    """Iterative DFS with 6-neighbor connectivity, restricted to mask."""
    size = grid.shape[0]
    stack = [(start_x, start_y, start_z)]
    directions = [(1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)]
    
    while stack:
        x, y, z = stack.pop()
        if (x < 0 or x >= size or
            y < 0 or y >= size or
            z < 0 or z >= size or
            visited[x,y,z] or
            not grid[x,y,z] or
            not mask[x,y,z]):
            continue
        visited[x,y,z] = True
        for dx, dy, dz in directions:
            stack.append((x+dx, y+dy, z+dz))

def percolates_to_core_3d_sphere(grid, mask, core_size=4):
    """Check if open sites connect the spherical boundary to the central cubic core."""
    size = grid.shape[0]
    visited = np.zeros_like(grid, dtype=bool)
    
    # Central cubic core
    c_start = size//2 - core_size//2
    c_end = c_start + core_size
    core = [(i,j,k) for i in range(c_start, c_end)
                     for j in range(c_start, c_end)
                     for k in range(c_start, c_end)]
    
    # Distance from center
    center = (size - 1) / 2
    X, Y, Z = np.indices(grid.shape)
    distances = np.sqrt((X-center)**2 + (Y-center)**2 + (Z-center)**2)
    
    # Flood fill from **surface boundary only**
    boundary_points = np.array(np.where(
        mask & ((distances >= n - 0.5) & (distances <= n + 0.5))
    ))
    for x, y, z in zip(*boundary_points):
        if grid[x,y,z] and not visited[x,y,z]:
            dfs_iterative(grid, visited, x, y, z, mask)
    
    # Check if any core site was reached
    reaches_core = any(visited[x,y,z] for x,y,z in core)
    return int(reaches_core), visited

def simulate_percolation(n, p_values, trials, core_size=4):
    """Simulate percolation probabilities on spherical lattice."""
    size = 2*n
    center = (size - 1) / 2
    X, Y, Z = np.indices((size, size, size))
    
    dist = np.sqrt((X-center)**2 + (Y-center)**2 + (Z-center)**2)
    mask = np.floor(dist + 0.5) <= n
    
    # mask = (X-center)**2 + (Y-center)**2 + (Z-center)**2 <= n**2
    
    results = []
    median_result = []
    for p in p_values:
        successes = 0
        for _ in range(trials):
            grid = (np.random.rand(size, size, size) < p) & mask
            success, _ = percolates_to_core_3d_sphere(grid, mask, core_size=core_size)
            successes += success
        results.append(successes / trials)
        median_result.append(np.median(np.append(np.ones(successes), np.zeros(trials-successes))))
    return results, median_result

def visualize_3d_realization_sphere(grid, visited, mask):
    """3D visualization of spherical lattice, visited sites, and cubic core."""
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111, projection='3d')
    
    # Sites inside sphere
    x_all, y_all, z_all = np.where(mask)
    ax.scatter(x_all, y_all, z_all, c='gray', s=5, alpha=0.1)
    
    # Visited sites
    x_vis, y_vis, z_vis = np.where(visited)
    ax.scatter(x_vis, y_vis, z_vis, c='gold', s=40, alpha=0.9)
    
    # Core cube
    size = grid.shape[0]
    c = size // 2
    core = [(i,j,k) for i in [c-1,c] for j in [c-1,c] for k in [c-1,c]]
    cx, cy, cz = zip(*core)
    ax.scatter(cx, cy, cz, c='red', s=120)
    
    ax.set_xlim(0,size-1); ax.set_ylim(0,size-1); ax.set_zlim(0,size-1)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    plt.show()


# Example usage
if __name__ == "__main__":
    n = 12
    p = 0.6
    trials = 400
    plt.figure(figsize=(3,3))
    p_values = np.linspace(0,1,200)
    
    percolation_probs = simulate_percolation(n, p_values, trials)
    
    found=False
    for i in np.arange(0,len(p_values)):
        if np.asarray(percolation_probs)[0][i] > 0 and found == False:
                if np.asarray(percolation_probs)[1][i] == 1:
                    plt.vlines(1-p_values[i], 0, 1.2, linestyle='--')
                    print(p_values[i])
                    found = True

    plt.plot(1-p_values, np.asarray(percolation_probs)[0])
    plt.scatter(1-p_values, np.asarray(percolation_probs)[0])
    plt.ylim(0,1.2)
    plt.xlim(0,1)
    plt.savefig('Fig1G.svg', format='svg')
    plt.show()
    #Visualize an example
    grid, mask = create_3d_spherical_grid(n, p)
    success, visited = percolates_to_core_3d_sphere(grid, mask)
    print("Reaches core:", success)
    visualize_3d_realization_sphere(grid, visited, mask)
    