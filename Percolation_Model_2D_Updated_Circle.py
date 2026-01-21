# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 13:04:12 2026

@author: holtj
"""

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def create_circular_grid(n, p):
    size = 2 * n
    grid = np.zeros((size, size), dtype=bool)

    center = (size - 1) / 2
    radius = n

    for i in range(size):
        for j in range(size):
            if (i - center)**2 + (j - center)**2 <= radius**2:
                grid[i, j] = np.random.rand() < p

    return grid

def percolates_to_core(grid, core_size=4):
    size = grid.shape[0]
    visited = np.zeros_like(grid, dtype=bool)
    
    # Define central square core
    center = size // 2
    core_start = center - core_size // 2
    core_indices = [(i, j) for i in range(core_start, core_start + core_size)
                            for j in range(core_start, core_start + core_size)]
    
    def dfs(x, y):
        if x < 0 or x >= size or y < 0 or y >= size:
            return
        if visited[x, y] or not grid[x, y]:
            return
        visited[x, y] = True
        dfs(x + 1, y)
        dfs(x - 1, y)
        dfs(x, y + 1)
        dfs(x, y - 1)
    
    # Start DFS from all boundary sites in the circle
    center_coord = (size - 1) / 2
    radius = n
    for i in range(size):
        for j in range(size):
            if grid[i, j] and not visited[i, j]:
                distance = np.sqrt((i - center_coord)**2 + (j - center_coord)**2)
                if distance >= radius - 1:  # near the circle edge
                    dfs(i, j)
    
    reaches_core = any(visited[x, y] for x, y in core_indices)
    return reaches_core, visited

def visualize_circular_percolation(grid, visited, p, core_size=4):
    size = grid.shape[0]
    plt.figure(figsize=(8, 8))
    
    # Mask points outside the circle
    center = (size - 1) / 2
    radius = n
    mask = np.ones_like(grid, dtype=bool)
    for i in range(size):
        for j in range(size):
            if (i - center)**2 + (j - center)**2 <= radius**2:
                mask[i, j] = False
    masked_grid = np.ma.masked_array(grid, mask=mask)
    
    
    core_start = center - core_size // 2
    rect = patches.Rectangle(
        (core_start, core_start),
        core_size, core_size,
        linewidth=4,
        edgecolor='black',
        facecolor='none',
        linestyle='--'
    )
    plt.gca().add_patch(rect)
    
    plt.title(f"packing = {p}, radius=6um, pixel=0.5um")
    plt.xticks([])
    plt.yticks([])
    plt.axis('equal')

    # Plot grid
    cmap_custom = matplotlib.colors.ListedColormap(['purple', 'grey'])
    plt.imshow(masked_grid, cmap=cmap_custom, origin='upper')
    plt.savefig('FigS41.svg', format='svg')
    plt.show()
    
    plt.figure(figsize=(8, 8))
    cmap_custom = matplotlib.colors.ListedColormap(['purple', 'grey'])
    plt.imshow(masked_grid, cmap=cmap_custom, origin='upper')
    # Plot visited sites
    for i in range(size):
        for j in range(size):
            if visited[i, j]:
                plt.scatter(j, i, color='yellow', s=200)
    
    # Draw square core
    core_start = center - core_size // 2
    rect = patches.Rectangle(
        (core_start, core_start),
        core_size, core_size,
        linewidth=3,
        edgecolor='black',
        facecolor='none',
        linestyle='--'
    )
    plt.gca().add_patch(rect)
    
    plt.title(f"packing = {p}, radius=6um, pixel=0.5um")
    plt.xticks([])
    plt.yticks([])
    plt.axis('equal')

# Example usage
n = 12        # Half-grid size
p = 0.6

grid = create_circular_grid(n, p)
reaches_core, visited = percolates_to_core(grid)

visualize_circular_percolation(grid, visited, p)
plt.savefig('FigS42.svg', format='svg')
plt.show()