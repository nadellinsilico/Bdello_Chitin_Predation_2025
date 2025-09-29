# -*- coding: utf-8 -*-
"""
Created on Fri May 31 10:38:49 2024

@author: holtj
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

np.random.seed(1889)

def create_grid(n, p):
    """
    Create an n x n grid where each cell is open with probability p.
    
    Args:
    n (int): Size of the grid (n x n).
    p (float): Probability that a cell is open.
    
    Returns:
    np.ndarray: n x n grid with open (True) and blocked (False) cells.
    """
    grid = np.random.rand(n, n) < p
    # grid[0, :] = True  # Ensure the top row is always open
    return grid

def percolates(grid):
    """
    Check if the grid percolates (i.e., if there's a path from the top to the bottom).
    
    Args:
    grid (np.ndarray): n x n grid with open (True) and blocked (False) cells.
    
    Returns:
    np.ndarray: n x n grid with the percolation path marked.
    """
    n = len(grid)
    visited = np.zeros_like(grid, dtype=bool)
    
    def dfs(x, y):
        if x < 0 or x >= n or y < 0 or y >= n or visited[x, y] or not grid[x, y]:
            return
        visited[x, y] = True
        dfs(x + 1, y)
        dfs(x - 1, y)
        dfs(x, y + 1)
        dfs(x, y - 1)
    
    for i in range(n):
        if grid[0, i] and not visited[0, i]:
            dfs(0, i)
    
    return visited

def visualize_percolation(grid, visited, p):
    """
    Visualize the percolation grid and path.
    
    Args:
    grid (np.ndarray): n x n grid with open (True) and blocked (False) cells.
    visited (np.ndarray): n x n grid with the percolation path marked.
    """
    n = len(grid)
    plt.figure(figsize=(8, 8))
    cmap_custom = matplotlib.colors.ListedColormap(['purple', 'black'])
    plt.imshow(grid, cmap=cmap_custom, origin='upper')
    
    # Highlight the percolation path
    for i in range(n):
        for j in range(n):
            if visited[i, j]:
                plt.scatter(j, i, color='yellow', s=800)
    
    plt.title("Percolation Path Visualization")
    plt.xlabel("Column")
    plt.ylabel("Row")

#%%
n = 4       # Size of the grid
p = 0.3      # Packing probability

grid = create_grid(n, p)
visited = percolates(grid)
visualize_percolation(grid, visited, p)
plt.savefig('FigureS3B.svg', format="svg")
plt.show()

#%%
n = 4       # Size of the grid
p = 0.5      # Packing probability

grid = create_grid(n, p)
visited = percolates(grid)
visualize_percolation(grid, visited, p)
plt.savefig('FigureS3C.svg', format="svg")
plt.show()
