# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 09:30:17 2024

@author: holtj
"""
import numpy as np
import matplotlib.pyplot as plt

def create_3d_grid(n, p):
    """Create an n x n x n 3D grid where each cell is open with probability p."""
    return np.random.rand(n, n, n) < p

def dfs(grid, visited, x, y, z):
    """Perform a depth-first search to check for percolation, allowing diagonal moves."""
    n = len(grid)
    if x < 0 or x >= n or y < 0 or y >= n or z < 0 or z >= n or visited[x, y, z] or not grid[x, y, z]:
        return False
    visited[x, y, z] = True
    if x == n - 1:
        return True
    
    directions = [
        (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1), 
        (1, 1, 0), (1, -1, 0), (-1, 1, 0), (-1, -1, 0),
        (1, 0, 1), (1, 0, -1), (-1, 0, 1), (-1, 0, -1),
        (0, 1, 1), (0, 1, -1), (0, -1, 1), (0, -1, -1),
        (1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1),
        (-1, 1, 1), (-1, 1, -1), (-1, -1, 1), (-1, -1, -1)
    ]
    
    return any(dfs(grid, visited, x + dx, y + dy, z + dz) for dx, dy, dz in directions)

def percolates_3d(grid):
    """Check if the 3D grid percolates from top to bottom."""
    n = len(grid)
    visited = np.zeros_like(grid, dtype=bool)
    for i in range(n):
        for j in range(n):
            if grid[0, i, j] and dfs(grid, visited, 0, i, j):
                return 1
    return 0

def simulate_percolation(n, p_values, trials):
    """Simulate the percolation process for a range of probabilities."""
    def run_trial(p):
        return percolates_3d(create_3d_grid(n, p))
    
    def run_for_p(p):
        return np.mean([run_trial(p) for _ in range(trials)])
    
    return list(map(run_for_p, 1-p_values))

def plot_percolation_probability(p_values, percolation_probabilities, n):
    """Plot the percolation probability as a function of the probability p."""
    plt.figure(figsize=(3,3))
    # plt.plot(np.log10(p_values), np.log10(np.asarray(percolation_probabilities)))
    plt.plot((p_values), (np.asarray(percolation_probabilities)), marker='o')
    found = False
    for i in np.arange(0,len(p_values)):
        while percolation_probabilities[i] < 0.5 and found == False:
            plt.vlines(p_values[i], 0, 1.2, linestyle='--')
            print(p_values[i])
            found = True
    plt.xlabel("Density")
    plt.ylabel("Percolation Probability (mean)")
    plt.ylim(0,1.2)
    plt.xlim(0,1)
    plt.grid(False)

def plot_percolation_probability_log(p_values, percolation_probabilities, n):
    """Plot the percolation probability as a function of the probability p."""
    plt.figure(figsize=(3,3))
    plt.plot(np.log10(p_values), np.log10(np.asarray(percolation_probabilities)))
    plt.plot((p_values), (np.asarray(percolation_probabilities)), marker='o')
    plt.xlabel("Density")
    plt.ylabel("Percolation Probability (mean)")
    plt.xlim(-1,0)
    plt.ylim(-1,0.1)
    plt.grid(False)
    
if __name__ == "__main__":
    
    p_values = np.linspace(0.0, 1.0, 100)  # Range of probabilities
    trials = 200  # Number of trials for each probability
    
    n = 3  # Size of the grid     
    percolation_probabilities = simulate_percolation(n, p_values, trials)
    plot_percolation_probability(p_values, percolation_probabilities, n)
    plt.savefig('Figure1G_Perc.svg', dpi=300, facecolor='w', edgecolor='b',
            orientation='portrait', format='svg',
            transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
    plt.show()

    n = 2  # Size of the grid     
    percolation_probabilities = simulate_percolation(n, p_values, trials)
    plot_percolation_probability_log(p_values, percolation_probabilities, n)
    plt.savefig('FigureS3D_L2.svg', dpi=300, facecolor='w', edgecolor='b',
            orientation='portrait', format='svg',
            transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
    plt.show()
    
    n = 3  # Size of the grid     
    percolation_probabilities = simulate_percolation(n, p_values, trials)
    plot_percolation_probability_log(p_values, percolation_probabilities, n)
    plt.savefig('FigureS3D_L3.svg', dpi=300, facecolor='w', edgecolor='b',
            orientation='portrait', format='svg',
            transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
    plt.show()
    
    n = 4  # Size of the grid     
    percolation_probabilities = simulate_percolation(n, p_values, trials)
    plot_percolation_probability_log(p_values, percolation_probabilities, n)
    plt.savefig('FigureS3D_L4.svg', dpi=300, facecolor='w', edgecolor='b',
            orientation='portrait', format='svg',
            transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
    plt.show()
    