# -*- coding: utf-8 -*-
"""
Created on Wed May 29 13:54:44 2024

@author: holtj
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import seaborn as sns
#%%
"""
Model case where Bdello rate of loss is modified proportionally to the fraction of accessible prey 
"""

def model(y, t, c, a, b, f, V, k):
    N, P, B = y
    dNdt = c*N*((300000-N)/300000) - a * P * V * N
    dPdt = k*B*b - V*a*P*N - f*P*V
    dBdt = V*a*P*N-k*B
    return [dNdt, dPdt, dBdt]

#r, maximal growth rate (this gets varied later)
c = 1.1
#a B. bacteriovorus adsorption rate
a=5.52*10**-7
#f B. bacteriovorus loss rate
f = 1.7952
#b prey to predator conversion factor
b = 0.438
#k bdelloplast maturation rate
k = 2.616

# Initial conditions: N0 (prey) and P0 (predator)
N0 = 100000
B0 = 0
P0 = 10**10*0.00067416

y0 = [N0, P0, B0]

t = np.linspace(0, 5, 10000)

# Define ranges and resolution for c and V
resolution = 100

c_values = np.arange(1, 1.6+0.6/resolution, 0.6/resolution)
V_values = np.arange(0+0.6/resolution, 1+0.6/resolution*2, 1/resolution)

# Initialize the matrix to store results
heatmap_data = np.zeros((len(c_values), len(V_values)))
heatmap_data2 = np.zeros((len(c_values), len(V_values)))

p_list = []
n_list = []

# Loop through combinations of c and V
for i, c in enumerate(c_values):
    for j, V in enumerate(V_values):
        # Solve ODE for predator-present case
        # print(c,V)
        solution = odeint(model, y0, t, args=(c, a, b, f, V, k))
        N = solution[:, 0]
        
        # Solve ODE for predator-absent case
        y0_null = [N0, 0, 0]
        solution_null = odeint(model, y0_null, t, args=(c, a, b, f, V, k))
        N_null = solution_null[:, 0]
        
        P1 = np.round(solution[:, 1][-1]+solution[:, 1][-2])
        N1 = np.mean(N_null - N)
        
        p_list.append(P1)
        n_list.append(N1)
        
        # find cumulative difference and store it in the matrix
        heatmap_data[i, j] = N1
        heatmap_data2[i, j] = np.log(P1+0.1)

plt.figure(figsize=(10, 8))
plt.imshow(heatmap_data, aspect='auto', origin='lower', 
           extent=[V_values.min(), V_values.max(), c_values.min(), c_values.max()],
           cmap='YlOrRd')
plt.colorbar(label='Cumulative Prey Difference')
plt.xlabel('V (Interaction Parameter)')
plt.ylabel('c (Intrinsic Growth Rate)')
plt.title('Heatmap of Prey Population Difference')
plt.savefig("FigureS5C.svg", format="svg")
plt.show()

plt.figure(figsize=(10, 8))
plt.imshow(heatmap_data2, aspect='auto', origin='lower', 
           extent=[V_values.min(), V_values.max(), c_values.min(), c_values.max()],
           cmap='BuPu')
plt.colorbar(label='Pred Abundance')
plt.xlabel('V (Interaction Parameter)')
plt.ylabel('c (Intrinsic Growth Rate)')
plt.title('Heatmap of Prey Population Difference')
# plt.savefig("FigureS5E.svg", format="svg")
plt.show()

sns.histplot(x=p_list, y=np.asarray(n_list), bins=25, log_scale=(False, False), cbar=True)
plt.ylabel('Effect on Prey')
plt.xlabel('Predator Abundance')
plt.savefig("FigureS5D.svg", format="svg")
plt.show()

#%%
"""
Null model case
"""
def model(y, t, c, a, b, f, V, k):
    N, P, B = y
    dNdt = c*N*((300000-N)/300000) - a * P * V * N
    dPdt = k*B*b - V*a*P*N - f*P
    dBdt = V*a*P*N-k*B
    return [dNdt, dPdt, dBdt]

# Initialize the matrix to store results
heatmap_data = np.zeros((len(c_values), len(V_values)))
heatmap_data2 = np.zeros((len(c_values), len(V_values)))

p_list = []
n_list = []

# Loop through combinations of c and V
for i, c in enumerate(c_values):
    for j, V in enumerate(V_values):
        # Solve ODE for predator-present case
        # print(c,V)
        solution = odeint(model, y0, t, args=(c, a, b, f, V, k))
        N = solution[:, 0]
        
        # Solve ODE for predator-absent case
        y0_null = [N0, 0, 0]
        solution_null = odeint(model, y0_null, t, args=(c, a, b, f, V, k))
        N_null = solution_null[:, 0]
        
        N1 = np.mean(N_null - N)
        P1 = (solution[:, 1][-1]+solution[:, 2][-1])   
        
        p_list.append(P1)
        n_list.append(N1)
        
        # Compute cumulative difference and store it in the matrix
        heatmap_data[i, j] = N1
        heatmap_data2[i, j] = P1

# Plot heatmap
plt.figure(figsize=(10, 8))
plt.imshow(heatmap_data, aspect='auto', origin='lower', 
           extent=[V_values.min(), V_values.max(), c_values.min(), c_values.max()],
           cmap='YlOrRd')
plt.colorbar(label='Cumulative Prey Difference')
plt.xlabel('V (Interaction Parameter)')
plt.ylabel('c (Intrinsic Growth Rate)')
plt.title('Heatmap of Prey Population Difference')
plt.savefig("FigureS5A.svg", format="svg")
plt.show()

# Plot heatmap
plt.figure(figsize=(10, 8))
plt.imshow(heatmap_data2, aspect='auto', origin='lower', 
           extent=[V_values.min(), V_values.max(), c_values.min(), c_values.max()],
           cmap='BuPu')
plt.colorbar(label='Pred Abundance')
plt.xlabel('V (Interaction Parameter)')
plt.ylabel('c (Intrinsic Growth Rate)')
plt.title('Heatmap of Prey Population Difference')
# plt.savefig("FigureS5B.svg", format="svg")
plt.show()

sns.histplot(x=p_list, y=np.asarray(n_list), bins=25, log_scale=(False, False), cbar=True)
plt.ylabel('Effect on Prey')
plt.xlabel('Predator Abundance')
plt.savefig("FigureS5B.svg", format="svg")
plt.show()

#%%
"""
Model case where predator efficiency is inversely proportion to fraction of accessible prey (V cancels out in dP, dB)
"""
def model(y, t, c, a, b, f, V, k):
    N, P, B = y
    dNdt = c*N*((300000-N)/300000) - a * P * V * N
    dPdt = k*B*b - a*P*N - f*P
    dBdt = a*P*N-k*B
    return [dNdt, dPdt, dBdt]

# Initialize the matrix to store results
heatmap_data = np.zeros((len(c_values), len(V_values)))
heatmap_data2 = np.zeros((len(c_values), len(V_values)))

p_list = []
n_list = []

# Loop through combinations of c and V
for i, c in enumerate(c_values):
    for j, V in enumerate(V_values):
        # Solve ODE for predator-present case
        # print(c,V)
        solution = odeint(model, y0, t, args=(c, a, b, f, V, k))
        N = solution[:, 0]
        
        # Solve ODE for predator-absent case
        y0_null = [N0, 0, 0]
        solution_null = odeint(model, y0_null, t, args=(c, a, b, f, V, k))
        N_null = solution_null[:, 0]
        
        N1 = np.mean(N_null - N)
        P1 = (solution[:, 1][-1]+solution[:, 2][-1])
        
        p_list.append(P1)
        n_list.append(N1)
        
        # Compute cumulative difference and store it in the matrix
        heatmap_data[i, j] = N1
        heatmap_data2[i, j] = P1

# Plot heatmap
plt.figure(figsize=(10, 8))
plt.imshow(heatmap_data, aspect='auto', origin='lower', 
           extent=[V_values.min(), V_values.max(), c_values.min(), c_values.max()],
           cmap='YlOrRd')
plt.colorbar(label='Cumulative Prey Difference')
plt.xlabel('V (Interaction Parameter)')
plt.ylabel('c (Intrinsic Growth Rate)')
plt.title('Heatmap of Prey Population Difference')
plt.savefig("FigureS5E.svg", format="svg")
plt.show()

# Plot heatmap
plt.figure(figsize=(10, 8))
plt.imshow(heatmap_data2, aspect='auto', origin='lower', 
           extent=[V_values.min(), V_values.max(), c_values.min(), c_values.max()],
           cmap='BuPu')
plt.colorbar(label='Pred Abundance')
plt.xlabel('V (Interaction Parameter)')
plt.ylabel('c (Intrinsic Growth Rate)')
plt.title('Heatmap of Prey Population Difference')
# plt.savefig("FigureS5H.svg", format="svg")
plt.show()

sns.histplot(x=p_list, y=np.asarray(n_list), bins=25, log_scale=(False, False), cbar=True)
plt.ylabel('Effect on Prey')
plt.xlabel('Predator Abundance')
plt.savefig("FigureS5F.svg", format="svg")
plt.show()