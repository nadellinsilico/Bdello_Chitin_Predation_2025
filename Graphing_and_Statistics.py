# -*- coding: utf-8 -*-
"""
Code used to generate graphs and stats. 

Graphs were aesthetically modified in Adobe Illustrator after being generated here. 
@author: holtj
"""
#%%
"""
Libraries 
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from scipy import stats
import seaborn as sns
from scipy.integrate import solve_ivp
#%%
"""
Functions
"""
"""
Used to make pop dynamics graphs in main text figures
"""
def average(df, subject, times):
    average_m = []
    average_std = []
    time_return = []
    for t2 in times:
        df2 = df.drop(df[(df.Day != t2)].index, inplace = False)
        print(df2)
        if len(df2) > 0:
            
            df2 = df2.apply(lambda x: pd.to_numeric(x, errors='ignore')).dropna()
            average_m.append(np.mean(df2[subject]))
            average_std.append(stats.sem(df2[subject]))
            time_return.append(t2)
    return((average_m, average_std, time_return))

"""
Used in the box and whisker plots to visualize individual data points
"""
def simple_beeswarm2(y, nbins=10, width=.25):
    """
    Returns x coordinates for the points in ``y``, so that plotting ``x`` and
    ``y`` results in a bee swarm plot.
    https://stackoverflow.com/questions/36153410/how-to-create-a-swarm-plot-with-matplotlib
    """
    y = np.asarray(y)
    print(y)
    if nbins is None:
        nbins = len(y) // 6

    # Get upper bounds of bins
    x = np.zeros(len(y))
    ylo = np.min(y)
    yhi = np.max(y)
    dy = (yhi - ylo) / nbins
    ybins = np.linspace(ylo + dy, yhi - dy, nbins - 1)

    # Divide indices into bins
    i = np.arange(len(y))
    ibs = [0] * nbins
    ybs = [0] * nbins
    nmax = 0
    for j, ybin in enumerate(ybins):
        f = y <= ybin
        ibs[j], ybs[j] = i[f], y[f]
        nmax = max(nmax, len(ibs[j]))
        f = ~f
        i, y = i[f], y[f]
    ibs[-1], ybs[-1] = i, y
    nmax = max(nmax, len(ibs[-1]))

    # Assign x indices
    try:
        dx = width / (nmax // 2)
    except:
        dx = 0.25
    for i, y in zip(ibs, ybs):
        if len(i) > 1:
            j = len(i) % 2
            i = i[np.argsort(y)]
            a = i[j::2]
            b = i[j+1::2]
            x[a] = (0.5 + j / 3 + np.arange(len(b))) * dx
            x[b] = (0.5 + j / 3 + np.arange(len(b))) * -dx
    return(x)

def weighted_histogram(data, bins, weighted=True):
    x = []
    y = []
    total = np.sum(data['Vol'])
    for i in np.arange(0,len(bins)-1):
        # print(len(data))
        data2 = data.drop(data[(data.Density >= bins[i+1])|(data.Density < bins[i])].index, inplace = False)
        # print(len(data2))
        x.append(bins[i])
        if weighted == True:
            if i == 0:
                y.append(np.sum(data2['Vol'])/total)
            else:
                y.append(np.sum(data2['Vol'])/total + y[i-1])
        else:
            y.append(np.sum(data2['Vol'])/total)
    return([x,y])
#%%
"""
Main Text Figures
"""
#%%
"""
Figure1A

Predation of chitin formed biofilms
"""
#load data from excel
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1A.1')
df_pred =  pd.read_excel('SI_Data.xlsx', sheet_name='Fig1A.2') 

times= np.arange(3,7,1)
v_mean, v_std, times = average(df_pred, 'Prey', times)
p_mean, p_std, times = average(df_pred, 'Pred', times)
v_mean_d, v_std_d, times = average(df, 'Prey', times)

print(np.sum(np.asarray(v_mean_d)-np.asarray(v_mean))/len(times))
#degrees of freedom for scipy.stats sem defaults to 1
print(stats.sem(np.asarray(v_mean_d)-np.asarray(v_mean)))
    
fig, ax1 = plt.subplots(figsize=(3, 4))
ax2 = ax1.twinx()

ax1.plot(times, p_mean, 'gold')
ax1.errorbar(x=times, y=p_mean, c='gold', fmt='p', yerr=p_std, capsize=5, elinewidth=1, markeredgewidth=1)

ax2.plot(times, v_mean, 'purple')
ax2.errorbar(x=times, y=v_mean, c='purple', fmt='p', yerr=v_std, capsize=5, elinewidth=1, markeredgewidth=1)

ax2.plot(times, v_mean_d, 'mediumorchid', linestyle='--') 
ax2.errorbar(x=times, y=v_mean_d, c='mediumorchid', fmt='p', yerr=v_std_d, capsize=5, elinewidth=1, markeredgewidth=1) 

ax2.set_xlim(2.75, 6.25)
ax1.set_xlim(2.75, 6.25)
ax1.set_ylim(0,16000)
ax2.set_ylim(0,600000)
ax2.set_ylabel('Prey (um)')
ax1.set_ylabel('Pred (um)')
ax1.set_xlabel('Time d')
plt.title('Chitin, WT')
plt.savefig('Figure1A.svg', format='svg')
plt.show()
#%%
"""
Figure1B

Predation of chitobiose formed biofilms
"""

df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1B.1')
df_pred =  pd.read_excel('SI_Data.xlsx', sheet_name='Fig1B.2') 

times= np.arange(3,6,1)
v_mean, v_std, times = average(df_pred, 'Prey', times)
p_mean, p_std, times = average(df_pred, 'Pred', times)
v_mean_d, v_std_d, times = average(df, 'Prey', times)

print(np.sum(np.asarray(v_mean_d)-np.asarray(v_mean))/len(v_mean_d))
print(stats.sem(np.asarray(v_mean_d)-np.asarray(v_mean)))

fig, ax1 = plt.subplots(figsize=(3, 4))
ax2 = ax1.twinx()

ax1.plot(times, p_mean, 'gold')
ax1.errorbar(x=times, y=p_mean, c='gold', fmt='s', yerr=p_std, capsize=5, elinewidth=1, markeredgewidth=1)

ax2.plot(times, v_mean, 'purple')
ax2.errorbar(x=times, y=v_mean, c='purple', fmt='s', yerr=v_std, capsize=5, elinewidth=1, markeredgewidth=1)

ax2.plot(times, v_mean_d, 'mediumorchid', linestyle='--') 
ax2.errorbar(x=times, y=v_mean_d, c='mediumorchid', fmt='s', yerr=v_std_d, capsize=5, elinewidth=1, markeredgewidth=1) 


ax2.set_xlim(2.75, 5.25)
ax1.set_xlim(2.75, 5.25)
ax1.set_ylim(0,16000)
ax2.set_ylim(0,600000)
ax2.set_ylabel('Prey')
ax1.set_ylabel('Pred')
ax1.set_xlabel('Time d')
plt.title('Chitobiose, WT')
plt.savefig('Figure1B.svg', format='svg')
plt.show()
#%%
"""
Figure 1C

Predation on GlcNAc
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1C.1')
df_pred =  pd.read_excel('SI_Data.xlsx', sheet_name='Fig1C.2')

times= np.arange(3,6,1)
v_mean, v_std, times = average(df_pred, 'Prey', times)
p_mean, p_std, times = average(df_pred, 'Pred', times)
v_mean_d, v_std_d, times = average(df, 'Prey', times)

print(np.sum(np.asarray(v_mean_d)-np.asarray(v_mean))/len(v_mean_d))
print(stats.sem(np.asarray(v_mean_d)-np.asarray(v_mean)))

fig, ax1 = plt.subplots(figsize=(3, 4))
ax2 = ax1.twinx()

ax1.plot(times, p_mean, 'gold')
ax1.errorbar(x=times, y=p_mean, c='gold', fmt='o', yerr=p_std, capsize=5, elinewidth=1, markeredgewidth=1)

ax2.plot(times, v_mean, 'purple')
ax2.errorbar(x=times, y=v_mean, c='purple', fmt='o', yerr=v_std, capsize=5, elinewidth=1, markeredgewidth=1)

ax2.plot(times, v_mean_d, 'mediumorchid', linestyle='--') 
ax2.errorbar(x=times, y=v_mean_d, c='mediumorchid', fmt='o', yerr=v_std_d, capsize=5, elinewidth=1, markeredgewidth=1) 


ax2.set_xlim(2.75, 5.25)
ax1.set_xlim(2.75, 5.25)
ax1.set_ylim(0,16000)
ax2.set_ylim(0,600000)
ax2.set_ylabel('Prey')
ax1.set_ylabel('Pred')
ax1.set_xlabel('Time d')
plt.title('GlcNAc, WT')
plt.savefig('Figure1C.svg', format='svg')
plt.show()
#%%
"""
Figure1D-F

images exported from Zen
"""
#%%
"""
Figure 1G 

Observed Predation Locations

See percolation model for simulation data
"""

df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1G')
df = df.fillna(0)
bins = np.arange(0,1,0.02)

y_data = []
bins_sample = []

counter = 0
for i in np.arange(0,len(bins)-1):
    data2 = df.drop(df[(df['Density'] != bins[i])].index, inplace=False)
    if len(data2) > 1:
        y_data.append(np.mean(data2['Predation']))
        bins_sample.append(bins[i])
    if counter == 0:
        if np.median(data2['Predation']) == 0:
            counter+=1
            vline = bins[i]

plt.figure(figsize=(3,3))
plt.plot(np.asarray(bins_sample[:]),y_data, color='black')
plt.scatter(bins_sample[:], y_data, color='black')
plt.grid(None)
plt.ylabel('Predation Probability')
plt.xlabel('Cell Packing')
plt.vlines(vline, 0, 1.2, color='black', linestyle='--')
plt.xlim(0,1)
plt.ylim(0,1.2)
plt.grid(None)
plt.savefig('Figure1G.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
Figure 1H 

local density histograms
"""

f1 = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1H_Chitobiose')
f2 = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1H_Chitin')
f3 = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1H_GlcNAc')

#%%
plt.figure(figsize=(3,3))
sns.histplot(data=f1, x='Density', bins=50, color='#FAA51A', alpha=0.5, stat='probability', label ='GlcNAc2, WT', cumulative=False, weights='Vol')
sns.histplot(data=f2, x='Density', bins=50, color='#6CCBD8', alpha=0.5, stat='probability', label ='Chitin, WT', cumulative=False, weights='Vol')
sns.histplot(data=f3, x='Density', bins=50, color='#7D2880', alpha=0.5, stat='probability', label ='GlcNAc, WT', cumulative=False, weights='Vol')
plt.ylabel('Volume Weighted Frequency')
plt.xlabel('Cell Packing (6um)')
plt.ylim(0,0.20)
plt.xlim(0,1)

plt.savefig('Figure1H.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
Figure 1I
 
See ODE Model for heat map and percolation model for perc prob and bins
"""

#pulled from percolation model
bins = [0.        , 0.00502513, 0.01005025, 0.01507538, 0.0201005 ,
       0.02512563, 0.03015075, 0.03517588, 0.04020101, 0.04522613,
       0.05025126, 0.05527638, 0.06030151, 0.06532663, 0.07035176,
       0.07537688, 0.08040201, 0.08542714, 0.09045226, 0.09547739,
       0.10050251, 0.10552764, 0.11055276, 0.11557789, 0.12060302,
       0.12562814, 0.13065327, 0.13567839, 0.14070352, 0.14572864,
       0.15075377, 0.15577889, 0.16080402, 0.16582915, 0.17085427,
       0.1758794 , 0.18090452, 0.18592965, 0.19095477, 0.1959799 ,
       0.20100503, 0.20603015, 0.21105528, 0.2160804 , 0.22110553,
       0.22613065, 0.23115578, 0.2361809 , 0.24120603, 0.24623116,
       0.25125628, 0.25628141, 0.26130653, 0.26633166, 0.27135678,
       0.27638191, 0.28140704, 0.28643216, 0.29145729, 0.29648241,
       0.30150754, 0.30653266, 0.31155779, 0.31658291, 0.32160804,
       0.32663317, 0.33165829, 0.33668342, 0.34170854, 0.34673367,
       0.35175879, 0.35678392, 0.36180905, 0.36683417, 0.3718593 ,
       0.37688442, 0.38190955, 0.38693467, 0.3919598 , 0.39698492,
       0.40201005, 0.40703518, 0.4120603 , 0.41708543, 0.42211055,
       0.42713568, 0.4321608 , 0.43718593, 0.44221106, 0.44723618,
       0.45226131, 0.45728643, 0.46231156, 0.46733668, 0.47236181,
       0.47738693, 0.48241206, 0.48743719, 0.49246231, 0.49748744,
       0.50251256, 0.50753769, 0.51256281, 0.51758794, 0.52261307,
       0.52763819, 0.53266332, 0.53768844, 0.54271357, 0.54773869,
       0.55276382, 0.55778894, 0.56281407, 0.5678392 , 0.57286432,
       0.57788945, 0.58291457, 0.5879397 , 0.59296482, 0.59798995,
       0.60301508, 0.6080402 , 0.61306533, 0.61809045, 0.62311558,
       0.6281407 , 0.63316583, 0.63819095, 0.64321608, 0.64824121,
       0.65326633, 0.65829146, 0.66331658, 0.66834171, 0.67336683,
       0.67839196, 0.68341709, 0.68844221, 0.69346734, 0.69849246,
       0.70351759, 0.70854271, 0.71356784, 0.71859296, 0.72361809,
       0.72864322, 0.73366834, 0.73869347, 0.74371859, 0.74874372,
       0.75376884, 0.75879397, 0.7638191 , 0.76884422, 0.77386935,
       0.77889447, 0.7839196 , 0.78894472, 0.79396985, 0.79899497,
       0.8040201 , 0.80904523, 0.81407035, 0.81909548, 0.8241206 ,
       0.82914573, 0.83417085, 0.83919598, 0.84422111, 0.84924623,
       0.85427136, 0.85929648, 0.86432161, 0.86934673, 0.87437186,
       0.87939698, 0.88442211, 0.88944724, 0.89447236, 0.89949749,
       0.90452261, 0.90954774, 0.91457286, 0.91959799, 0.92462312,
       0.92964824, 0.93467337, 0.93969849, 0.94472362, 0.94974874,
       0.95477387, 0.95979899, 0.96482412, 0.96984925, 0.97487437,
       0.9798995 , 0.98492462, 0.98994975, 0.99497487, 1.        ]

#pulled from percolation model
perc_prob = [1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    ,
       1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    ,
       1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    ,
       1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    ,
       1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    ,
       1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    ,
       1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    ,
       1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    ,
       1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    ,
       1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    ,
       1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    ,
       1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    ,
       1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    ,
       1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    ,
       1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 1.    ,
       1.    , 1.    , 1.    , 1.    , 1.    , 1.    , 0.9975, 1.    ,
       0.9975, 0.9975, 0.995 , 0.99  , 0.9825, 0.9725, 0.945 , 0.9225,
       0.935 , 0.8575, 0.81  , 0.775 , 0.665 , 0.61  , 0.535 , 0.435 ,
       0.395 , 0.3   , 0.25  , 0.1925, 0.165 , 0.105 , 0.07  , 0.06  ,
       0.05  , 0.04  , 0.0175, 0.025 , 0.01  , 0.0025, 0.0025, 0.0025,
       0.0025, 0.    , 0.    , 0.0025, 0.0025, 0.    , 0.    , 0.    ,
       0.    , 0.    , 0.    , 0.    , 0.    , 0.    , 0.    , 0.    ,
       0.    , 0.    , 0.    , 0.    , 0.    , 0.    , 0.    , 0.    ,
       0.    , 0.    , 0.    , 0.    , 0.    , 0.    , 0.    , 0.    ,
       0.    , 0.    , 0.    , 0.    , 0.    , 0.    , 0.    , 0.    ]

chitin_predated = 0 
glcnac_predated = 0

for i in np.arange(0,len(bins)-1):
    data_chitin = f2.drop(f2[(f2['Density'] < bins[i])|(f2['Density'] >= bins[i+1])].index, inplace=False)
    data_glcnac = f3.drop(f3[(f3['Density'] < bins[i])|(f3['Density'] >= bins[i+1])].index, inplace=False)
    chitin_predated+= perc_prob[i+1]*np.sum(data_chitin['Vol'])
    glcnac_predated+= perc_prob[i+1]*np.sum(data_glcnac['Vol'])

#these printed values are used in Figure 1 
print(chitin_predated/np.sum(f2['Vol']))    
print(glcnac_predated/np.sum(f3['Vol']))    
#%%
"""
Figure 2A

Cumulative histogram of densities 
"""
f2 = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1H_Chitin')
f3 = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1H_GlcNAc')
f4 =pd.read_excel('SI_Data.xlsx', sheet_name='Fig2A_rbmA')
fv = pd.read_excel('SI_Data.xlsx', sheet_name='Fig2A_vpsL')

#%%
plt.figure(figsize=(3,4))

data = weighted_histogram(f2, bins=np.arange(0,1+1/50,1/50))
x = data[0]
y = data[1]

plt.plot(x,y,color='#6CCBD8', alpha=1, label ='Chitin, WT')

data = weighted_histogram(f3, bins=np.arange(0,1+1/50,1/50))
x = data[0]
y = data[1]

plt.plot(x,y,color='purple', alpha=1, label ='Chitin, WT')

data = weighted_histogram(f4, bins=np.arange(0,1+1/50,1/50))
x = data[0]
y = data[1]

plt.plot(x,y,color='black', alpha=1, label ='Chitin, WT')

data = weighted_histogram(fv, bins=np.arange(0,1+1/50,1/50))
x = data[0]
y = data[1]

plt.plot(x,y,color='firebrick', alpha=1, label ='Chitin, WT')

plt.ylabel('Volume Weighted Frequency')
plt.xlabel('Cell Packing (6um)')
plt.ylim(0,1.01)
plt.xlim(-0.01,1.01)

plt.savefig('Figure2A.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
Figure 2B

biofilm growth curves of strains in GlcNAc
"""

df_rbma = pd.read_excel('SI_Data.xlsx', sheet_name='Fig2B_rbmA')
df_vpsL = pd.read_excel('SI_Data.xlsx', sheet_name='Fig2B_vpsL') 
df_WT = pd.read_excel('SI_Data.xlsx', sheet_name='Fig2B_GlcNAc') 
df_Chitin = pd.read_excel('SI_Data.xlsx', sheet_name='Fig2B_Chitin') 

times= np.arange(0,4,1)

v_mean, v_std, times = average(df_rbma, 'Biovolume', times=times)
v_mean_d, v_std_d, times = average(df_vpsL, 'Biovolume', times=times)
v_mean_wt, v_std_wt, times = average(df_WT, 'Biovolume', times=times)
v_mean_R, v_std_R, times = average(df_Chitin, 'Biovolume', times=times)

fig, ax1 = plt.subplots(figsize=(4, 4))
ax1.plot(times, v_mean, 'black', label='rbmA')
ax1.errorbar(x=times, y=v_mean, c='black', fmt='o', yerr=v_std, capsize=5, elinewidth=1, markeredgewidth=1)

ax1.plot(times, v_mean_d, 'firebrick', label='vpsL') 
ax1.errorbar(x=times, y=v_mean_d, c='firebrick', fmt='o', yerr=v_std_d, capsize=5, elinewidth=1, markeredgewidth=1) 

ax1.plot(times, v_mean_R, 'cyan', label='WT_Chitin') 
ax1.errorbar(x=times, y=v_mean_R, c='cyan', fmt='p', yerr=v_std_R, capsize=5, elinewidth=1, markeredgewidth=1) 

ax1.plot(times, v_mean_wt, 'purple', label='WT') 
ax1.errorbar(x=times, y=v_mean_wt, c='purple', fmt='o', yerr=v_std_wt, capsize=5, elinewidth=1, markeredgewidth=1) 

plt.legend()
plt.ylim(0,700000)
plt.xlim(-0.25,3.25)
plt.title('Chitin Trap, M9+GlcNac')
plt.xlabel('Time (d)')
plt.ylabel('Biovolume (um^3)')
plt.savefig('Figure2B.svg', format='svg')
plt.show()
#%%
"""
Figure 2C and D

Images exported from Zen
"""
#%%
"""
Figure 2E

Predation of delta rbmA grown in GlcNAc
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig2E.1')
df_pred =  pd.read_excel('SI_Data.xlsx', sheet_name='Fig2E.2') 

times= np.arange(3,7,1)
v_mean, v_std, times = average(df_pred, 'Prey', times)
p_mean, p_std, times = average(df_pred, 'Pred', times)
v_mean_d, v_std_d, times = average(df, 'Prey', times)

#effect of predation calculation used in later figure
print(np.sum(np.asarray(v_mean_d)-np.asarray(v_mean))/len(v_mean_d))
print(stats.sem(np.asarray(v_mean_d)-np.asarray(v_mean)))

fig, ax1 = plt.subplots(figsize=(3, 4))
ax2 = ax1.twinx()

ax1.plot(times, p_mean, 'gold')
ax1.errorbar(x=times, y=p_mean, c='gold', fmt='o', yerr=p_std, capsize=5, elinewidth=1, markeredgewidth=1)

ax2.plot(times, v_mean, 'purple')
ax2.errorbar(x=times, y=v_mean, c='purple', fmt='o', yerr=v_std, capsize=5, elinewidth=1, markeredgewidth=1)

ax2.plot(times, v_mean_d, 'mediumorchid', linestyle='--') 
ax2.errorbar(x=times, y=v_mean_d, c='mediumorchid', fmt='o', yerr=v_std_d, capsize=5, elinewidth=1, markeredgewidth=1) 

ax2.set_xlim(2.75, 6.25)
ax1.set_xlim(2.75, 6.25)
ax1.set_ylim(0,16000)
ax2.set_ylim(0,700000)
ax2.set_ylabel('Prey')
ax1.set_ylabel('Pred')
ax1.set_xlabel('Time d')
plt.title('GlcNAc, rbmA')
plt.savefig('Figure2E.svg', format='svg')
plt.show()
#%%
"""
Figure 2F

predation of vpsL deletion in GlcNAc
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig2F.1')
df_pred =  pd.read_excel('SI_Data.xlsx', sheet_name='Fig2F.2') 

times= np.arange(3,7,1)

v_mean, v_std, times = average(df_pred, 'Prey', times)
p_mean, p_std, times = average(df_pred, 'Pred', times)
v_mean_d, v_std_d, times = average(df, 'Prey', times)

#effect of predation calculation used in later figure
print(np.sum(np.asarray(v_mean_d)-np.asarray(v_mean))/len(v_mean_d))
print(stats.sem(np.asarray(v_mean_d)-np.asarray(v_mean)))

fig, ax1 = plt.subplots(figsize=(3, 4))
ax2 = ax1.twinx()

ax1.plot(times, p_mean, 'gold')
ax1.errorbar(x=times, y=p_mean, c='gold', fmt='o', yerr=p_std, capsize=5, elinewidth=1, markeredgewidth=1)

ax2.plot(times, v_mean, 'purple')
ax2.errorbar(x=times, y=v_mean, c='purple', fmt='o', yerr=v_std, capsize=5, elinewidth=1, markeredgewidth=1)

ax2.plot(times, v_mean_d, 'mediumorchid', linestyle='--') 
ax2.errorbar(x=times, y=v_mean_d, c='mediumorchid', fmt='o', yerr=v_std_d, capsize=5, elinewidth=1, markeredgewidth=1) 

ax2.set_xlim(2.75, 6.25)
ax1.set_xlim(2.75, 6.25)
ax1.set_ylim(0,32000)
ax2.set_ylim(0,700000)
ax2.set_ylabel('Prey')
ax1.set_ylabel('Pred')
ax1.set_xlabel('Time d')
plt.title('GlcNAc, vpsL del')
plt.savefig('Figure2F.svg', format='svg')
plt.show()
#%%
"""
Figure 2G

Phase diagram

Median density and effect of predation

See percoaltion model for dashed line value
"""
#raw data used for density calculations
#raw data used for effect of predation calculations is in each of the respective figures
chitobiose = [0.715468716, 0.812510808, 0.895028367, 0.78315755]
chitin = [0.2833054348326858, 0.4482471141513467, 0.15643408215958812, 0.17718253055396546]
rbmA = [0.733278967474043, 0.5167090702159072, 0.7959625494446243, 0.5927073377271085, 0.7130987387526068]
rbma_Chitobiose = [0.637125353354137, 0.735207754047549, 0.7284998286889839, 0.6921422569066131, 0.660676874949215, 0.7021804236493288, 0.5187146015272898]
vpsL = [0.6809196348105482, 0.59830947439377, 0.6112050721236861, 0.38178179596390716, 0.4398097589991775]
glcnac = [0.8535945964655037,0.9018552048102313,0.7962787474879456 ,0.9678083301526497, 0.8073513915576768]

#chitobiose rbmA data calc
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig2G_rbmA_Chitobiose_1')
df_pred =  pd.read_excel('SI_Data.xlsx', sheet_name='Fig2G_rbmA_Chitobiose_2')

times= np.arange(3,6,1)
v_mean, v_std, times = average(df_pred, 'Prey', times)
p_mean, p_std, times = average(df_pred, 'Pred', times)
v_mean_d, v_std_d, times = average(df, 'Prey', times)

print(np.sum(np.asarray(v_mean_d)-np.asarray(v_mean))/len(v_mean_d))
print(stats.sem(np.asarray(v_mean_d)-np.asarray(v_mean)))

df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig2G')

fig = plt.figure(figsize=(3, 3))
for i in np.arange(0,len(df)):
    plt.errorbar(x=df['Density'][i], y=df['PredEffect'][i],  ls = "None", color=df['Color'][i], yerr=df['PredEffectErr'][i], xerr=df['DensitySem'][i], alpha=0.5,  elinewidth=1)
    plt.scatter(x=df['Density'][i], y=df['PredEffect'][i], marker=df['Marker'][i], color=df['Color'][i],  s=50, alpha=0.75)
plt.ylabel('Pred Effect (um)')
plt.xlabel('Median Density')

plt.xlim(0,1)
plt.ylim(-50000,250000)
plt.axvline(x=0.7142857142857143, linestyle='--', color='blue')
plt.axhline(y=0, linestyle='--', color='red')
plt.savefig('Figure2G.svg', format='svg')
plt.show()
#%%
"""
Figure 3A

VPS Staining
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig3A')

plot_dict = {
             'Chitin_VPS': df['Chitin_VPS'],
             'GlcNAc_VPS': df['GlcNAc_VPS'] 
             }

plt.figure(figsize=(2,3))
plt.boxplot(plot_dict.values(), showfliers=False, widths=0.5)
plt.xticks(ticks=[1,2], labels=plot_dict.keys())

plt.scatter((simple_beeswarm2((plot_dict['Chitin_VPS']), nbins=2))+1, plot_dict['Chitin_VPS'], alpha=0.5, color='#59BA57', s=75, marker='p')
plt.scatter((simple_beeswarm2((plot_dict['GlcNAc_VPS'])))+2, plot_dict['GlcNAc_VPS'], alpha=0.5, color='#59BA57', s=75, marker='o')
plt.xticks(rotation=45)
plt.ylabel(r'Shell Median', fontsize=12)
plt.xlabel(r'Treatment', fontsize = 12)
plt.savefig('Figure3A.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(mannwhitneyu(plot_dict['GlcNAc_VPS'], plot_dict['Chitin_VPS']))
#%%
"""
Figure 3B

Bap1 Immunostaining
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig3B')

plot_dict = {
             'Chitin_Bap1': df['Chitin_Bap1'].dropna(),
             'GlcNAc_Bap1': df['GlcNAc_Bap1'].dropna() 
             }

plt.figure(figsize=(2,3))
plt.boxplot(plot_dict.values(), showfliers=False, widths=0.5)
plt.xticks(ticks=[1,2], labels=plot_dict.keys())

plt.scatter((simple_beeswarm2(plot_dict['Chitin_Bap1']))+1, plot_dict['Chitin_Bap1'], alpha=0.5, color='#59BA57', s=75, marker='p')
plt.scatter((simple_beeswarm2(plot_dict['GlcNAc_Bap1']))+2, plot_dict['GlcNAc_Bap1'], alpha=0.5, color='#59BA57', s=75)
plt.xticks(rotation=45)
plt.ylabel(r'Shell Median', fontsize=12)
plt.xlabel(r'Treatment', fontsize = 12)
plt.savefig('Figure3B.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(mannwhitneyu(plot_dict['GlcNAc_Bap1'], plot_dict['Chitin_Bap1']))
#%%
"""
Figure 3C

RbmC Immunostaining
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig3C')

plot_dict = {
             'Chitin_RbmC': df['Chitin_RbmC'].dropna(),
             'GlcNAc_RbmC': df['GlcNAc_RbmC'].dropna() 
             }

plt.figure(figsize=(2,3))
plt.boxplot(plot_dict.values(), showfliers=False, widths=0.5)
plt.xticks(ticks=[1,2], labels=plot_dict.keys())

plt.scatter((simple_beeswarm2(plot_dict['Chitin_RbmC']))+1, plot_dict['Chitin_RbmC'], alpha=0.5, color='#59BA57', s=75, marker='p')
plt.scatter((simple_beeswarm2(plot_dict['GlcNAc_RbmC']))+2, plot_dict['GlcNAc_RbmC'], alpha=0.5, color='#59BA57', s=75)
plt.xticks(rotation=45)
plt.ylabel(r'Shell Median', fontsize=12)
plt.xlabel(r'Treatment', fontsize = 12)
plt.savefig('Figure3C.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(mannwhitneyu(plot_dict['GlcNAc_RbmC'], plot_dict['Chitin_RbmC']))
#%%
"""
Figure 3D

RbmA Immunostaining
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig3D')

plot_dict = {
             'Chitin_RbmA': df['Chitin_RbmA'].dropna(),
             'GlcNAc_RbmA': df['GlcNAc_RbmA'].dropna() 
             }

plt.figure(figsize=(2,3))
plt.boxplot(plot_dict.values(), showfliers=False, widths=0.5)
plt.xticks(ticks=[1,2], labels=plot_dict.keys())

plt.scatter((simple_beeswarm2(plot_dict['Chitin_RbmA']))+1, plot_dict['Chitin_RbmA'], alpha=0.5, color='#59BA57', s=75)
plt.scatter((simple_beeswarm2(plot_dict['GlcNAc_RbmA']))+2, plot_dict['GlcNAc_RbmA'], alpha=0.5, color='#59BA57', s=75, marker='p')
plt.xticks(rotation=45)
plt.ylabel(r'Shell Median', fontsize=12)
plt.xlabel(r'Treatment', fontsize = 12)
plt.savefig('Figure3D.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(mannwhitneyu(plot_dict['GlcNAc_RbmA'], plot_dict['Chitin_RbmA']))
#%%
"""
Figure 3E

Biovolume of WT, vpsL, and tcpA on chitin at day 3
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig3E')

plot_dict = {'Chitin_WT': df['WT'].dropna(), 
             'Chitin_vpsL': df['vpsL'].dropna(),
             'Chitin_tcpA': df['tcpA'].dropna(),
             }

plt.figure(figsize=(3,3))
plt.boxplot(plot_dict.values(), showfliers=False, widths=0.5)
plt.xticks(ticks=[1,2,3], labels=plot_dict.keys(), fontsize=12)

plt.scatter((simple_beeswarm2(plot_dict['Chitin_WT']))+1, plot_dict['Chitin_WT'], alpha=0.5, color='purple', s=75, marker='p')
plt.scatter((simple_beeswarm2(plot_dict['Chitin_vpsL']))+2, plot_dict['Chitin_vpsL'], alpha=0.5, color='red', s=75, marker='p')
plt.scatter((simple_beeswarm2(plot_dict['Chitin_tcpA']))+3, plot_dict['Chitin_tcpA'], alpha=0.5, color='blue', s=75, marker='p')
plt.xticks(rotation=45)
plt.ylim(0,212*212*15/2)
plt.ylabel(r'Biovolume', fontsize=12)
plt.xlabel(r'Treatment', fontsize = 12)
plt.title('Day 3', fontsize = 12)
plt.savefig('Figure3E.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(mannwhitneyu(plot_dict['Chitin_WT'], plot_dict['Chitin_vpsL']))
print(mannwhitneyu(plot_dict['Chitin_WT'], plot_dict['Chitin_tcpA']))
#%%
"""
Figure 3F

Immunostaining of Rugose and WT grown on chitin
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig3F')

plot_dict = {'GlcNAc_WT': df['WT'] ,
             'GlcNAc_tcpA': df['tcpA']
             }

plt.figure(figsize=(2,3))
plt.boxplot(plot_dict.values(), showfliers=False, widths=0.5)
plt.xticks(ticks=[1,2], labels=plot_dict.keys(), fontsize=12)

plt.scatter((simple_beeswarm2(plot_dict['GlcNAc_WT'], nbins=1))+1, plot_dict['GlcNAc_WT'], alpha=0.5, color='purple', s=75, marker='o')
plt.scatter((simple_beeswarm2(plot_dict['GlcNAc_tcpA']))+2, plot_dict['GlcNAc_tcpA'], alpha=0.5, color='blue', s=75, marker='o')
plt.xticks(rotation=45)
plt.ylim(0,212*212*15)
plt.ylabel(r'Biovolume', fontsize=12)
plt.xlabel(r'Treatment', fontsize = 12)
plt.title('Day 3', fontsize = 12)
plt.savefig('Figure3F.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(mannwhitneyu(plot_dict['GlcNAc_WT'], plot_dict['GlcNAc_tcpA']))
#%%
"""
Figure 4D

Immunostaining of Rugose and WT grown on chitin
"""
df_new = pd.read_excel('SI_Data.xlsx', sheet_name='Fig4D')
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig3D')

plot_dict = {
             'Chitin_RbmA': df['Chitin_RbmA'].dropna(),
             'Chitin_r_RbmA': df_new['Chitin_Rugose'].dropna(),
             'GlcNAc_RbmA': df['GlcNAc_RbmA'].dropna()
             }

plt.figure(figsize=(2.5,3.5))
plt.boxplot(plot_dict.values(), showfliers=False, widths=0.5)
plt.xticks(ticks=[1,2,3], labels=plot_dict.keys())

plt.scatter((simple_beeswarm2((plot_dict['Chitin_RbmA']), nbins=2))+1, plot_dict['Chitin_RbmA'], alpha=0.5, color='#59BA57', s=75, marker='p')
plt.scatter((simple_beeswarm2((plot_dict['Chitin_r_RbmA'])))+2, plot_dict['Chitin_r_RbmA'], alpha=0.5, color='#59BA57', s=75, marker='p')
plt.scatter((simple_beeswarm2((plot_dict['GlcNAc_RbmA'])))+3, plot_dict['GlcNAc_RbmA'], alpha=0.5, color='#59BA57', s=75, marker='o')
plt.xticks(rotation=45)
plt.ylabel(r'Shell Median', fontsize=12)
plt.xlabel(r'Treatment', fontsize = 12)
plt.title('Immunostaining')
plt.savefig('Figure4D.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(mannwhitneyu(plot_dict['Chitin_RbmA'], plot_dict['Chitin_r_RbmA']))
#%%
"""
Figure 4E

Rugose, WT, and Rugose-vpsL biofilm local density on chitin
"""
#load data
f1 = pd.read_excel('SI_Data.xlsx', sheet_name='Fig4E_WT')
f2 = pd.read_excel('SI_Data.xlsx', sheet_name='Fig4E_Rugose_vpsL')
f3 = pd.read_excel('SI_Data.xlsx', sheet_name='Fig4E_Rugose')

#%%
sns.histplot(data=f1, x='Density', bins=50, color='#70CDDD', alpha=0.5, stat='probability', label ='Chitin, WT', weights='Vol')
sns.histplot(data=f2, x='Density', bins=50, color='black', alpha=0.5, stat='probability', label ='Chitin, Rugose, del vpsL', weights='Vol')
sns.histplot(data=f3, x='Density', bins=50, color='purple', alpha=0.5, stat='probability', label ='Chitin, Rugose', weights='Vol')
    
plt.ylabel('Volume Weighted Frequency')
plt.xlabel('Cell Packing (6um)')
# plt.ylim(0,0.1)
plt.xlim(0,1)
plt.legend()
plt.savefig('Figure4E.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
Figure 4H

Predation of rugose strain grown on chitin substrate
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig4H.1')
df_pred =  pd.read_excel('SI_Data.xlsx', sheet_name='Fig4H.2') 

times= np.arange(3,7,1)

v_mean, v_std, times = average(df_pred, 'Prey', times)
p_mean, p_std, times = average(df_pred, 'Pred', times)
v_mean_d, v_std_d, times = average(df, 'Prey', times)

print(np.sum(np.asarray(v_mean_d)-np.asarray(v_mean))/len(v_mean_d))
print(stats.sem(np.asarray(v_mean_d)-np.asarray(v_mean)))

fig, ax1 = plt.subplots(figsize=(4, 4))
ax2 = ax1.twinx()

ax1.plot(times, p_mean, 'gold')
ax1.errorbar(x=times, y=p_mean, c='gold', fmt='p', yerr=p_std, capsize=5, elinewidth=1, markeredgewidth=1)

ax2.plot(times, v_mean, 'purple')
ax2.errorbar(x=times, y=v_mean, c='purple', fmt='p', yerr=v_std, capsize=5, elinewidth=1, markeredgewidth=1)

ax2.plot(times, v_mean_d, 'mediumorchid', linestyle='--') 
ax2.errorbar(x=times, y=v_mean_d, c='mediumorchid', fmt='p', yerr=v_std_d, capsize=5, elinewidth=1, markeredgewidth=1) 

ax2.set_xlim(2.75, 6.25)
ax1.set_xlim(2.75, 6.25)
ax1.set_ylim(0,16000)
ax2.set_ylim(0,500000)
ax2.set_ylabel('Prey')
ax1.set_ylabel('Pred')
ax1.set_xlabel('Time d')
plt.title('Chitin, Rugose')
plt.savefig('Figure4H.svg', format='svg')
plt.show()
#%%
"""
Figure 4I
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig4I.1')
df_pred =  pd.read_excel('SI_Data.xlsx', sheet_name='Fig4I.2') 

times= np.arange(3,7,1)

v_mean, v_std, times = average(df_pred, 'Prey', times)
p_mean, p_std, times = average(df_pred, 'Pred', times)
v_mean_d, v_std_d, times = average(df, 'Prey', times)

print(np.sum(np.asarray(v_mean_d)-np.asarray(v_mean))/len(v_mean_d))
print(stats.sem(np.asarray(v_mean_d)-np.asarray(v_mean)))

fig, ax1 = plt.subplots(figsize=(4, 4))
ax2 = ax1.twinx()

ax1.plot(times, p_mean, 'gold')
ax1.errorbar(x=times, y=p_mean, c='gold', fmt='p', yerr=p_std, capsize=5, elinewidth=1, markeredgewidth=1)

ax2.plot(times, v_mean, 'purple')
ax2.errorbar(x=times, y=v_mean, c='purple', fmt='p', yerr=v_std, capsize=5, elinewidth=1, markeredgewidth=1)

ax2.plot(times, v_mean_d, 'mediumorchid', linestyle='--') 
ax2.errorbar(x=times, y=v_mean_d, c='mediumorchid', fmt='p', yerr=v_std_d, capsize=5, elinewidth=1, markeredgewidth=1) 

ax2.set_xlim(2.75, 6.25)
ax1.set_xlim(2.75, 6.25)
ax1.set_ylim(0,16000)
ax2.set_ylim(0,500000)
ax2.set_ylabel('Prey')
ax1.set_ylabel('Pred')
ax1.set_xlabel('Time d')
plt.title('Chitin, Rugose, del vpsL')
plt.savefig('Figure4I.svg', format='svg')
plt.show()
#%%
"""
SI Figures
"""
#%%
"""
SI Figure 2A
Related to Figure 1

Box and whisker plot of predation on chitin
"""

df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1A.1')
df_pred =  pd.read_excel('SI_Data.xlsx', sheet_name='Fig1A.2') 

plot_dict = {'2hPI P-':list(df.drop(df[df.Day > 3].index)['Prey']), 
             '2hPI P+':list(df_pred.drop(df_pred[df_pred.Day > 3].index)['Prey']),
             '72hPI P-':list(df.drop(df[df.Day != 5].index)['Prey']), 
              '72hPI P+':list(df_pred.drop(df_pred[df_pred.Day != 5].index)['Prey']),}

plt.figure(figsize=(3,4))
plt.boxplot(plot_dict.values(), showfliers=False)
plt.xticks(ticks=[1,2, 3, 4], labels=plot_dict.keys(), fontsize=20)

plt.scatter(simple_beeswarm2(plot_dict['2hPI P-'])+1, plot_dict['2hPI P-'], alpha=0.5, color='black', s=75, marker='p')
plt.scatter(simple_beeswarm2(plot_dict['2hPI P+'], nbins=2)+2, plot_dict['2hPI P+'], alpha=0.5, color='red', s=75, marker='p')
plt.scatter(simple_beeswarm2(plot_dict['72hPI P-'])+3, plot_dict['72hPI P-'], alpha=0.5, color='black', s=75, marker='p')
plt.scatter(simple_beeswarm2(plot_dict['72hPI P+'])+4, plot_dict['72hPI P+'], alpha=0.5, color='red', s=75, marker='p')
plt.xticks(rotation=45)
plt.hlines(y=1, xmin=0.5, xmax=4.5, linewidth=1, linestyles='--', color='r')
plt.ylabel(r'Biovolume', fontsize=20)
plt.xlabel(r'Treatment', fontsize = 20)
plt.ylim(-50000,600000)
plt.title('WT', fontsize = 20)
plt.savefig('FigureS2A.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(mannwhitneyu(plot_dict['72hPI P-'], plot_dict['72hPI P+']))
print(mannwhitneyu(plot_dict['2hPI P-'], plot_dict['2hPI P+']))
#%%
"""
SI Figure S2B
Related to Figure 1

Box and whisker plot of predation on chitobiose
"""

plt.figure(figsize=(3,4))

df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1B.1')
df_pred =  pd.read_excel('SI_Data.xlsx', sheet_name='Fig1B.2') 

plot_dict = {'2hPI P-':list(df.drop(df[df.Day > 3].index)['Prey']), 
             '2hPI P+':list(df_pred.drop(df_pred[df_pred.Day > 3].index)['Prey']),
             '48hPI P-':list(df.drop(df[df.Day != 5].index)['Prey']), 
              '48hPI P+':list(df_pred.drop(df_pred[df_pred.Day != 5].index)['Prey']),}

plt.boxplot(plot_dict.values(), showfliers=False)
plt.xticks(ticks=[1,2, 3, 4], labels=plot_dict.keys(), fontsize=20)

plt.scatter((simple_beeswarm2(plot_dict['2hPI P-']))+1, plot_dict['2hPI P-'], alpha=0.5, color='black', s=75, marker='s')
plt.scatter((simple_beeswarm2(plot_dict['2hPI P+']))+2, plot_dict['2hPI P+'], alpha=0.5, color='red', s=75, marker='s')
plt.scatter((simple_beeswarm2(plot_dict['48hPI P-']))+3, plot_dict['48hPI P-'], alpha=0.5, color='black', s=75, marker='s')
plt.scatter((simple_beeswarm2(plot_dict['48hPI P+']))+4, plot_dict['48hPI P+'], alpha=0.5, color='red', s=75, marker='s')

plt.xticks(rotation=45)
plt.ylim(-50000, 700000)
plt.hlines(y=1, xmin=0.5, xmax=4.5, linewidth=1, linestyles='--', color='r')
plt.ylabel(r'Biovolume', fontsize=20)
plt.xlabel(r'Treatment', fontsize = 20)
plt.title('WT', fontsize = 20)
plt.savefig('FigureS2B.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(mannwhitneyu(plot_dict['48hPI P-'], plot_dict['48hPI P+']))
print(mannwhitneyu(plot_dict['2hPI P-'], plot_dict['2hPI P+']))
#%%
"""
SI Figure S2C
Related to Figure 1

Box and whisker plot of predation on GlcNAc
"""

df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1C.1')
df_pred =  pd.read_excel('SI_Data.xlsx', sheet_name='Fig1C.2') 

plot_dict = {'2hPI P-':list(df.drop(df[df.Day > 3].index)['Prey']), 
             '2hPI P+':list(df_pred.drop(df_pred[df_pred.Day > 3].index)['Prey']),
             '48hPI P-':list(df.drop(df[df.Day != 5].index)['Prey']), 
              '48hPI P+':list(df_pred.drop(df_pred[df_pred.Day != 5].index)['Prey']),}

plt.figure(figsize=(3,4))
plt.boxplot(plot_dict.values(), showfliers=False)
plt.xticks(ticks=[1,2, 3, 4], labels=plot_dict.keys(), fontsize=20)

plt.scatter((simple_beeswarm2(plot_dict['2hPI P-']))+1, plot_dict['2hPI P-'], alpha=0.5, color='black', s=75)
plt.scatter((simple_beeswarm2(plot_dict['2hPI P+']))+2, plot_dict['2hPI P+'], alpha=0.5, color='red', s=75)
plt.scatter((simple_beeswarm2(plot_dict['48hPI P-']))+3, plot_dict['48hPI P-'], alpha=0.5, color='black', s=75)
plt.scatter((simple_beeswarm2(plot_dict['48hPI P+']))+4, plot_dict['48hPI P+'], alpha=0.5, color='red', s=75)

plt.xticks(rotation=45)
plt.ylim(-50000, 700000)
plt.hlines(y=1, xmin=0.5, xmax=4.5, linewidth=1, linestyles='--', color='r')
plt.ylabel(r'Biovolume', fontsize=20)
plt.xlabel(r'Treatment', fontsize = 20)
plt.title('WT', fontsize = 20)
plt.savefig('FigureS2C.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(mannwhitneyu(plot_dict['48hPI P-'], plot_dict['48hPI P+']))
print(mannwhitneyu(plot_dict['2hPI P-'], plot_dict['2hPI P+']))
#%%
"""
SI Figure S2D
Related to Figure 1

Chitin loss curve 
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='SI_Figure_S2D')
plt.figure(figsize=(3,4))
plt.step(df['Day'], df['Pred'], color='red', label='Pred+')
plt.step(df['Day'], df['No Pred'], color='black', label='Pred-', linestyle='--')
plt.legend()
plt.ylim(0,1.1)
plt.xlabel('Time (d)')
plt.ylabel('Fraction Chitin Remaining')
plt.savefig('FigureS2D.svg', format='svg')
plt.show()

#%%
"""
SI Figure S3A
Related to Figure 1

Box and whisker plot comparing biovolume accumulation at time of predator addition
"""

df_chitin = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1A.1')
df_glcnac2 = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1B.1')
df_glcnac = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1C.1')

plot_dict = {'Chitin':list(df_chitin.drop(df_chitin[df_chitin.Day > 3].index)['Prey']), 
             'GlcNAc2':list(df_glcnac2.drop(df_glcnac2[df_glcnac2.Day > 3].index)['Prey']),
             'GlcNAc':list(df_glcnac.drop(df_glcnac[df_glcnac.Day > 3].index)['Prey'])}

plt.figure(figsize=(3,4))
plt.boxplot(plot_dict.values(), showfliers=False)
plt.xticks(ticks=[1,2, 3], labels=plot_dict.keys(), fontsize=20)

plt.scatter((simple_beeswarm2(plot_dict['Chitin']))+1, plot_dict['Chitin'], alpha=0.5, color='cyan', marker='p', s=75)
plt.scatter((simple_beeswarm2(plot_dict['GlcNAc2']))+2, plot_dict['GlcNAc2'], alpha=0.5, color='orange', marker='s', s=75)
plt.scatter((simple_beeswarm2(plot_dict['GlcNAc']))+3, plot_dict['GlcNAc'], alpha=0.5, color='purple', s=75)

plt.xticks(rotation=45)
plt.ylim(-50000, 700000)
plt.ylabel(r'Biovolume', fontsize=20)
plt.xlabel(r'Treatment', fontsize = 20)
plt.title('WT', fontsize = 20)
plt.savefig('FigureS3A.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(mannwhitneyu(plot_dict['Chitin'], plot_dict['GlcNAc']))
print(mannwhitneyu(plot_dict['Chitin'], plot_dict['GlcNAc2']))
#%%
"""
SI Figure S3B

Comparison of predator dynamics
"""
df_pred =  pd.read_excel('SI_Data.xlsx', sheet_name='Fig1A.2') 
plt.figure(figsize=(3,4))

times= np.arange(3,7,1)
p_mean, p_std, times = average(df_pred, 'Pred', times)

plt.plot(times, p_mean, '#6CCBD8')
plt.errorbar(x=times, y=p_mean, c='#6CCBD8', fmt='p', yerr=p_std, capsize=5, elinewidth=1, markeredgewidth=1)

df_pred =  pd.read_excel('SI_Data.xlsx', sheet_name='Fig1B.2') 

times= np.arange(3,6,1)
p_mean, p_std, times = average(df_pred, 'Pred', times)
plt.plot(times, p_mean, '#FAA51A')
plt.errorbar(x=times, y=p_mean, c='#FAA51A', fmt='s', yerr=p_std, capsize=5, elinewidth=1, markeredgewidth=1)

df_pred =  pd.read_excel('SI_Data.xlsx', sheet_name='Fig1C.2') 

times= np.arange(3,6,1)
p_mean, p_std, times = average(df_pred, 'Pred', times)
plt.plot(times, p_mean, '#7D2880')
plt.errorbar(x=times, y=p_mean, c='#7D2880', fmt='o', yerr=p_std, capsize=5, elinewidth=1, markeredgewidth=1)

plt.ylabel('Bdello biovolume')
plt.xlabel('Time (d)')
plt.savefig('FigureS3B.svg', format='svg')
plt.show()
#%%
"""
SI Figure S3C 
Related to Figure 1

Cell packing at day 5 with no predation
"""
df1 = pd.read_excel("SI_Data.xlsx", sheet_name='SI_Figure_S3C_1')

df2 = pd.read_excel("SI_Data.xlsx", sheet_name='SI_Figure_S3C_2')

df3 = pd.read_excel("SI_Data.xlsx", sheet_name='SI_Figure_S3C_3')

plt.figure(figsize=(3,4))
sns.histplot(data=df1, x='Density', bins=50, color='#6CCBD8', alpha=0.5, stat='probability', label ='Chitin', cumulative=False, weights='Vol')
sns.histplot(data=df2, x='Density', bins=50, color='#FAA51A', alpha=0.5, stat='probability', label ='GlcNAc2', cumulative=False, weights='Vol')
sns.histplot(data=df3, x='Density', bins=50, color='#7D2880', alpha=0.5, stat='probability', label ='GlcNAc, WT', cumulative=False, weights='Vol')
plt.title('Day 5, P-')
plt.legend()
plt.ylim(0,0.3)
plt.xlim(0,1)
plt.xlabel('Cell packing')
plt.ylabel('Volume Weighted Frequency')
plt.savefig('FigureS3C.svg', format='svg')
plt.show()
#%%
"""
SI Figure S4D

Related to Figure 1

See percolation model scripts for simulation data.
"""

df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1G')
df = df.fillna(0)
bins = np.arange(0,1,0.02)

y_data = []
bins_sample = []

counter = 0
for i in np.arange(0,len(bins)-1):
    data2 = df.drop(df[(df['Density'] != bins[i])].index, inplace=False)
    if len(data2) > 1:
        y_data.append(np.mean(data2['Predation']))
        bins_sample.append(bins[i])
    if counter == 0:
        if np.median(data2['Predation']) == 0:
            counter+=1
            vline = bins[i]

plt.figure(figsize=(3,3))
plt.plot(np.log10(np.asarray(bins_sample[:])),np.log10(y_data), color='black')
plt.grid(None)
plt.ylabel('Predation Probability')
plt.xlabel('Cell Packing')
plt.xlim(-1,0)
plt.ylim(-1,0.1)
plt.grid(None)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.savefig('SI_Figure_S4D.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S5A
"""
#load data from excel
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig2B_Chitin') 
df['Biovolume'] = np.log(df['Biovolume'])

data_reg = df.drop(df[(df.Day == 0) | (df.Day > 3)].index, inplace = False)
slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(data_reg['Day'], data_reg['Biovolume'])
sns.regplot(data=data_reg, x='Day', y='Biovolume', marker='p', line_kws={'label':"y={0:.2f}x+{1:0.1f}".format(np.round(slope1, decimals=2), np.round(intercept1, decimals=0))}, scatter=False, color='orange')

#r for chitin
print(slope1, intercept1, r_value1, p_value1, std_err1)

plt.scatter(df['Day'], df['Biovolume'], marker='p', alpha=0.5, color='purple', s=80)
plt.ylabel('ln(Biovolume)')
plt.title('Chitin, WT')
plt.xlim(-0.25,3.25)
plt.ylim(5,15)

plt.savefig('SI_FigureS5A_Chitin.svg', format='svg')
plt.show()
#%%
#load data from excel
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig2B_GlcNAc') 
df['Biovolume'] = np.log(df['Biovolume'])

data_reg = df.drop(df[(df.Day > 2)].index, inplace = False)
slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(data_reg['Day'], data_reg['Biovolume'])
sns.regplot(data=data_reg, x='Day', y='Biovolume', marker='p', line_kws={'label':"y={0:.2f}x+{1:0.1f}".format(np.round(slope1, decimals=2), np.round(intercept1, decimals=0))}, scatter=False, color='orange')

#r for glcnac
print(slope1, intercept1, r_value1, p_value1, std_err1)

plt.scatter(df['Day'], df['Biovolume'], marker='o', alpha=0.5, color='purple', s=80)
plt.ylabel('ln(Biovolume)')
plt.title('GlcNAc, WT')
plt.xlim(-0.25,3.25)
plt.ylim(5,15)

plt.savefig('SI_FigureS5A_GlcNAc.svg', format='svg')
plt.show()
#%%
"""
SI Figure S5B
Related to Figure 2

carrying capacity of WT
"""
df_chitin = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1A.1')

df_glcnac = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1C.1')

df_chitobiose = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1B.1')

df_glcnac.drop(df_glcnac[(df_glcnac.Day != 5)].index, inplace = True)
df_chitin.drop(df_chitin[(df_chitin.Day != 5)].index, inplace = True)
df_chitobiose.drop(df_chitobiose[(df_chitobiose.Day != 5)].index, inplace = True)

plot_dict = {'5d_chitin':df_chitin['Prey'], '5d_glcnac':df_glcnac['Prey'], '5d_chitobiose':df_chitobiose['Prey']}

plt.figure(figsize=(3,3))
plt.boxplot(plot_dict.values(), showfliers=False, widths=0.5)
plt.xticks(ticks=[1,2,3], labels=plot_dict.keys(), fontsize=20)

plt.scatter((simple_beeswarm2(plot_dict['5d_chitin']))+1, plot_dict['5d_chitin'], alpha=0.5, color='purple', s=75, marker='p')
plt.scatter((simple_beeswarm2(plot_dict['5d_glcnac']))+2, plot_dict['5d_glcnac'], alpha=0.5, color='purple', s=75, marker='o')
plt.scatter((simple_beeswarm2(plot_dict['5d_chitobiose']))+3, plot_dict['5d_chitobiose'], alpha=0.5, color='purple', s=75, marker='s')
plt.xticks(rotation=45)

plt.xlabel(r'Carbon Source', fontsize = 20)
plt.ylim(0,700000)

plt.title('Biovolume', fontsize = 20)
plt.savefig('SI_FigureS5B.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(np.median(plot_dict['5d_glcnac']))
print(np.median(plot_dict['5d_chitin']))
print(np.median(plot_dict['5d_chitobiose']))
#%%
"""
SI Figure S5C
Related to Figure 2

Death curve of Bdello monoculture 
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='SI_Figure_S5C') 

# Parameters
D = 1.7952   # dilution rate, 1/h
x0 = 15000*2 # initial concentration of Bdello
t_span = (0, 24)  # time (hours)
t_eval = np.linspace(t_span[0], t_span[1], 200)

# ODE definition
def dxdt(t, x, D):
    return [-D * x[0]]

# Solve
sol = solve_ivp(dxdt, t_span, [x0], args=(D,), t_eval=t_eval)

# Plot
fig, ax1 = plt.subplots(figsize=(3, 3))
plt.plot(sol.t, sol.y[0], label="Dilution", lw=2, linestyle='--')
plt.scatter(df['Time'], df['Bdello'], c='gold', alpha=0.75, marker='p')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2))
plt.xticks(np.arange(0,25,12))
plt.xlim(-1,25)
ax1.set_ylabel('Bdello')
ax1.set_xlabel('Time (h)')
plt.savefig('SI_FigureS5C.svg', format='svg')
plt.show()
#%%
"""
SI Figure S6
Related to Figure 1

See model scripts. 
"""
#%%
"""
SI Figure S7A
Related to Figure 3

Box and whisker plot of delta rbmA predation in GlcNAc
"""

df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig2E.1')
df_pred =  pd.read_excel('SI_Data.xlsx', sheet_name='Fig2E.2') 

plot_dict = {'2hPI P-':list(df.drop(df[df.Day > 3].index)['Prey']), 
             '2hPI P+':list(df_pred.drop(df_pred[df_pred.Day > 3].index)['Prey']),
             '48hPI P-':list(df.drop(df[df.Day != 5].index)['Prey']), 
              '48hPI P+':list(df_pred.drop(df_pred[df_pred.Day != 5].index)['Prey']),}

plt.figure(figsize=(3,4))
plt.boxplot(plot_dict.values(), showfliers=False)
plt.xticks(ticks=[1,2, 3, 4], labels=plot_dict.keys(), fontsize=20)

plt.scatter((simple_beeswarm2(plot_dict['2hPI P-']))+1, plot_dict['2hPI P-'], alpha=0.5, color='black', s=75)
plt.scatter((simple_beeswarm2(plot_dict['2hPI P+']))+2, plot_dict['2hPI P+'], alpha=0.5, color='red', s=75)
plt.scatter((simple_beeswarm2(plot_dict['48hPI P-']))+3, plot_dict['48hPI P-'], alpha=0.5, color='black', s=75)
plt.scatter((simple_beeswarm2(plot_dict['48hPI P+']))+4, plot_dict['48hPI P+'], alpha=0.5, color='red', s=75)
plt.xticks(rotation=45)
plt.hlines(y=1, xmin=0.5, xmax=4.5, linewidth=1, linestyles='--', color='r')
plt.ylabel(r'Biovolume', fontsize=20)
plt.xlabel(r'Treatment', fontsize = 20)
plt.title('rbmA del', fontsize = 20)
plt.ylim(-50000, 700000)

plt.savefig('FigureS7A.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(mannwhitneyu(plot_dict['48hPI P-'], plot_dict['48hPI P+']))
print(mannwhitneyu(plot_dict['2hPI P-'], plot_dict['2hPI P+']))
#%%
"""
SI Figure S7B
Related to Figure 3

Box and whisker plot of delta vpsL predation in GlcNAc
"""

df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig2F.1')
df_pred =  pd.read_excel('SI_Data.xlsx', sheet_name='Fig2F.2') 

plot_dict = {'2hPI P-':list(df.drop(df[df.Day > 3].index)['Prey']), 
             '2hPI P+':list(df_pred.drop(df_pred[df_pred.Day > 3].index)['Prey']),
             '48hPI P-':list(df.drop(df[df.Day != 5].index)['Prey']), 
              '48hPI P+':list(df_pred.drop(df_pred[df_pred.Day != 5].index)['Prey']),}

plt.figure(figsize=(3,4))
plt.boxplot(plot_dict.values(), showfliers=False)
plt.xticks(ticks=[1,2, 3, 4], labels=plot_dict.keys(), fontsize=20)
plt.scatter((simple_beeswarm2(plot_dict['2hPI P-']))+1, plot_dict['2hPI P-'], alpha=0.5, color='black', s=75)
plt.scatter((simple_beeswarm2(plot_dict['2hPI P+']))+2, plot_dict['2hPI P+'], alpha=0.5, color='red', s=75)
plt.scatter((simple_beeswarm2(plot_dict['48hPI P-']))+3, plot_dict['48hPI P-'], alpha=0.5, color='black', s=75)
plt.scatter((simple_beeswarm2(plot_dict['48hPI P+']))+4, plot_dict['48hPI P+'], alpha=0.5, color='red', s=75)
plt.xticks(rotation=45)
plt.hlines(y=1, xmin=0.5, xmax=4.5, linewidth=1, linestyles='--', color='r')
plt.ylabel(r'Biovolume', fontsize=20)
plt.xlabel(r'Treatment', fontsize = 20)
plt.title('vpsL del', fontsize = 20)
plt.ylim(-50000, 700000)
plt.savefig('FigureS7B.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(mannwhitneyu(plot_dict['48hPI P-'], plot_dict['48hPI P+']))
print(mannwhitneyu(plot_dict['2hPI P-'], plot_dict['2hPI P+']))

#%%
"""
SI Figure S7E
Related to Figure 3

Box and whisker plot of delta vpsL and WT bdello plaques from a 10**9 stock of bdello
"""

plot_dict = {'WT':[142e8/100,42e8/100,192e7/100],
             'vpsL':[140e7/100, 146e7/100, 400e8/100]}

plt.figure(figsize=(3,4))
plt.boxplot(plot_dict.values(), showfliers=False)
plt.xticks(ticks=[1,2], labels=plot_dict.keys(), fontsize=20)

plt.scatter((simple_beeswarm2(plot_dict['WT'], nbins=3))+1, plot_dict['WT'], alpha=0.5, color='purple', s=75)
plt.scatter((simple_beeswarm2(plot_dict['vpsL'], nbins=3))+2, plot_dict['vpsL'], alpha=0.5, color='firebrick', s=75)

plt.xticks(rotation=45)
plt.ylim(1,1e9)
plt.ylabel(r'PFU/uL', fontsize=12)
plt.title('Plaque Assay', fontsize = 20)
plt.yscale("log")
plt.savefig('FigureS7E.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(mannwhitneyu(plot_dict['WT'], plot_dict['vpsL']))
#%%
"""
SI Figure S8
Related to Figure 2

Images of vpsL deletion mutant grown in different carbon sources
"""
#%%
"""
SI Figure S9
Related to Figure 3

Images of VPS, Bap1, RbmC, and RbmA stain. 
"""
#%%
"""
SI Figure S10C
Related to Figure 3

Dual carbon source matrix staining experiment 

Control is the same expeirmental setup minus the matrix stain
"""


df =  pd.read_excel('SI_Data.xlsx', sheet_name='SI_Figure_S10C')

plot_dict = {'Chitin':list(df['Chitin'].dropna()), 
             'Control Chitin':list(df['Chitin_Control'].dropna()),
             'Microcolony':list(df['Microcolony'].dropna()), 
             'Control Micro':list(df['Microcolony_Control'].dropna()),}

plt.figure(figsize=(3,3))
plt.boxplot(plot_dict.values(), showfliers=False)
plt.xticks(ticks=[1,2,3,4], labels=plot_dict.keys(), fontsize=20)

plt.scatter((simple_beeswarm2(plot_dict['Chitin'], nbins=5))+1, plot_dict['Chitin'], alpha=0.5, color='#59BA57', s=75, marker='p')
plt.scatter((simple_beeswarm2(plot_dict['Control Chitin'], nbins=5))+2, plot_dict['Control Chitin'], alpha=0.5, color='grey', s=75, marker='p')
plt.scatter((simple_beeswarm2(plot_dict['Microcolony'], nbins=5))+3, plot_dict['Microcolony'], alpha=0.5, color='#59BA57', s=75, marker='o')
plt.scatter((simple_beeswarm2(plot_dict['Control Micro'], nbins=5))+4, plot_dict['Control Micro'], alpha=0.5, color='grey', s=75, marker='o')

plt.xticks(rotation=45)
plt.ylim(0,450)
plt.ylabel(r'Intensity (A.U.)', fontsize=12)
plt.title('Dual Carbon Source Immunostaining', fontsize = 20)
plt.savefig('FigureS10C.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(mannwhitneyu(plot_dict['Chitin'], plot_dict['Control Chitin']))
print(mannwhitneyu(plot_dict['Microcolony'], plot_dict['Control Micro']))
print(mannwhitneyu(plot_dict['Chitin'], plot_dict['Microcolony']))

#%%
"""
SI Figure S11A
Related to Figure 3

Competition between vpsL deletion and WT in GlcNAc
"""

df = pd.read_excel('SI_Data.xlsx', sheet_name='SI_Figure_S11B')
x = df['RA_t0']                 
y = df['RA_Change']

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(x, y)
                
slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(x, y)
print(r_value1**2)
print(p_value1)

plt.figure(figsize=(6,3))
sns.regplot(x=x, y=y, color='red', ci=None, scatter=False)
plt.scatter(x=x, y=y, color='red', s=200, alpha=0.5)
plt.plot(np.arange(0,1.1,0.1), 1-np.arange(0,1.1,0.1), color='blue', linestyle='--')
plt.plot(np.arange(0,1.1,0.1), -np.arange(0,1.1,0.1), color='blue', linestyle='--')
plt.axhline(y=0, color='r', linestyle='--')
plt.ylim(-1, 1)
plt.xlim(0,1)

plt.xlabel('vpsL RA 0 d', fontsize=12)
plt.ylabel('vpsL RA 5 d', fontsize=12)
plt.savefig('FigureS11A.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S11B
Related to Figure 3

Competition between vpsL deletion and WT on Chitin
"""

df = pd.read_excel('SI_Data.xlsx', sheet_name='SI_Figure_S11A')
x = df['RA_t0']                 
y = df['RA_Change']

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(x, y)

plt.figure(figsize=(6,3))
sns.regplot(x=x, y=y, color='red', ci=None, scatter=False)
plt.scatter(x=x, y=y, color='red', marker='p', s=200, alpha=0.5)
plt.plot(np.arange(0,1.1,0.1), 1-np.arange(0,1.1,0.1), color='blue', linestyle='--')
plt.plot(np.arange(0,1.1,0.1), -np.arange(0,1.1,0.1), color='blue', linestyle='--')
plt.axhline(y=0, color='r', linestyle='--')
plt.ylim(-1, 1)
plt.xlim(0,1)

plt.xlabel('vpsL RA 0 d', fontsize=12)
plt.ylabel('vpsL RA 5 d', fontsize=12)
plt.savefig('FigureS11B.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S11C
Related to Figure 3

Competition between tcpA deletion and WT in GlcNAc
"""

df = pd.read_excel('SI_Data.xlsx', sheet_name='SI_Figure_S11C')
x = df['RA_t0']                 
y = df['RA_Change']
                
slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(x, y)
print(r_value1**2)
print(p_value1)

plt.figure(figsize=(6,3))
sns.regplot(x=x, y=y, color='blue', ci=None, scatter=False)
plt.scatter(x=x, y=y, color='blue', marker='o', s=200, alpha=0.5)
plt.plot(np.arange(0,1.1,0.1), 1-np.arange(0,1.1,0.1), color='blue', linestyle='--')
plt.plot(np.arange(0,1.1,0.1), -np.arange(0,1.1,0.1), color='blue', linestyle='--')
plt.axhline(y=0, color='r', linestyle='--')
plt.ylim(-1, 1)
plt.xlim(0,1)

plt.xlabel('tcpA RA 0 d', fontsize=12)
plt.ylabel('tcpA RA 5 d', fontsize=12)
plt.savefig('FigureS11C.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S11D
Related to Figure 3

Competition between tcpA deletion and WT on Chitin
"""

df = pd.read_excel('SI_Data.xlsx', sheet_name='SI_Figure_S11D')
x = df['RA_t0']                 
y = df['RA_Change']

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(x, y)
print(r_value1**2)
print(p_value1)

plt.figure(figsize=(6,3))
sns.regplot(x=x, y=y, color='blue', ci=None, scatter=False)
plt.scatter(x=x, y=y, color='blue', marker='p', s=200, alpha=0.5)
plt.plot(np.arange(0,1.1,0.1), 1-np.arange(0,1.1,0.1), color='blue', linestyle='--')
plt.plot(np.arange(0,1.1,0.1), -np.arange(0,1.1,0.1), color='blue', linestyle='--')
plt.axhline(y=0, color='r', linestyle='--')
plt.ylim(-1, 1)
plt.xlim(0,1)

plt.xlabel('tcpA RA 0 d', fontsize=12)
plt.ylabel('tcpA RA 5 d', fontsize=12)
plt.savefig('FigureS11D.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S12D
"""
data = pd.read_excel('SI_Data.xlsx', sheet_name='SI_Figure_S12D')['Pearson_r']

plt.figure(figsize=(2,3))
plt.boxplot(data, showfliers=False, widths=0.5)
plt.scatter(np.ones(len(data)), data)
plt.ylabel('Pearson r')
plt.xlabel('GlcNAc Z-Stacks')
plt.ylim(0,1)
plt.savefig('SI_Figure_S12D.svg', format='svg')
plt.show()
#%%
"""
SI Figure S13A
Related to Figure 1

Box and whisker plot of predation on GlcNAc
"""

df_R = pd.read_excel('SI_Data.xlsx', sheet_name='Fig4H.1')
df_R_vpsL = pd.read_excel('SI_Data.xlsx', sheet_name='Fig4I.1')

plot_dict = {'R':list(df_R.drop(df_R[df_R.Day > 3].index)['Prey']), 
             'R_vpsL':list(df_R_vpsL.drop(df_R_vpsL[df_R_vpsL.Day > 3].index)['Prey'])}

plt.figure(figsize=(3,4))
plt.boxplot(plot_dict.values(), showfliers=False)
plt.xticks(ticks=[1,2], labels=plot_dict.keys(), fontsize=20)

plt.scatter((simple_beeswarm2(plot_dict['R']))+1, plot_dict['R'], alpha=0.5, color='black', s=75, marker='p')
plt.scatter((simple_beeswarm2(plot_dict['R_vpsL']))+2, plot_dict['R_vpsL'], alpha=0.5, color='black', s=75, marker='p')

plt.xticks(rotation=45)
plt.ylim(-50000, 700000)
plt.ylabel(r'Biovolume', fontsize=20)
plt.xlabel(r'Treatment', fontsize = 20)
plt.title('WT', fontsize = 20)
plt.savefig('FigureS13A.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(mannwhitneyu(plot_dict['R'], plot_dict['R_vpsL']))
#%%
"""
SI Figure S13B
Related to Figure 5

Box and whisker plot of rugose predation on chitin
"""

df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig4H.1')
df_pred =  pd.read_excel('SI_Data.xlsx', sheet_name='Fig4H.2') 

plot_dict = {'2hPI P-':list(df.drop(df[df.Day > 3].index)['Prey']), 
             '2hPI P+':list(df_pred.drop(df_pred[df_pred.Day > 3].index)['Prey']),
             '48hPI P-':list(df.drop(df[df.Day != 5].index)['Prey']), 
              '48hPI P+':list(df_pred.drop(df_pred[df_pred.Day != 5].index)['Prey']),}

plt.figure(figsize=(3,4))
plt.boxplot(plot_dict.values(), showfliers=False)
plt.xticks(ticks=[1,2, 3, 4], labels=plot_dict.keys(), fontsize=20)

plt.scatter((simple_beeswarm2(plot_dict['2hPI P-']))+1, plot_dict['2hPI P-'], alpha=0.5, color='black', s=75, marker='p')
plt.scatter((simple_beeswarm2(plot_dict['2hPI P+']))+2, plot_dict['2hPI P+'], alpha=0.5, color='red', s=75, marker='p')
plt.scatter((simple_beeswarm2(plot_dict['48hPI P-']))+3, plot_dict['48hPI P-'], alpha=0.5, color='black', s=75, marker='p')
plt.scatter((simple_beeswarm2(plot_dict['48hPI P+']))+4, plot_dict['48hPI P+'], alpha=0.5, color='red', s=75, marker='p')
plt.xticks(rotation=45)
plt.hlines(y=1, xmin=0.5, xmax=4.5, linewidth=1, linestyles='--', color='r')
plt.ylabel(r'Biovolume', fontsize=20)
plt.xlabel(r'Treatment', fontsize = 20)
plt.title('Rugose', fontsize = 20)

plt.ylim(-50000, 700000)
plt.savefig('FigureS13B.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(mannwhitneyu(plot_dict['48hPI P-'], plot_dict['48hPI P+']))
print(mannwhitneyu(plot_dict['2hPI P-'], plot_dict['2hPI P+']))
#%%
"""
SI Figure S13C
Related to Figure 5

Box and whisker plot of rugose, vpsL predation on chitin
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig4I.1')
df_pred =  pd.read_excel('SI_Data.xlsx', sheet_name='Fig4I.2') 

plot_dict = {'2hPI P-':list(df.drop(df[df.Day != 3].index)['Prey']), 
             '2hPI P+':list(df_pred.drop(df_pred[df_pred.Day != 3].index)['Prey']),
             '48hPI P-':list(df.drop(df[df.Day != 5].index)['Prey']), 
              '48hPI P+':list(df_pred.drop(df_pred[df_pred.Day != 5].index)['Prey']),}

plt.figure(figsize=(3,4))
plt.boxplot(plot_dict.values(), showfliers=False)
plt.xticks(ticks=[1,2, 3, 4], labels=plot_dict.keys(), fontsize=20)

plt.scatter((simple_beeswarm2(plot_dict['2hPI P-']))+1, plot_dict['2hPI P-'], alpha=0.5, color='black', s=75, marker='p')
plt.scatter((simple_beeswarm2(plot_dict['2hPI P+']))+2, plot_dict['2hPI P+'], alpha=0.5, color='red', s=75, marker='p')
plt.scatter((simple_beeswarm2(plot_dict['48hPI P-']))+3, plot_dict['48hPI P-'], alpha=0.5, color='black', s=75, marker='p')
plt.scatter((simple_beeswarm2(plot_dict['48hPI P+']))+4, plot_dict['48hPI P+'], alpha=0.5, color='red', s=75, marker='p')
plt.xticks(rotation=45)

plt.ylim(-50000, 700000)
plt.hlines(y=1, xmin=0.5, xmax=4.5, linewidth=1, linestyles='--', color='r')
plt.ylabel(r'Biovolume', fontsize=20)
plt.xlabel(r'Treatment', fontsize = 20)

plt.title('Rugose del vpsl', fontsize = 20)
plt.savefig('FigureS13C.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(mannwhitneyu(plot_dict['48hPI P-'], plot_dict['48hPI P+']))
print(mannwhitneyu(plot_dict['2hPI P-'], plot_dict['2hPI P+']))

