# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 15:39:22 2025

@author: holtj
"""

import glob
import re
import tifffile
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from collections import defaultdict
import seaborn as sns

# Directory with TIFFs
# Collect all TIFF files
tif_files = glob.glob("tiff_files/*.tif")

# Group files by sample (everything before "_chX_")
groups = defaultdict(dict)
pattern = re.compile(r"(.*)_ch(\d+)_")

for f in tif_files:
    #remove images that aren't full stacks
    if 'snap' in f:
        continue
    else:
        fname = f.split("/")[-1]  # just the filename
        match = pattern.match(fname)
        if match:
            prefix, ch = match.groups()
            groups[prefix][ch] = f

data = []

# Loop through groups and compare channels
for prefix, ch_dict in groups.items():
    if "1" in ch_dict and "3" in ch_dict:  # compare ch1 vs ch2 (change if needed)
        f1 = ch_dict["1"]
        f2 = ch_dict["3"]
        print(f"\nProcessing {prefix}: ch1 vs ch2")

        stack1 = tifffile.imread(f1).astype(float)
        stack2 = tifffile.imread(f2).astype(float)

        if stack1.shape != stack2.shape:
            print(f"Skipping {prefix}, shapes differ: {stack1.shape} vs {stack2.shape}")
            continue

        vals1 = np.log10(stack1.ravel()+1)
        vals2 = np.log10(stack2.ravel()+1)
        
        plt.imshow(np.log10(stack1[5]+0.0001), cmap='Blues', vmin=-4, vmax=np.max(np.log10(stack1[5]+0.0001)))
        plt.title('Ch1, VPS')
        plt.colorbar()
        plt.show()
        plt.imshow(np.log10(stack2[5]+0.0001), cmap='Greens', vmin=-4, vmax=np.max(np.log10(stack2[5]+0.0001)))
        plt.title('Ch3, RbmA')
        plt.colorbar()
        plt.show()      
        sns.jointplot(x=np.log10(stack1[5].ravel()+1), y=np.log10(stack2[5].ravel()+0.0001))
        plt.show()
        
        plt.hist(vals1, bins=100, label='ch1, VPS', alpha=0.5)
        plt.hist(vals2, bins=100, label='ch3, RbmA', alpha=0.5)
        plt.ylabel('Counts')
        plt.xlabel('log10(Fluorescence Intensity (A.U.))')
        plt.legend()
        plt.title(f"Pixel correlation: {prefix} (ch1 vs ch3)")
        plt.grid(False)
        plt.show()       
        # Subsample for plotting
        n_points = 5000
        if len(vals1) > n_points:
            idx = np.random.choice(len(vals1), size=n_points, replace=False)
            vals1 = vals1[idx]
            vals2 = vals2[idx]

        # Scatter plot
        plt.figure(figsize=(6, 6))
        # sns.jointplot(x=vals1, y=vals2)
        plt.scatter(vals1, vals2, alpha=0.5)
        plt.xlabel("Channel 1 fluorescence, VPS")
        plt.ylabel("Channel 3 fluorescence, RbmA")
        plt.title(f"Pixel correlation: {prefix} (ch1 vs ch3)")
        plt.grid(False)
        plt.show()

        # Pearson correlation
        r,p = stats.pearsonr(vals1, vals2)
        data.append(r.statistic)
        print(f"Pearson r = {r:.3f}, p = {p:.2e}")
#%%
print(data)