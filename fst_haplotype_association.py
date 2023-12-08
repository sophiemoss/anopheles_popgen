# Identify haplotypes within each of the windows which were identified using fst

# Store the haplotypes as numpy arrays
# https://anopheles-genomic-surveillance.github.io/workshop-7/module-3-haplotype-clustering.html

# %%
# use the dendogram function from malariagen

import numpy as np
import scipy.cluster
import matplotlib.pyplot as plt
import string


def dendrogram(
    haplotypes, 
    linkage_method='single', 
    metric='hamming', 
    orient='right', 
    size=(7,5)
):
    """Takes a 2D numpy array of values and performs hierarchical clustering 
    and plots dendrogram."""
    
    # perform clustering
    linkage_matrix = scipy.cluster.hierarchy.linkage(
        haplotypes, 
        method=linkage_method, 
        metric=metric
    )
    if metric == 'hamming': 
        # convert hamming to number of snp differences
        linkage_matrix[:,2] *= haplotypes.shape[1] 
    
    # plot a dendrogram
    figsize = size if orient == 'right' else size[::-1]
    fig, ax = plt.subplots(figsize=size) 
    z = scipy.cluster.hierarchy.dendrogram(
        linkage_matrix,
        orientation=orient,
        leaf_rotation=0 if orient=='right' else 90,
        # label leaf nodes using letters
        labels=[' '.join(map(str, row)) + " | " + str(letter) 
                for letter, row in zip(string.ascii_uppercase, haplotypes)], 
        ax=ax
    )
    
    # tidy up the plot
    if orient == "right":
        ax.set_xlabel("Distance (no. SNPs)") 
        ax.set_ylabel("Haplotypes") 
        # change the limit so we can easily see identical haplotypes
        ax.set_xlim(-0.05, np.max(z['dcoord']) + 0.2)
    else:
        ax.set_xlabel("Haplotypes")
        ax.set_ylabel("Distance (no. SNPs)")
        # change the limit so we can easily see identical haplotypes
        ax.set_ylim(-0.05, np.max(z['dcoord']) + 0.2)

## create a numpy array of haplotypes, eg. 

haplotypes = np.array(
    [[0, 1, 0, 1, 0, 0, 0, 0, 0, 0],  # hap A
     [0, 1, 0, 0, 0, 0, 0, 0, 0, 0],  # hap B
     [0, 1, 0, 1, 0, 0, 0, 0, 0, 0],  # hap C
     [0, 0, 1, 1, 0, 0, 1, 1, 1, 1],  # hap D 
     [0, 1, 0, 1, 0, 0, 0, 0, 0, 0],  # hap E
     [0, 0, 1, 1, 0, 1, 1, 1, 1, 1],  # hap F
     [0, 1, 0, 1, 0, 0, 1, 1, 1, 1],  # hap G
     [0, 1, 0, 1, 0, 0, 0, 0, 0, 0]],  # hap H
) 

# %% perform hierarchical clustering

dendrogram(haplotypes, linkage_method='single', metric='hamming', orient='top')
# %%
