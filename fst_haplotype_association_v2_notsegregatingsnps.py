# Identify haplotypes within each of the windows which were identified using fst

# Store the haplotypes as numpy arrays
# https://anopheles-genomic-surveillance.github.io/workshop-7/module-3-haplotype-clustering.html

# %%
# import packages

import numpy as np
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import string
import os
import argparse
import zarr
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
import pandas as pd
import allel
import sys
from datetime import datetime

## Create genotype array using the same code that was used when doing the fst calculations

# %%
working_directory = '/mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf_agamP4/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022_filtering'
callset_file = '/mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf_agamP4/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022_filtering/2022_gambiae.zarr'
chromosome = "2L"

# %% open callset file
callset = zarr.open(callset_file, mode='r')

# %% filter by chromosome. callset;'variants/CHROM'][:] accesses all entries in the
# CHROM column of the variant group within the Zarr datasaet. Then '== chromosome' compares
# each of the entries to the 'chromosome variable defined at the start of the script, and 
# creates a boolean array where each entry is True or False for the chromosome filter.
# np.where(chromosome_filter)[0] finds the indices where 'chromosome_filter' is True and then
# selects the corresponding POS column of the variants group in the Zarr dataset, and stores them in the pos_all array.
chromosome_filter = callset['variants/CHROM'][:] == chromosome
pos_all = callset['variants/POS'][np.where(chromosome_filter)[0]]
print("Chromosome being analysed:", chromosome)
print("Number of variants in chrom:", len(pos_all))

# %%  CONVERT ZARR FILE TO GENOTYPE ARRAY
# whole genome
# genotype_all = allel.GenotypeChunkedArray(callset['calldata/GT'])
# genotype_all
# create genotype array for just chrom
genotype_all = allel.GenotypeChunkedArray(callset['calldata/GT'][np.where(chromosome_filter)[0]])
# check length of pos_all and genotype_all are the same
print("Length of pos_all variable:",len(pos_all))
print("Length of genotype_all variable:",len(genotype_all))

if len(pos_all)==len(genotype_all):
    print("Length of positions and genotypes in the genotype array are the same, script continuing")
else:
    print("Something is wrong with the genotype_all array as the length of pos_all and genotype_all are different. Stopping script.")
    sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line

# %%  IMPORT METADATA
df_samples= pd.read_csv('metadata_gambiae_2022.csv',sep=',',usecols=['sample','year','country','island','phenotype'])
df_samples.head()
df_samples.groupby(by=['phenotype']).count()
print("Imported metadata")

# %% Choose sample populations to work with
pop1 = 'resistant'
pop2 = 'susceptible'
n_samples_pop1 = np.count_nonzero(df_samples.phenotype == pop1)
n_samples_pop2 = np.count_nonzero(df_samples.phenotype == pop2)
print("Population 1:", pop1, "Number of samples in pop1:", n_samples_pop1, "Population 2:", pop2, "Number of samples in pop2:", n_samples_pop2)

# %% dictonary mapping population names to sample indices
subpops = {
    pop1: df_samples[df_samples.phenotype == pop1].index,
    pop2: df_samples[df_samples.phenotype == pop2].index,
}
# %% get allele counts
acs = genotype_all.count_alleles_subpops(subpops)
acs

# %% filter out multiallelic variants (these should already be removed during filtering of vcf but just to check)
acu = allel.AlleleCountsArray(acs[pop1][:] + acs[pop2][:])
#flt = acu.is_segregating() & (acu.max_allele() == 1)
flt = (acu.max_allele() == 1)
print('Filtered out variants that are not biallelic. Now retaining', np.count_nonzero(flt), 'SNPs')

# %% create the new genotype array with the variants that passed the filters
pos = pos_all.compress(flt)
ac1 = allel.AlleleCountsArray(acs[pop1].compress(flt, axis=0)[:, :2])
ac2 = allel.AlleleCountsArray(acs[pop2].compress(flt, axis=0)[:, :2])
genotype = genotype_all.compress(flt, axis=0)
genotype
print("Created genotype array")

# %% Bring in windows of significance as a pandas dataframe 
significant_windows=pd.read_csv(f'significant_fst_values_hist_99_resistant_susceptible_{chromosome}.csv',sep=',')

# %% Subset the genotype array for each window of interest

# Split the window ranges into start and end and convert them to integers in case they're floats
significant_windows['start'], significant_windows['end'] = zip(*significant_windows['Window'].str.split('-').tolist())
significant_windows['start'] = significant_windows['start'].astype(float).astype(int)
significant_windows['end'] = significant_windows['end'].astype(float).astype(int)

# %% Initialize a dictionary to hold haplotypes for each window
window_haplotypes = {}

# %% Iterate over each window, subset genotype array to make one for just 
# that window, then convert it to a haplotype array, so you have a haplotype array for each window
for _, row in significant_windows.iterrows():
    start = row['start']
    end = row['end']

    # Subset the genotype data for this window
    window_mask = (pos >= start) & (pos <= end)
    genotype_window = genotype.compress(window_mask, axis=0)

    # Convert to haplotypes
    haplotypes_window = genotype_window.to_haplotypes()

    # Make a dictionary called window_haplotypes which stores each window and it's corresponding haplotype array
    window_haplotypes[f"{start}-{end}"] = haplotypes_window

# To select a particular array:
# Specify the window range of interest
# window_of_interest = "5811206-5812205"
# Access the haplotype array for this window
# specific_haplotype_array = window_haplotypes[window_of_interest]
# Now, 'specific_haplotype_array' contains the haplotype data for the specified window

# %% Check the number of haplotypes within the haplotype array for each window of interest:

# Iterate through each window in the window_haplotypes dictionary
for window_range, hap_array in window_haplotypes.items():
    # Get the number of haplotypes (rows) in the haplotype array
    num_haplotypes = hap_array.shape[0]
    # Print the window range and the number of haplotypes
    print(f"Window {window_range} has {num_haplotypes} haplotypes.")

# If there is a window with over 20 haplotypes, make a dendrogram and save the linkage matrix

# %% Create dendrogram function

def dendrogram(haplotypes, linkage_method='single', metric='hamming', orient='right', size=(7,5)):
    """
    Takes a 2D numpy array of values, performs hierarchical clustering, 
    and plots dendrogram. Returns the linkage matrix.
    """
    # perform clustering
    linkage_matrix = scipy.cluster.hierarchy.linkage(haplotypes, method=linkage_method, metric=metric)
    if metric == 'hamming': 
        # convert hamming to number of snp differences
        linkage_matrix[:,2] *= haplotypes.shape[1] 

    # plot a dendrogram
    figsize = size if orient == 'right' else size[::-1]
    fig, ax = plt.subplots(figsize=size) 
    scipy.cluster.hierarchy.dendrogram(
        linkage_matrix,
        orientation=orient,
        leaf_rotation=0 if orient=='right' else 90,
        labels=[' '.join(map(str, row)) for row in haplotypes], 
        ax=ax
    )

    # tidy up the plot
    if orient == "right":
        ax.set_xlabel("Distance (no. SNPs)") 
        ax.set_ylabel("Haplotypes") 
        ax.set_xlim(-0.05, np.max(linkage_matrix[:, 2]) + 0.2)
    else:
        ax.set_xlabel("Haplotypes")
        ax.set_ylabel("Distance (no. SNPs)")
        ax.set_ylim(-0.05, np.max(linkage_matrix[:, 2]) + 0.2)

    return linkage_matrix

# %% Dictionary to store linkage matrices for each window
linkage_matrices = {}

import scipy

# %% Iterate through each window in the window_haplotypes dictionary
for window_range, hap_array in window_haplotypes.items():
    # Check if the number of haplotypes is 20 or more
    if hap_array.shape[0] >= 20:
        print(f"Creating dendrogram for window {window_range} with {hap_array.shape[0]} haplotypes.")
        
        # If there are more than 20 haplotypes, create a dendrogram of haplotype clusters, 
        # with hierarchical clustering based on pairwise distances, and save the linkage matrix
        linkage_matrix = dendrogram(hap_array, linkage_method='single', metric='hamming', orient='top')
        
        # Store the linkage matrix with the window range as key
        linkage_matrices[window_range] = linkage_matrix

# %% If you want to access a specific linkage matrix
#specific_window = "5811206-5812205"
#specific_linkage_matrix = linkage_matrices.get(specific_window)

# Check if the linkage matrix exists for this window and then use it
#if specific_linkage_matrix is not None:
 #   print(f"Linkage matrix for window {specific_window}:")
  #  print(specific_linkage_matrix)
#else:
#    print(f"No linkage matrix found for window {specific_window}")


# %%  Go over each linkage matrix to determine if there are clusters of haplotypes 
# when the dendrogram is cut at 0.001, and check if any of these clusters contain 
# over 20 haplotypes

import scipy.cluster.hierarchy as sch
import numpy as np

height_cutoff = 0.001  # Height cutoff for clustering

# Iterate over each linkage matrix
for window_range, matrix in linkage_matrices.items():
    # Form flat clusters from the hierarchical clustering
    clusters = sch.fcluster(matrix, height_cutoff, criterion='distance')

    # Determine unique clusters and their sizes
    unique_clusters = np.unique(clusters)
    clustered_haplotypes = {cluster: [] for cluster in unique_clusters}
    
    for i, cluster_label in enumerate(clusters):
        clustered_haplotypes[cluster_label].append(i)

    # Check for clusters with more than 20 haplotypes and print information
    large_clusters = {cluster: haplotypes for cluster, haplotypes in clustered_haplotypes.items() if len(haplotypes) > 10}
    
    if large_clusters:
        print(f"Window {window_range} has the following large clusters:")
        for cluster, haplotype_indices in large_clusters.items():
            num_haplotypes = len(haplotype_indices)
            print(f"Cluster {cluster} contains {num_haplotypes} haplotypes: {haplotype_indices}")
    else:
        print(f"There are no clusters containing over 20 haplotypes in window {window_range}")

# %%
