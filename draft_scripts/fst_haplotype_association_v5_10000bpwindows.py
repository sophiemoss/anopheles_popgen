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
import scipy
import statsmodels.api as sm
import statsmodels.formula.api as smf
print('scikit-allel', allel.__version__)

## Create genotype array using the same code that was used when doing the fst calculations

# %%
working_directory = '/mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf_agamP4/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022_filtering'
callset_file = '/mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf_agamP4/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022_filtering/2022_gambiae.zarr'
chromosome = "3L"
os.chdir(working_directory)

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

# %% filter out non-segregating and filter out multiallelic variants (these should already be removed during filtering of vcf but just to check)
acu = allel.AlleleCountsArray(acs[pop1][:] + acs[pop2][:])
flt = acu.is_segregating() & (acu.max_allele() == 1) # makes a boolean array indicating if each position is TRUE or FALSE for these criteria
print('Filtered out variants that are not biallelic and not segregating. Now retaining', np.count_nonzero(flt), 'SNPs')

# %% 
# create the new pos array, the allele count arrays for each population, and 
# the new genotype array, all containing only the variants that passed the flt criteria

pos = pos_all.compress(flt)
ac1 = allel.AlleleCountsArray(acs[pop1].compress(flt, axis=0)[:, :2])
ac2 = allel.AlleleCountsArray(acs[pop2].compress(flt, axis=0)[:, :2])
genotype = genotype_all.compress(flt, axis=0)
genotype
print("Created genotype array")

# %% Check length of pos and genotype are the same
print("Length of pos variable:",len(pos))
print("Length of genotype variable:",len(genotype))

if len(pos_all)==len(genotype_all):
    print("Length of positions and genotypes in the genotype array are the same, script continuing")
else:
    print("Something is wrong with the genotype_all array as the length of pos_all and genotype_all are different. Stopping script.")
    sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line

# %% Bring in windows of significance as a pandas dataframe 
significant_windows=pd.read_csv(f'/mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf_agamP4/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022_filtering/fst_patterson/significant_fst_values_hist_99_resistant_susceptible_{chromosome}.csv',sep=',')
print("Brought in the windows that were found significant during the fst analysis")

# Split the window ranges into start and end and convert them to integers in case they're floats
significant_windows['start'], significant_windows['end'] = zip(*significant_windows['Window'].str.split('-').tolist())
significant_windows['start'] = significant_windows['start'].astype(float).astype(int)
significant_windows['end'] = significant_windows['end'].astype(float).astype(int)

# Make a new window range which is 10,000bp wide, with the window of interest in the middle
significant_windows['wide_start'], significant_windows['wide_end'] = zip(*significant_windows['Window'].str.split('-').tolist())
significant_windows['wide_start'] = significant_windows['start'].astype(float).astype(int)-4500
significant_windows['wide_end'] = significant_windows['end'].astype(float).astype(int)+4500
print("Made the windows bigger to look at haplotypes including and around those windows")

# %% Subset the genotype array by POS in the pos array to make a separate 
# Dictionary to store subset genotype arrays for each significant window
subsetted_genotype_arrays_by_window = {}

# Iterate through each row in the significant windows DataFrame
for _, row in significant_windows.iterrows():
    # Get wide_start and wide_end for the current window
    start = row['wide_start']
    end = row['wide_end']
    # Identify indices of positions within this window
    # This gives us the indices where the conditions (pos >= start and pos <= end) are True
    window_indices = np.where((pos >= start) & (pos <= end))[0]
    # Check if window_indices is empty (meaning no positions fall within this window)
    if window_indices.size == 0:
        print(f"No positions found within the window {start}-{end}. Skipping this window.")
        continue
    # Subset the genotype array using these indices
    subsetted_genotype_array = genotype.subset(sel0=window_indices)

    # Save the subset genotype array in the dictionary
    # Use a key that represents the window (e.g., 'wide_start-wide_end')
    subsetted_genotype_arrays_by_window[f"{start}-{end}"] = subsetted_genotype_array
    length = len(subsetted_genotype_array)
    # Print confirmation
    print(f"Subset genotype array created for window {start}-{end}, length = {length}")

print("Subsetted the genotype array to make individual genotype arrays for each of the significant windows")

## Note that the subsetted genotype arrays for each window are different lengths
## This is because there are a different number of SNPs within each window. 
## The length of a subset genotype array corresponds to the number of SNP positions that fall within the specified genomic window.

# Groovy! 

# %% Now make corresponding position arrays for each of the subsets. 
# Do this using the same 'window indices' that found the SNPs in the genotype array in each window.

# Dictionary to store subset position arrays for each significant window
window_position_arrays = {}

# Iterate through each row in the significant windows DataFrame
for _, row in significant_windows.iterrows():
    # Get wide_start and wide_end for the current window
    start = row['wide_start']
    end = row['wide_end']
    # Identify indices of positions within this window
    window_indices = np.where((pos >= start) & (pos <= end))[0]
    # Subset the pos array using these indices
    pos_window = pos[window_indices]
    # Save the subset position array in the dictionary
    # Use a key that represents the window (e.g., 'wide_start-wide_end')
    window_position_arrays[f"{start}-{end}"] = pos_window
    length = len(pos_window)
    # Print confirmation
    print(f"Subset position array created for window {start}-{end}, length = {length}")

print("Subsetted the pos array to make individual position arrays corresponding to those genotype arrays")

# %% Check that the lengths of the subsetted genotype arrays and subsetted pos arrays are all the same

# Iterate over each significant window and check the lengths
for window in subsetted_genotype_arrays_by_window.keys():
    # Retrieve the genotype array for the window
    genotype_array = subsetted_genotype_arrays_by_window[window]
    # Retrieve the position array for the window
    pos_array = window_position_arrays[window]

    # Check if lengths are the same
    if len(genotype_array) == len(pos_array):
        print(f"Window {window} has matching lengths for genotype and position arrays. Length: {len(genotype_array)}")
    else:
        print(f"Length mismatch in window {window}. Genotype length: {len(genotype_array)}, Position length: {len(pos_array)}")

print("You have subset the genotype array and the pos array for each significant window. The length of the subsetted genotype arrays and their corresponding position arrays are the same. Hurrah! Crack on")

## Now I want to look at haplotypes, so I want to convert each of these genotype arrays to a haplotype array.
## Haplotype arrays represent the specific combination of alleles along a single chromosome.
## They simplify the representation by considering each allele separately, thus converting the 3D genotype array 
## into a 2D array. The dimensions are (number of variants, number of chromosome copies across all samples) 
## which essentially doubles the second dimension of the genotype array since each sample contributes two haplotypes.

# %% Create a dictionary to store haplotype arrays for each window
window_haplotype_arrays = {}

# Iterate over each significant window and convert genotype arrays to haplotype arrays
for window, genotype_array in subsetted_genotype_arrays_by_window.items():
    # Convert the genotype array to a haplotype array
    haplotype_array = genotype_array.to_haplotypes()

    # Store the haplotype array in the dictionary
    window_haplotype_arrays[window] = haplotype_array
    shape = haplotype_array.shape

    print(f"Converted genotype array to haplotype array for window {window}. Shape = {shape}")

print("Converted all of the individual subsetted genotype arrays to haplotype arrays")

# %% Check how many unique haplotypes there are in each haplotype array
# This should be fewer than double your sample number (for diploid organisms, because each sample has two haplotypes at each window)
unique_haplotype_counts = {}

for window, hap_array in window_haplotype_arrays.items():
    # Convert the array to a NumPy array and then transpose it
    numpy_hap_array = np.array(hap_array)
    transposed_haplotypes = numpy_hap_array.T
    # Convert each haplotype to a tuple to make them hashable
    tuple_haplotypes = [tuple(haplotype) for haplotype in transposed_haplotypes]
    # Use a set to find unique haplotypes
    unique_haplotypes = set(tuple_haplotypes)
    # Store the count of unique haplotypes for this window
    unique_haplotype_counts[window] = len(unique_haplotypes)
    print(f"Window {window} has {unique_haplotype_counts[window]} unique haplotypes.")

###  Make a dendrogram for each haplotype array and save the linkage matrix to use later

# %% First, define the dendrogram function:
# Global dictionary to store linkage matrices
window_linkage_matrices = {}

def dendrogram(haplotypes, window_range, linkage_method='single', metric='hamming', orient='right', size=(7,5), title_font_size=10):
    """Takes a 2D numpy array of values, performs hierarchical clustering, 
    and plots dendrogram, storing the linkage matrix in a global dictionary."""
    
    # perform clustering
    linkage_matrix = sch.linkage(haplotypes, method=linkage_method, metric=metric)
    if metric == 'hamming': 
        # convert hamming to number of snp differences
        linkage_matrix[:,2] *= haplotypes.shape[1]
    
    # Store linkage matrix in the global dictionary
    window_linkage_matrices[window_range] = linkage_matrix
    
    # generate labels for each haplotype
    labels = [str(i+1) for i in range(haplotypes.shape[0])]
    
    # plot a dendrogram
    fig, ax = plt.subplots(figsize=size) 
    sch.dendrogram(
        linkage_matrix,
        orientation=orient,
        leaf_rotation=0 if orient=='right' else 90,
        labels=labels, 
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

    # Remove x-axis labels if orientation is 'top'
    #if orient == 'top':
     #   ax.set_xticklabels([])

     # Add title to the plot
    plt.title(f"Dendrogram for Chromosome: {chromosome} Window: {window}", fontsize = title_font_size)

    # Save the plot as PNG file
    filename = f"dendrogram_{chromosome}_{window}.png"
    plt.savefig(filename, bbox_inches='tight')
    plt.show()
    return linkage_matrix  # Optionally, return the linkage matrix if needed

# %% Plot the dendrograms. Loop over each haplotype array, which are stored in the window_haplotype_arrays dictionary, and plot the dendrogram.
for window, hap_array in window_haplotype_arrays.items():
    print(f"Plotting dendrogram for {window}")
    # Transpose the haplotype array to get haplotypes as rows
    transposed_hap_array = np.transpose(hap_array)
    # Call the dendrogram function and pass the window range
    dendrogram(transposed_hap_array, window, metric='hamming', linkage_method='single', orient='top')

print("Finished making dendrograms")

# %% Identify clusters of haplotypes with distance 1% of the maximum distance
# Most of those 'clusters' only contain one snp, or a couple of snps, as they need to be close together in distance
# Filter through them to only pull out clusters that have more than 20 haplotypes
# Assuming window_linkage_matrices is a dictionary containing your linkage matrices

window_clusters_with_many_haplotypes = {}
for window, linkage_matrix in window_linkage_matrices.items():
    # Determine the maximum distance in the linkage matrix
    max_distance = np.max(linkage_matrix[:, 2])

    # Set the height cutoff to 1% of the maximum distance
    height_cutoff = 0.01 * max_distance

    # Form flat clusters from the hierarchical clustering
    clusters = sch.fcluster(linkage_matrix, height_cutoff, criterion='distance')

    # Dictionary to store haplotypes in each cluster
    clustered_haplotypes = {}

    # Iterate over each cluster label and find corresponding haplotypes
    for i, cluster_label in enumerate(clusters):
        if cluster_label not in clustered_haplotypes:
            clustered_haplotypes[cluster_label] = []
        clustered_haplotypes[cluster_label].append(i)

    # Print out the clusters with more than 20 haplotypes or a message if no such clusters are found
    clusters_with_more_than_five = {k: v for k, v in clustered_haplotypes.items() if len(v) >= 20}

    # Only add to the dictionary if there are clusters with more than five haplotypes
    if clusters_with_more_than_five:
        window_clusters_with_many_haplotypes[window] = clusters_with_more_than_five
        print(f"Clusters for window {window} with many haplotypes:")
        for cluster_label, haplotype_indices in clusters_with_more_than_five.items():
            print(f" - Cluster {cluster_label} contains haplotypes: {haplotype_indices}")

# Check the dictionary after the loop
print("Finished processing all windows.")
print("Saved clusters with many haplotypes as a dictionary:")
print(window_clusters_with_many_haplotypes)

#%% write a function to identify which haplotypes index in the haplotype array are from which mosquito

def hap_to_mosq(hap_idx):
    '''convert from haplotype index back to mosquitoes index to discover which mosquitoes have which haplotypes'''
    return int(hap_idx/2)

# %% write a function to take the selected haplotypes that we think are important and use the hap_to_mosq function to identify which mosquitoes these relate to
def selected_hap_array_to_mosq_array(selected_haps):
    mosq_array = [0]*42
    for hap_idx in selected_haps:
        mosq_idx = hap_to_mosq(hap_idx)
        mosq_array[mosq_idx]+=1
    return mosq_array
    
# %%  write a function to take the mosquito index and identify which phenotype these mosquitoes have
def map_mosq_index_to_phenotype(mosq_idx):
    phenotype = df_samples.iloc[mosq_idx]["phenotype"]
    return phenotype

# %% create a dictionary called map_mosquitoes and use the selected_hap_array_to_mosq_array to identify which mosquitoes have the selected haplotypes
map_mosquitoes = {}
for window,cluster_dict in window_clusters_with_many_haplotypes.items():
    map_mosquitoes[window]={}
    for cluster,selected_haps in cluster_dict.items():
        map_mosquitoes[window][cluster] = selected_hap_array_to_mosq_array(selected_haps)
print(map_mosquitoes)

## make a dictionary of phenotype_scores with window, cluster, phenotype. Use the map_mosq_index_to_phenotype function to identify the halotypes of the mosquitoes and then give a score
## mosquitoes score 1 if one of their haplotypes matches one of the 'selected haplotypes' that we care about
## mosquitoes score 2 if both of their haplotypes match 'selected haplotypes' that we care about

# %% 
phenotype_scores = {}
for window,cluster_dict in map_mosquitoes.items():
    phenotype_scores[window] = {}
    for cluster, haplotypes_mapped_per_mosquito in cluster_dict.items():
        phenotype_scores[window][cluster] = {}
        for mosq_idx, haplotypes_mapped_per_mosquito_score in enumerate(haplotypes_mapped_per_mosquito):
            mosq_phenotype = map_mosq_index_to_phenotype(mosq_idx)
            if mosq_phenotype in phenotype_scores[window][cluster]:
                phenotype_scores[window][cluster][mosq_phenotype]+= haplotypes_mapped_per_mosquito_score
            else:
                phenotype_scores[window][cluster][mosq_phenotype] = haplotypes_mapped_per_mosquito_score

print(phenotype_scores)

# %% Convert the phenotype_scores to a pandas dataframe

phenotype_score_data = phenotype_scores
rows = []
for key1, nested_dict in phenotype_score_data.items():
    for key2, values in nested_dict.items():
        row = values.copy()
        row['Significant Window'] = key1
        row['Haplotype Cluster'] = key2
        rows.append(row)

# Convert to DataFrame
phenotype_scores_dataframe = pd.DataFrame(rows)

# Rearranging the columns order
phenotype_scores_dataframe = phenotype_scores_dataframe[['Significant Window', 'Haplotype Cluster', 'control', 'resistant', 'susceptible']]

# %% Make separate tables for each haplotype cluster
df = phenotype_scores_dataframe

# Iterate through each unique combination
for window, cluster in df[['Significant Window', 'Haplotype Cluster']].drop_duplicates().values:
    # Filter data for the current combination
    cluster_data = df[(df['Significant Window'] == window) & (df['Haplotype Cluster'] == cluster)]

    # Create a new DataFrame for the phenotype and score
    phenotype_data_per_haplotype_cluster = pd.DataFrame({
        'phenotype': ['resistant', 'susceptible'],
        'score': [cluster_data['resistant'].iloc[0], cluster_data['susceptible'].iloc[0]]
    })
    # File name
    file_name = f'haplotype_cluster_{chromosome}_{window}_{cluster}_phenotype_table.csv'
    # Save to CSV
    #phenotype_data_per_haplotype_cluster.to_csv(file_name, index=False)

# %% Run a GLM (logistic regression) with logit link function for each haplotype cluster

# Loop over each file in the directory
#for filename in os.listdir(working_directory):
#    if filename.endswith('_phenotype_table.csv'):
#        # Load the CSV file into a DataFrame
#        filepath = os.path.join(working_directory, filename)
#        df = pd.read_csv(filepath)
#        # Convert phenotype to a binary variable
#        df['response'] = df['phenotype'].map({'resistant': 1, 'susceptible': 0})
#        # Define the model formula
#        formula = 'response ~ score'
#        # Fit the GLM model with a logit link function
#        model = smf.glm(formula=formula, data=df, family=sm.families.Binomial()).fit()
#        # Prepare the output filename
#        output_filename = f'{filename.replace(".csv", "")}_glm_summary.txt'
#        # Save the model summary to a text file
#        with open(output_filename, 'w') as file:
#            file.write(model.summary().as_text())
#        print(f"GLM summary for {filename} saved to {output_filename}")
#
# %%
