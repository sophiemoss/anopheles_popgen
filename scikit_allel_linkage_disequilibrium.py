# Documentation https://scikit-allel.readthedocs.io/en/stable/stats/ld.html

# This script is for analysing linkage disequilibrium between SNPs, to see if they are associated
# allel.rogers_huff_r(gn)
# Estimate the linkage disequilibrium parameter r for each pair of variants using the method of Rogers and Huff (2008)


# Note that the zarr file should have been made from a phased vcf
# You make the zarr file with allel.vcf_to_zarr('phased_vcf_file_name.vcf.gz', 'output_name.zarr', fields='*', overwrite=True)

# %% import packages
import os
import zarr
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
import pandas as pd
import allel; print('scikit-allel', allel.__version__)

# %% set working directory
os.chdir('/mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf_agamP4/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022_combinedvcf_filteringsteps')
os.getcwd()

# %%
callset = zarr.open('2022gambiaevcfphased.zarr', mode='r')
#callset.tree(expand=True)

# %%
## import metadata
df_samples=pd.read_csv('metadata_gambiae_2022.csv',sep=',',usecols=['sample','year','country','island','phenotype'])
df_samples.head()
df_samples.groupby(by=['phenotype']).count

# %%
# create genotype array
genotype_all = allel.GenotypeArray(callset['calldata/GT'])

# %% subset genotype array for area of interest
# 1. Get the positions and chromosomes from your callset
all_chromosomes = callset['variants/CHROM'][:]
all_positions = callset['variants/POS'][:]

# 2. Filter for chromosome 2L
mask_2L = all_chromosomes == '2L'
positions_2L = all_positions[mask_2L]

# 3. Find the indices of desired positions in 2L
desired_positions = [2390177,
2391228,
2399997,
2400071,
2402466,
2407967,
2416980,
2422651,
2422652,
2429556,
2429617,
2429745,
2429897,
2429915,
2430424,
2430817,
2430863,
2430880,
2430881,
2431061,
2431079,
25429236,
25429235] # extend the list with all your positions

indices_2L = np.where(np.isin(positions_2L, desired_positions))[0]
# this number may be smaller than the number of desired positions, 
# as it will be the number of these positions that are found in the callset

# 4. Subset the genotype array by adding the start index of 2L variants

start_index_2L = np.where(mask_2L)[0][0]
corrected_indices = indices_2L + start_index_2L
subset_genotype_all = genotype_all.take(corrected_indices, axis=0)

# %% convert the genotype array into an array of alternate allele counts
gn = subset_genotype_all.to_n_alt(fill=-1)
gn

# %%
r = allel.rogers_huff_r(gn)
r  # doctest: +ELLIPSIS

# %%
r ** 2  # doctest: +ELLIPSIS

from scipy.spatial.distance import squareform
squareformr = squareform(r ** 2)

# %% Plot a matrix of genotype linkage disequilibrium values between all pairs of variants.
m = squareformr

# Create a filtered list of positions that are actually present in your dataset
filtered_positions = [positions_2L[i] for i in indices_2L]

# 1. Import the metadata
snp_metadata = pd.read_csv('/mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf_agamP4/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022_combinedvcf_filteringsteps/allele_freqs/snp_positions.txt', sep='\t')

# 2. Create a dictionary mapping the POS to the SNP
pos_to_snp = dict(zip(snp_metadata['POS'], snp_metadata['SNP']))

# 3. Generate the desired labels based on the filtered_positions
labels = [f"{pos} ({pos_to_snp[pos]})" for pos in filtered_positions if pos in pos_to_snp]

# Plot the matrix
fig, ax = plt.subplots(figsize=(10, 10))  # Adjust the figsize as needed
cax = ax.imshow(m, interpolation='none', cmap='viridis')  # Specify the colourmap
cbar = fig.colorbar(cax)

# Label the axes
ax.set_xticks(range(len(labels)))
ax.set_yticks(range(len(labels)))
ax.set_xticklabels(labels, rotation=90)  # Rotate x-labels for better visibility
ax.set_yticklabels(labels)
ax.set_xlabel("Variant Positions and SNPs")
ax.set_ylabel("Variant Positions and SNPs")
ax.set_title("Pairwise LD")

plt.show()

# %% look at the different matplotlib colour maps, all can be inverted using _r, for example plasma_r
import numpy as np
import matplotlib.pyplot as plt

cmaps = ['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'cubehelix']
n = len(cmaps)

data = np.random.random((10,10))
fig, axs = plt.subplots(1, n, figsize=(15, 3), constrained_layout=True)
for ax, cmap in zip(axs, cmaps):
    psm = ax.pcolormesh(data, cmap=cmap, rasterized=True, vmin=0, vmax=1)
    fig.colorbar(psm, ax=ax)
    ax.set_title(cmap)
plt.show()


# %% Sense check results of linkage disequilibrium by extracting the genotype data for the specific SNPs of interest (in this case, 2416980 and 2430424) and then manually inspecting the relationship between them. I
# Extract genomic data for SNPs you are checking
# Assuming you have subsetted for '2L' chromosome and have positions_2L and genotype_all arrays
index_2416980 = np.where(positions_2L == 2416980)[0][0]
index_2430424 = np.where(positions_2L == 2430424)[0][0]

genotype_2416980 = genotype_all[index_2416980]
genotype_2430424 = genotype_all[index_2430424]

# Display the genotypes for both SNPs side-by-side
for sample_idx in range(genotype_2416980.shape[0]):
    print(f"Sample {sample_idx}: 2416980-{genotype_2416980[sample_idx]}, 2430424-{genotype_2430424[sample_idx]}")

# This prints the genotypes for each sample for the two SNPs, you can check manually if they are in linkage disequilibrium

# or you can check with bcftools
bcftools view -r 2L:2416980-2416980,2L:2430424-2430424 VGSC_F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz
