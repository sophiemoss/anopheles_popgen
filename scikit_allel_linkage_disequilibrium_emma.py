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
os.chdir('/mnt/storage11/sophie/emma')
os.getcwd()

# %% make zarr file

allel.vcf_to_zarr('phased_tmp_norm_snps_filtered.vcf.gz', 'phased_tmp_norm_snps_filtered.zarr', fields='*', overwrite=True)

# %%
callset = zarr.open('phased_tmp_norm_snps_filtered.zarr', mode='r')
callset.tree(expand=True)

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

# %% 2. Filter for chromosome 3
mask_3 = all_chromosomes == '3'
positions_3 = all_positions[mask_3]

# %% 3. Find the indices of desired positions in 2L

desired_positions =[161500150,
315931756,
315931943,
315932144,
315932210,
315939224,
315983762,
315983763,
315984130,
315998386,
315998453,
315998530,
316080722]

# %% 
indices_3 = np.where(np.isin(positions_3, desired_positions))[0]

# this number may be smaller than the number of desired positions, 
# as it will be the number of these positions that are found in the callset

# %% 4. Subset the genotype array by adding the start index of 2L variants
start_index_3 = np.where(mask_3)[0][0]
corrected_indices = indices_3 + start_index_3
subset_genotype_all = genotype_all[mask_3][indices_3]

# %% convert the genotype array into an array of alternate allele counts
gn = subset_genotype_all.to_n_alt(fill=-1)
gn

# %%
r = allel.rogers_huff_r(gn)
r  # doctest: +ELLIPSIS
r ** 2

# %%
from scipy.spatial.distance import squareform
squareform(r ** 2)

# %%
import allel.stats.ld
from allel.util import asarray_ndim
gn = asarray_ndim(gn,2,dtype='i1')

# %% check inputs
from allel.compat import memoryview_safe
gn = memoryview_safe(gn)

# %% compute correlation coefficients
from allel.opt.stats import gn_pairwise_corrcoef_int8, gn_pairwise2_corrcoef_int8, \
    gn_locate_unlinked_int8
r = gn_pairwise_corrcoef_int8(gn)

# %% convenience for singletons
if r.size == 1:
    r = r[0]
return r

# %% Plot a matrix of genotype linkage disequilibrium values between all pairs of variants.
m = squareform(r ** 2)

############# m is not correct at the mo, identical SNPs are 0 instead of 1 #####################
# %% Create a filtered list of positions that are actually present in your dataset
filtered_positions = [positions_3[i] for i in indices_3]

# 1. Import the metadata
snp_metadata = pd.read_csv('snp_positions.txt', sep='\t')

# 2. Create a dictionary mapping the POS to the SNP
pos_to_snp = dict(zip(snp_metadata['POS'], snp_metadata['SNP']))

# 3. Generate the desired labels based on the filtered_positions
labels = [f"{pos} ({pos_to_snp[pos]})" for pos in filtered_positions if pos in pos_to_snp]

# %% Plot the matrix
fig, ax = plt.subplots(figsize=(10, 10))  # Adjust the figsize as needed
cax = ax.imshow(m, interpolation='none', cmap='plasma')  # Specify the colourmap
cbar = fig.colorbar(cax)

# Label the axes
ax.set_xticks(range(len(labels)))
ax.set_yticks(range(len(labels)))
ax.set_xticklabels(labels, rotation=90)  # Rotate x-labels for better visibility
ax.set_yticklabels(labels)
ax.set_xlabel("Variant Positions and SNPs, Chromosome 2L")
ax.set_ylabel("Variant Positions and SNPs, Chromosome 2L")
ax.set_title("Pairwise LD")

plt.show()

# %% Create dataframe of SNPs at positions present in my data
# Use filtered_positions to create labels
labels = [f"{pos} ({pos_to_snp.get(pos, 'SNP not found')})" for pos in filtered_positions]

# Generate a DataFrame from the LD square matrix with the proper labels
ld_df = pd.DataFrame(squareform(r), index=labels, columns=labels)

# Display the DataFrame (in a Jupyter Notebook, this would render as a HTML table)
ld_df

# Optionally, save the DataFrame to a CSV file
#ld_df.to_csv('filtered_ld_values_table.csv')

# %% Extract the sample IDs of the heterozygous samples

# This assumes you have a DataFrame df_samples with a column 'sample' 
# corresponding to the samples in the order they are in the genotype array.
sample_ids = df_samples['sample'].values

# This list will hold the heterozygous samples for each SNP
heterozygous_samples_per_snp = []

# Loop over each SNP in the subset
for i, snp_genotypes in enumerate(subset_genotype_all):
    # Find the indices of heterozygous samples for this SNP
    heterozygous_indices = np.where((snp_genotypes[:, 0] != snp_genotypes[:, 1]) & (snp_genotypes[:, 0] >= 0))[0]
    
    # Retrieve the sample IDs corresponding to the heterozygous samples
    heterozygous_samples = sample_ids[heterozygous_indices]
    
    # Store the heterozygous sample IDs in the list
    heterozygous_samples_per_snp.append((filtered_positions[i], heterozygous_samples.tolist()))

# Convert the results to a DataFrame for easier viewing and analysis
df_heterozygous_samples = pd.DataFrame(heterozygous_samples_per_snp, columns=['Position', 'Heterozygous_Samples'])

# Display the DataFrame
print(df_heterozygous_samples)

# Optionally, save to a CSV
#df_heterozygous_samples.to_csv('heterozygous_samples_per_snp.csv')



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

# %% Try another set of SNPs
# Extract genomic data for SNPs you are checking
# Assuming you have subsetted for '2L' chromosome and have positions_2L and genotype_all arrays
index_2416980 = np.where(positions_2L == 2416980)[0][0]
index_2429745 = np.where(positions_2L == 2429745)[0][0]

genotype_2422652 = genotype_all[index_2422652]
genotype_2429745 = genotype_all[index_2429745]

# Display the genotypes for both SNPs side-by-side
for sample_idx in range(genotype_2416980.shape[0]):
    print(f"Sample {sample_idx}: 2416980-{genotype_2416980[sample_idx]}, 2416980-{genotype_2416980[sample_idx]}")

# This prints the genotypes for each sample for the two SNPs, you can check manually if they are in linkage disequilibrium

# or you can check with bcftools
# bcftools view -r 2L:2416980-2416980,2L:2430424-2430424 VGSC_F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz

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
