
######################## SELECTION STATISTICS #########################

# some selection tests only support biallelic variants, not multiallelic. 
# This filtered VCF should already be biallelic SNPs only.
# Note that the zarr file should have been made from a phased vcf
# You make the zarr file with allel.vcf_to_zarr('phased_vcf_file_name.vcf.gz', 'output_name.zarr', fields='*', overwrite=True)

# %%
import os
os.chdir('/mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf_agamP4/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022_filtering')
os.getcwd()

# %%
import numpy as np
import allel
import zarr
import pandas as pd
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import random

## convert phased, filtered, VCF file to zarr file
# %%
# allel.vcf_to_zarr('2022gambiaevcfphased.vcf.gz', '2022gambiaevcfphased.zarr', fields='*', overwrite=True)

# %%
callset = zarr.open('2022_gambiae.zarr', mode='r')
#callset.tree(expand=True)

# %%
## convert zarr file to genotype array
gt = allel.GenotypeDaskArray(callset['calldata/GT'])
print(gt.shape)

# %%
## import metadata
df_samples=pd.read_csv('metadata_gambiae_2022.csv',sep=',',usecols=['sample','year','country','island','phenotype'])
df_samples.head()
df_samples.groupby(by=['phenotype']).count

# %% partition samples by phenotype using metadata, split by index value and store as array
res_samples = df_samples[df_samples['phenotype'] == 'resistant'].index.values
sus_samples = df_samples[df_samples['phenotype'] == 'susceptible'].index.values
all_samples = df_samples.index.values

# %% select genotypes for variants and store as a genotype dask array
# this array contains 23 columns (1 for each mosquito), 16 million rows for each variant, and stores the genotype for each.
gt_res_samples = gt.take(res_samples, axis=1)
gt_sus_samples = gt.take(sus_samples, axis=1)
gt_all_samples = gt.take(all_samples, axis=1)

# %% compute allele counts for samples and create allele counts array
# this array stores all of the varaints as rows and the count of alleles as columns (0, 1 etc)
ac_res = gt_res_samples.count_alleles(max_allele=8).compute()
ac_sus = gt_sus_samples.count_alleles(max_allele=8).compute()
ac_all = gt_all_samples.count_alleles(max_allele=8).compute()

# %% filter for those that are biallelic and store as a boolean array
# is.segregating() just finds variants where more than one allele is observed. is_non_segregating() finds non-segregating variants (where at most one allele is observed)
# I don't actually need to use either the biallelic or the segregating filter becuase my vcf split multiallelic sites to biallelic during filtering
res_seg_variants = ac_res.is_biallelic_01()
sus_seg_variants = ac_sus.is_biallelic_01()
all_seg_variants = ac_all.is_biallelic_01()

# %% remove variants that are on Y_unplaced using a boolean mask
chrom = callset['variants/CHROM'][:]
exclude_chrom = 'Y_unplaced'
res_seg_variants = res_seg_variants & (chrom != exclude_chrom)
sus_seg_variants = sus_seg_variants & (chrom != exclude_chrom)
all_seg_variants = all_seg_variants & (chrom != exclude_chrom)

# %% make an allele counts array from this
ac_res_seg = ac_res.compress(res_seg_variants, axis=0)
ac_sus_seg = ac_sus.compress(sus_seg_variants, axis=0)
ac_all_seg = ac_all.compress(all_seg_variants, axis=0)

# %% also make a genotype dask array of the segregating variants
gt_res_seg = gt_res_samples.compress(res_seg_variants, axis = 0)
gt_sus_seg = gt_sus_samples.compress(sus_seg_variants, axis = 0)
gt_all_seg = gt_all_samples.compress(all_seg_variants, axis = 0)

# %% convert this genotype array to haplotype array (we can do this because the original data was phased)
# the haplotype array is similar to the genotype array, but there are two columns per mosquito, one for each haplotype
h_res_seg = gt_res_seg.to_haplotypes().compute()
h_sus_seg = gt_sus_seg.to_haplotypes().compute()
h_all_seg = gt_all_seg.to_haplotypes().compute()

# %% also store chromosome of each variant as we need this for shading the plots later
chrom = callset['variants/CHROM'][:]
chrom_res_seg = chrom.compress(res_seg_variants, axis=0)
chrom_sus_seg = chrom.compress(sus_seg_variants, axis=0)
chrom_all_seg = chrom.compress(all_seg_variants, axis=0)

# %% we need variant positions
pos = callset['variants/POS'][:]
pos_res_seg = pos.compress(res_seg_variants, axis=0)
pos_sus_seg = pos.compress(sus_seg_variants, axis = 0)
pos_all_seg = pos.compress(all_seg_variants, axis = 0)

# %% some variants in 1000 genomes project have multiple variants at the same genomic position, 
# which causes problems for some selection tests in scikit-allel. 
# Let's check if there any of these.
count_multiple_variants = np.count_nonzero(np.diff(pos_res_seg == 0))

# check resistant samples

if count_multiple_variants == 0:
    print("No cases where there are multiple variants at the same genomic position, script will continue")
else:
    print("There are multiple variants at the same genomic position. This causes problems with some selection tests using sci-kit allel.")
    sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line

# check susceptible samples 

count_multiple_variants = np.count_nonzero(np.diff(pos_sus_seg == 0))

if count_multiple_variants == 0:
    print("No cases where there are multiple variants at the same genomic position, script will continue")
else:
    print("There are multiple variants at the same genomic position. This causes problems with some selection tests using sci-kit allel.")
    sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line

# %% H12 was calculated using phased biallelic SNPs in 1000 bp windows along the genome
# SNP windows, using the garuds_h function in scikit-allel.

# %% Calculate h values for resistant and susceptible samples
real_res_h1, real_res_h12, real_res_123, real_res_h2_h1 = allel.moving_garud_h(h_res_seg, 1000)

real_sus_h1, real_sus_h12, real_sus_h123, real_sus_h2_h1 = allel.moving_garud_h(h_sus_seg, 1000)

# %% The above variables are stored as numpy.nd arrays

max_real_res_h12 = np.max(real_res_h12)
print("Maximum res H12 value:", max_real_res_h12)

max_real_sus_h12 = np.max(real_sus_h12)
print("Maximum sus H12 value:", max_real_sus_h12)


##### Plotting #####

# %% Check the number of windows for each array 
num_windows_res = real_res_h12.shape[0]
num_windows_sus = real_sus_h12.shape[0]

print("Number of windows for resistant samples:", num_windows_res)
print("Number of windows for susceptible samples:", num_windows_sus)

# %% The h12 values are calculated in windows of 1000 SNPs. Each SNP has a POS value which is in the POS array.
window_size = 1000  # Define the window size as 1000 SNPs
num_windows = len(pos_res_seg) // window_size  # Calculate the number of windows in total

# %% Calculate the median genomic position for each window
# for each iteration (i) a segment of pos_res_seg is processed
median_positions_res = [np.median(pos_res_seg[i * window_size: (i + 1) * window_size]) for i in range(num_windows)]

# Plotting
if len(median_positions_res) == len(real_res_h12):
    plt.figure(figsize=(10, 6))
    plt.scatter(median_positions_res, real_res_h12, alpha=0.6)
    plt.xlabel('Median Genomic Position of 1000 SNP Windows')
    plt.ylabel('H12 Value in Resistant Samples')
    plt.title('H12 Values Against Median Genomic Position of SNP Windows')
    plt.show()
else:
    print("Mismatch in the lengths of position and H12 arrays. Further debugging needed.")

# %% 
# Identify windows with H12 values over 0.3
h12_threshold_mask = np.array(real_res_h12) > 0.2
# Extract median genomic positions of these windows
high_res_h12_positions = np.array(median_positions_res)[h12_threshold_mask]
# Print or use these positions as needed
print("Median genomic positions of windows with H12 > 0.2:", high_res_h12_positions)

# %%
window_size = 1000  # Define the window size as 1000 SNPs
num_windows = len(pos_sus_seg) // window_size  # Calculate the number of windows

# Calculate the median genomic position for each window
median_positions_sus = [np.median(pos_sus_seg[i * window_size: (i + 1) * window_size]) for i in range(num_windows)]

# Plotting
if len(median_positions_sus) == len(real_sus_h12):
    plt.figure(figsize=(10, 6))
    plt.scatter(median_positions_sus, real_sus_h12, alpha=0.6)
    plt.xlabel('Median Genomic Position of 1000 SNP Windows')
    plt.ylabel('H12 Value in Susceptible Samples')
    plt.title('H12 Values Against Median Genomic Position of SNP Windows')
    plt.show()
else:
    print("Mismatch in the lengths of position and H12 arrays. Further debugging needed.")

# %% Identify windows with H12 values over 0.3
h12_threshold_mask = np.array(real_sus_h12) > 0.3
# Extract median genomic positions of these windows
high_sus_h12_positions = np.array(median_positions_sus)[h12_threshold_mask]
# Print or use these positions as needed
print("Median genomic positions of windows with H12 > 0.3:", high_sus_h12_positions)

###### Done! Amazing! Now do permutations of H12 with any random samples from the population ######

# %%
################################ Permutations ##################################

# to calculate for all samples it would be:
# h1, h12, h123, h2_h1 = allel.moving_garud_h(h_all_seg, 1000)
# Use all samples to get an idea of what selective sweeps there are in the whole population

# Number of permutations
num_permutations = 5
# Initialize a list to store the H12 values for each permutation
permuted_h12_values = []

for i in range(num_permutations):
    # Use previously made haplotype array for all samples 
    # Calculate h12 values for the selected samples
    _, h12_score_selected, _, _ = allel.moving_garud_h(h_all_seg, size=1000)
    #print(f"Permutation {i}, H12 values: {h12_score_selected}")
    # Pairing H12 values with their corresponding SNP window indices
    snp_windows = np.arange(len(h12_score_selected))  # Generating SNP window indices
    permuted_h12_pairs = list(zip(snp_windows, h12_score_selected))  # Pairing indices with H12 values
    # Append the pairs to the list
    permuted_h12_values.append(permuted_h12_pairs)
    # Notify
    print(f"Calculated H12 values for permutation {i}")

print("Permutations all calculated")

# %% Plotting the permuted values
window_size = 1000
plt.figure(figsize=(12, 6))

for i, permuted_pairs in enumerate(permuted_h12_values):
    # Extracting SNP window indices and H12 values
    snp_windows, h12_values = zip(*permuted_pairs)
    # Calculating the median genomic position for each window
    median_positions_permuted = [np.median(pos_all_seg[snp_window * window_size: (snp_window + 1) * window_size]) for snp_window in snp_windows if len(pos_all_seg[snp_window * window_size: (snp_window + 1) * window_size]) > 0]
    # Plotting
    plt.scatter(median_positions_permuted, h12_values, color='plum', alpha=0.6, label=f'Permutation {i+1}')

plt.xlabel('Median Genomic Position of 1000 SNP Windows')
plt.ylabel('Permuted H12 Value')
plt.title('Permuted H12 Values Against Median Genomic Position of SNP Windows')
plt.show()

# %% Plot the permuted windows and real_h12 on the same graph with the same x axis

plt.figure(figsize=(12, 6))

# Plot permuted values
for i, permuted_pairs in enumerate(permuted_h12_values):
    # Extracting SNP window indices and H12 values
    snp_windows, h12_values = zip(*permuted_pairs)
    # Calculating the median genomic position for each window
    median_positions_permuted = [np.median(pos_all_seg[snp_window * window_size: (snp_window + 1) * window_size]) for snp_window in snp_windows if len(pos_all_seg[snp_window * window_size: (snp_window + 1) * window_size]) > 0]
    # Plotting
    plt.scatter(median_positions_permuted, h12_values, color='plum', alpha=0.6, label=f'Permutation {i+1}' if i == 0 else "")

# Plot real_res_h12 values
median_positions_res = [np.median(pos_res_seg[i * window_size: (i + 1) * window_size]) for i in range(num_windows_res)]
plt.scatter(median_positions_res, real_res_h12, color='blue', alpha=0.6, label='Real Resistant Samples')

# Setting up the plot
plt.xlabel('Median Genomic Position of 1000 SNP Windows')
plt.ylabel('H12 Value')
plt.title('Comparison of Permuted H12 Values and Real Resistant Sample H12 Values')
plt.legend()
plt.show()

# %% Calculate the 99th percentile value for each of the median_positions_permuted
# Number of SNP windows is the same as the length of the first permutation's data
#num_windows = len(permuted_h12_values[0])
# Initialize a list to collect H12 values for each window
#h12_values_per_window = [[] for i in range(num_windows)]
# Collecting H12 values for each window from each permutation
#for permuted_h12 in permuted_h12_values:
#    for snp_window_number, h12_value in permuted_h12:
#        h12_values_per_window[snp_window_number].append(h12_value)
# Now calculate the 99th percentile for each window
#h12_99th_percentiles = [np.percentile(h12_values, 99) for h12_values in h12_values_per_window]

# %% Identify the genomic positions of the peaks

# %% Identify genomic position of H12 values over 0.3 for the resistant samples
h12_threshold_mask = np.array(real_res_h12) > 0.3
# Extract median genomic positions of these windows
high_res_h12_positions = np.array(median_positions_res)[h12_threshold_mask]
# Print or use these positions as needed
print("Median genomic positions of windows with H12 > 0.3:", high_res_h12_positions)

# %% Identify genomic position of H12 values over 0.3 for the permuted values (all samples)

# take average of permuted values for each window
# Number of SNP windows is the same as the length of the first permutation's data
num_windows = len(permuted_h12_values[0])
# Initialize a list to collect H12 values for each window
h12_values_per_window = [[] for i in range(num_windows)]
# Collecting H12 values for each window from each permutation
for permuted_h12 in permuted_h12_values:
    for snp_window_number, h12_value in permuted_h12:
        h12_values_per_window[snp_window_number].append(h12_value)
# Now calculate the average for each window
h12_permuted_averages = [np.average(permuted_h12_values) for permuted_h12_values in h12_values_per_window]
# Identify genomic positions for values above threshold
h12_permuted_threshold_mask = np.array(h12_permuted_averages) > 0.3
# Extract median genomic positions of these windows
high_permuted_h12_positions = np.array(median_positions_permuted)[h12_permuted_threshold_mask]
# Print or use these positions as needed
print("Median genomic positions of windows with permuted H12 > 0.3:", high_permuted_h12_positions)

# identify which 

# %% PBS
# PBS uses FST to identify genomic regions showing greater evolutionary change in one group 
# (here, the resistant samples) 
# relative to a closely related group (susceptible samples) and an outgroup. While originally designed to detect
# positive selection, it has also been used to detect phenotypic association (Grau-Bové et al., 2021).
# Note, For both H12 and PBS, phenotype permutations were performed as for FST to filter out false positives
# caused by the presence of extended swept haplotypes.
# calculate in 1000 bp windows and plot against the genome. What is classed as significant?
# can I use the control samples as ac3 here?

# %% create ac3 from control samples - use a separate population
# select samples
#con_samples = df_samples[df_samples['phenotype'] == 'control'].index.values
#con_samples

# select genotypes for samples
#gt_con_samples = gt.take(con_samples, axis=1)
#gt_con_samples

# create allele counts array
#ac_con_samples = gt_con_samples.count_alleles()
#ac_con_samples

# %% compute PBS

# allel.pbs(ac_sus_samples, ac_res_samples, ac_con_samples, 1000)

