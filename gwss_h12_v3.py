
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
## convert zarr file to genotype array (GenotypeDaskArray)
genotype_array = allel.GenotypeDaskArray(callset['calldata/GT'])
print(genotype_array.shape)

# %%
## import metadata
df_samples=pd.read_csv('metadata_gambiae_2022.csv',sep=',',usecols=['sample','year','country','island','phenotype'])
df_samples.head()
df_samples.groupby(by=['phenotype']).count

# %% Filter genotype array
# 1. Make allele count array so that we can filter for biallelic and segregating variants
ac_array = genotype_array.count_alleles(max_allele=8).compute() # ac_array is an AlleleCountsArray

# %% Filter for biallelic and segregating and store as boolean array
ac_array_filter = ac_array.is_biallelic_01()

# %% remove variants that are on Y_unplaced chromosome
chrom = callset['variants/CHROM'][:]
exclude_chrom = 'Y_unplaced'
ac_array_filter = ac_array_filter & (chrom != exclude_chrom)

# %% filter the allele counts array using the filter
filtered_ac_array = ac_array.compress(ac_array_filter, axis=0)

# %% also make a genotype dask array
filtered_gt = genotype_array.compress(ac_array_filter, axis = 0)

# %% partition samples by phenotype using metadata, split by index value and store as array
res_samples = df_samples[df_samples['phenotype'] == 'resistant'].index.values
sus_samples = df_samples[df_samples['phenotype'] == 'susceptible'].index.values
notres_samples = df_samples[df_samples['phenotype'] != 'resistant'].index.values
all_samples = df_samples.index.values

# %% select genotypes for variants and store as a genotype dask array
# this array contains 23 columns (1 for each mosquito), 16 million rows for each variant, and stores the genotype for each.
gt_res_samples = filtered_gt.take(res_samples, axis=1)
gt_sus_samples = filtered_gt.take(sus_samples, axis=1)
gt_all_samples = filtered_gt.take(all_samples, axis=1)
gt_notres_samples = filtered_gt.take(notres_samples, axis=1)

# %% convert this genotype array to haplotype array (we can do this because the original data was phased)
# the haplotype array is similar to the genotype array, but there are two columns per mosquito, one for each haplotype
h_res_seg = gt_res_samples.to_haplotypes().compute()
h_sus_seg = gt_sus_samples.to_haplotypes().compute()
h_all_seg = gt_all_samples.to_haplotypes().compute()
h_notres_seg = gt_notres_samples.to_haplotypes().compute()

# %% also store chromosome of each variant as we need this for shading the plots later
chrom = callset['variants/CHROM'][:]
chrom_filtered = chrom.compress(ac_array_filter, axis=0)

# %% we need variant positions
pos = callset['variants/POS'][:]
pos_filtered = pos.compress(ac_array_filter, axis=0)

# %% some variants in 1000 genomes project have multiple variants at the same genomic position, 
# which causes problems for some selection tests in scikit-allel. 
# Let's check if there any of these.
count_multiple_variants = np.count_nonzero(np.diff(pos_filtered== 0))

if count_multiple_variants == 0:
    print("No cases where there are multiple variants at the same genomic position, script will continue")
else:
    print("There are multiple variants at the same genomic position. This causes problems with some selection tests using sci-kit allel.")
    sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line

# %% H12 was calculated using phased biallelic SNPs in 1000 bp windows along the genome
# SNP windows, using the garuds_h function in scikit-allel.
# A.miles:
# To calibrate the window sizes I ran the H12 scans with a range of different window sizes, and chose
# the smallest window size for which the mean value of H1 over all windows was below 0.01.
# Lucas et al (2023) to identify regions in which swept haplotypes are more frequent in resistant compared to susceptible individuals, they calculated
# the difference in H12 value between groups, deltaH12.

# %% Calculate h values for resistant and susceptible samples
real_res_h1, real_res_h12, real_res_123, real_res_h2_h1 = allel.moving_garud_h(h_res_seg, 1000)

real_sus_h1, real_sus_h12, real_sus_h123, real_sus_h2_h1 = allel.moving_garud_h(h_sus_seg, 1000)

real_notres_h1, real_notres_h12, real_notres_123, real_notres_h2_h1 = allel.moving_garud_h(h_notres_seg, 1000)

# %% The above variables are stored as numpy.nd arrays

max_real_res_h12 = np.max(real_res_h12)
print("Maximum res H12 value:", max_real_res_h12)

max_real_sus_h12 = np.max(real_sus_h12)
print("Maximum sus H12 value:", max_real_sus_h12)

max_real_notres_h12 = np.max(real_notres_h12)
print("Maximum not resistant H12 value:", max_real_notres_h12)

########## Plotting ##########

# %% Check the number of windows for each array 
num_windows_res = real_res_h12.shape[0]
num_windows_sus = real_sus_h12.shape[0]
num_windows_notres = real_notres_h12.shape[0]

print("Number of windows for resistant samples:", num_windows_res)
print("Number of windows for not resistant:", num_windows_notres)
print("Number of windows for susceptible:", num_windows_sus)

# %% The h12 values are calculated in windows of 1000 SNPs. Each SNP has a POS value which is in the POS array.
window_size = 1000  # Define the window size as 1000 SNPs
num_windows = len(pos_filtered) // window_size  # Calculate the number of windows in total

# %% Calculate the median genomic position for each window
# for each iteration (i) a segment of pos_res_seg is processed
median_positions = [np.median(pos_filtered[i * window_size: (i + 1) * window_size]) for i in range(num_windows)]

# Plotting
if len(median_positions) == len(real_res_h12):
    plt.figure(figsize=(10, 6))
    plt.scatter(median_positions, real_res_h12, alpha=0.6)
    plt.xlabel('Median Genomic Position of 1000 SNP Windows')
    plt.ylabel('H12 Value in Resistant Samples')
    plt.title('H12 Values Against Median Genomic Position of SNP Windows')
    plt.show()
else:
    print("Mismatch in the lengths of position and H12 arrays. Further debugging needed.")

# %% 
# Identify windows with H12 values over 0.3
h12_threshold_mask = np.array(real_res_h12) > 0.3
# Extract median genomic positions of these windows
high_res_h12_positions = np.array(median_positions)[h12_threshold_mask]
# Print or use these positions as needed
print("Resistant samples: Median genomic positions of windows with H12 > 0.2:", high_res_h12_positions)

# %% Plot susceptible samples 

if len(median_positions) == len(real_sus_h12):
    plt.figure(figsize=(10, 6))
    plt.scatter(median_positions, real_sus_h12, alpha=0.6)
    plt.xlabel('Median Genomic Position of 1000 SNP Windows')
    plt.ylabel('H12 Value in Susceptible Samples')
    plt.title('H12 Values Against Median Genomic Position of SNP Windows')
    plt.show()
else:
    print("Mismatch in the lengths of position and H12 arrays. Further debugging needed.")

# %% Identify windows with H12 values over 0.3
h12_threshold_mask = np.array(real_sus_h12) > 0.3
# Extract median genomic positions of these windows
high_sus_h12_positions = np.array(median_positions)[h12_threshold_mask]
# Print or use these positions as needed
print("Susceptible samples - Median genomic positions of windows with H12 > 0.3:", high_sus_h12_positions)

# %% Plot both resistant and susceptible samples on the same graph
if len(median_positions) == len(real_res_h12) == len(real_sus_h12):
    plt.figure(figsize=(10, 6))
    plt.title('H12 Values Against Median Genomic Position of SNP Windows')

    # Plotting resistant samples
    plt.scatter(median_positions, real_res_h12, alpha=0.6, color='blue', label='Resistant')

    # Plotting susceptible samples
    plt.scatter(median_positions, real_sus_h12, alpha=0.6, color='plum', label='Susceptible')

    plt.xlabel('Median Genomic Position of 1000 SNP Windows')
    plt.ylabel('H12 Value')
    plt.legend()
    plt.show()
else:
    print("Mismatch in the lengths of position and H12 arrays. Further debugging needed.")

# %% Plot the not resistant samples and the resistant samples on the same graph

if len(median_positions) == len(real_notres_h12):
    plt.figure(figsize=(10, 6))
    plt.title('Genome Wide Selection Scan: H12 values across the genome')

    # Plotting resistant samples
    plt.scatter(median_positions, real_res_h12, alpha=0.6, color='blue', label='Resistant mosquitoes')

    # Plotting resistant samples
    plt.scatter(median_positions, real_notres_h12, alpha=0.6, color='plum', label='Non-resistant mosquitoes')

    plt.xlabel('Median Genomic Position of 1000 SNP Windows')
    plt.ylabel('H12 Value')
    plt.legend()
    plt.show()
else:
    print("Mismatch in the lengths of position and H12 arrays. Further debugging needed.")

# %% Identify genomic position of H12 values over 0.2 for the resistant samples
h12_threshold_mask = np.array(real_notres_h12) > 0.2
# Extract median genomic positions of these windows
h12_peaks = np.array(median_positions)[h12_threshold_mask]
# Print or use these positions as needed
print("Median genomic positions of windows with H12 > 0.2:", h12_peaks)

# %% What is heighest H12 value genomic position for each peak?
# peak b is the heighest value of all so can just find the highest value of notres and look for genomic position
max_real_notres_h12 = np.max(real_notres_h12)
max_value_peak_b = np.array(real_notres_h12) == max_real_notres_h12
peak_b_genomic_position = np.array(median_positions)[max_value_peak_b]

# peak a is more complicated
# need to find the highest res mosquito window position but only for windows lower than 2*1e7
# It will be one of these: 15130159.  15160300.  15200374. 15291513.5 15483510.5 15521589.5
# Iterate through the high H12 positions and print both the H12 value and the median genomic position of the SNP window

h12_threshold_mask = np.array(real_notres_h12) > 0.2

print("Median genomic positions and corresponding H12 values for windows with H12 > 0.2:")
for i in range(len(h12_threshold_mask)):
    if h12_threshold_mask[i]:
        print(f"Position: {median_positions[i]}, H12 Value: {real_notres_h12[i]}")

# so peak a is at the highest H12 value out of those positions, 15291513.5


# %%
while IFS=',' read -r chr start end; do
    echo "Processing: $chr from $start to $end" # Debugging line
    awk -v chr="$chr" -v start="$start" -v end="$end" \
        '$1 == chr && $4 <= end && $5 >= start' \
        /mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.56.chr.gff3 >> h12_peaks_genes.txt
done < h12_peaks_locations.txt

# %% Find the genomic position of the actual peak and then look at x distance around the peak for genes within the sweep.
# Or do what i've done and look at genes within the region of the peak (think this is fine)

# AGAP000818 = cytochrome P450, cyp9k1

