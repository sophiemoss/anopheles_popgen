
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

# %% select resistant and susceptible samples from metadata by index value and store as array
res_samples = df_samples[df_samples['phenotype'] == 'resistant'].index.values
sus_samples = df_samples[df_samples['phenotype'] == 'susceptible'].index.values
all_samples = df_samples.index.values

# %% select genotypes for variants and store as a genotype dask array
gt_res_samples = gt.take(res_samples, axis=1)
gt_sus_samples = gt.take(sus_samples, axis=1)
gt_all_samples = gt.take(all_samples, axis=1)

# %% compute allele counts for samples and create allele counts array
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

if count_multiple_variants == 0:
    print("No cases where there are multiple variants at the same genomic position, script will continue")
else:
    print("There are multiple variants at the same genomic position. This causes problems with some selection tests using sci-kit allel.")
    sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line

count_multiple_variants = np.count_nonzero(np.diff(pos_sus_seg == 0))

if count_multiple_variants == 0:
    print("No cases where there are multiple variants at the same genomic position, script will continue")
else:
    print("There are multiple variants at the same genomic position. This causes problems with some selection tests using sci-kit allel.")
    sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line

# %% H12 was calculated using phased biallelic SNPs in 1000 bp windows along the genome
# SNP windows, using the garuds_h function in scikit-allel. 200 permutations in which penhotype labels randomly permuted and value recalculated?
# Calculate in 1000bp windows, look at the difference in H12
# A.miles:
# To calibrate the window sizes I ran the H12 scans with a range of different window sizes, and chose
# the smallest window size for which the mean value of H1 over all windows was below 0.01.
# Lucas et al (2023) to identify regions in which swept haplotypes are more frequent in resistant compared to susceptible individuals, they calculated
# the difference in H12 value between groups, deltaH12.

# %% Calculate h values for resistant and susceptible samples
real_res_h1, real_res_h12, real_res_123, real_res_h2_h1 = allel.moving_garud_h(h_res_seg, 1000)

real_sus_h1, real_sus_h12, real_sus_h123, real_sus_h2_h1 = allel.moving_garud_h(h_sus_seg, 1000)

# %% Calculate h values for all samples

# h1, h12, h123, h2_h1 = allel.moving_garud_h(h_all_seg, 1000)

# %% The above variables are stored as numpy.nd arrays

max_real_res_h12 = np.max(real_res_h12)
print("Maximum res H12 value:", max_real_res_h12)

max_real_sus_h12 = np.max(real_sus_h12)
print("Maximum sus H12 value:", max_real_sus_h12)


##### Plotting #####

# Check the number of windows for each array 
num_windows_res = real_res_h12.shape[0]
num_windows_sus = real_sus_h12.shape[0]

print("Number of windows for resistant samples:", num_windows_res)
print("Number of windows for susceptible samples:", num_windows_sus)

# %%
window_size = 1000  # Define the window size as 1000 SNPs
num_windows = len(pos_res_seg) // window_size  # Calculate the number of windows

# %% Calculate the median genomic position for each window
median_positions = [np.median(pos_res_seg[i * window_size: (i + 1) * window_size]) for i in range(num_windows)]

# Plotting
if len(median_positions) == len(real_res_h12):
    plt.figure(figsize=(10, 6))
    plt.scatter(median_positions, real_res_h12, alpha=0.6)
    plt.xlabel('Median Genomic Position of 1000 SNP Windows')
    plt.ylabel('H12 Value')
    plt.title('H12 Values Against Median Genomic Position of SNP Windows')
    plt.show()
else:
    print("Mismatch in the lengths of position and H12 arrays. Further debugging needed.")

# %% 
# Identify windows with H12 values over 0.2
h12_threshold_mask = np.array(real_res_h12) > 0.2
# Extract median genomic positions of these windows
high_h12_positions = np.array(median_positions)[h12_threshold_mask]
# Print or use these positions as needed
print("Median genomic positions of windows with H12 > 0.2:", high_h12_positions)

# %%
window_size = 1000  # Define the window size as 1000 SNPs
num_windows = len(pos_sus_seg) // window_size  # Calculate the number of windows

# Calculate the median genomic position for each window
median_positions = [np.median(pos_sus_seg[i * window_size: (i + 1) * window_size]) for i in range(num_windows)]

# Plotting
if len(median_positions) == len(real_sus_h12):
    plt.figure(figsize=(10, 6))
    plt.scatter(median_positions, real_sus_h12, alpha=0.6)
    plt.xlabel('Median Genomic Position of 1000 SNP Windows')
    plt.ylabel('H12 Value')
    plt.title('H12 Values Against Median Genomic Position of SNP Windows')
    plt.show()
else:
    print("Mismatch in the lengths of position and H12 arrays. Further debugging needed.")

# %% Identify windows with H12 values over 0.2
h12_threshold_mask = np.array(real_sus_h12) > 0.2
# Extract median genomic positions of these windows
high_h12_positions = np.array(median_positions)[h12_threshold_mask]
# Print or use these positions as needed
print("Median genomic positions of windows with H12 > 0.2:", high_h12_positions)

# %% Compute delta_h12, (H12 in resistant subset minus H12 in susceptible subset)

delta_h12 = real_res_h12 - real_sus_h12 
plt.plot(delta_h12)
plt.xlabel("Window index")
plt.ylabel("Delta H12 value")
plt.title("Delta H12 values across genomic windows")
plt.savefig('delta_h12.png', bbox_inches='tight')
plt.show()



# %% Now do 200 permutations for phenotype 

# Resistant samples first

# Initialize an array to store the permuted H12 values for each window
permuted_res_h12_values = []

for i in range(5):
    # Randomly permute the sample indices
    np.random.shuffle(all_samples)
    
    # Choose the indices for each group
    permuted_pop1_indices = all_samples[:23]
    permuted_pop2_indices = all_samples[23:33]

    # Calculate the allele counts for the permuted groups
    ac_permuted_pop1 = gt_all_seg.take(permuted_pop1_indices, axis=1).count_alleles().compute()
    ac_permuted_pop2 = gt_all_seg.take(permuted_pop2_indices, axis=1).count_alleles().compute()

    # Combine the counts and convert to haplotypes
    ac_combined = ac_permuted_pop1 + ac_permuted_pop2
    h_combined = allel.HaplotypeArray(gt_all_seg.compress(ac_combined.is_segregating(), axis=0).to_haplotypes().compute())

    # Calculate H12 for the combined haplotypes
    _, permuted_h12, _, _ = allel.moving_garud_h(h_combined, size=1000)

    # Store the permuted H12 values
    permuted_res_h12_values.append(permuted_h12)

# Convert the list of arrays into a 2D numpy array for easier percentile calculation later
permuted_res_h12_matrix = np.array(permuted_res_h12_values)

# Create dataframe of permuted resistance h12 values and save as csv
permuted_res_h12_values_df = pd.DataFrame(permuted_res_h12_values)
permuted_res_h12_values_df.to_csv(f'permuted_res_h12_values.csv', index=False)

# Plot the permuted_res_h12_values on the same plot as the real_res_h12 values.

# %% Plot all fst values on the same graph from each of the permutations
fig, ax = plt.subplots(figsize=(10,4))
sns.despine(ax=ax, offset=5)
# Plot each set of permuted Fst values in grey
for windows,fst in permuted_res_h12_values:
    x = [np.mean(w) for w in windows]
    ax.plot(x, fst, 'k-', lw=.5, color='grey', alpha=0.5, label='Permuted Fst')  # grey color
# Plot real fst values in blue
ax.plot(real_x, real_fst, color='#04B8D2', lw=1.5, label='Real Fst')  # 'r-' for red line


plt.legend()
plt.savefig(f'combined_fst_uneven_permutations_plot_{chromosome}.png')
plt.show()

if len(median_positions) == len(real_sus_h12):
    plt.figure(figsize=(10, 6))
    plt.scatter(median_positions, real_sus_h12, alpha=0.6)
    plt.xlabel('Median Genomic Position of 1000 SNP Windows')
    plt.ylabel('H12 Value')
    plt.title('H12 Values Against Median Genomic Position of SNP Windows')
    plt.show()
else:
    print("Mismatch in the lengths of position and H12 arrays. Further debugging needed.")
plt.legend()
plt.savefig(f'h12_res_samples_real_and_permuted.png')
plt.show()


# Calculate the 99th percentile value for each of the permuted fst windows









# %% Plot H12 values for Resistant Samples
plt.figure(figsize=(10, 6))
plt.plot(midpoints_res, resh12, color='blue')
plt.xlabel('Genomic Position (Midpoint of 1000bp Window)')
plt.ylabel('H12 Value')
plt.title('H12 Values Across Genomic Windows for Resistant Samples')
plt.show()

# Plot H12 values for Susceptible Samples
plt.figure(figsize=(10, 6))
plt.plot(midpoints_sus, sush12, color='red')
plt.xlabel('Genomic Position (Midpoint of 1000bp Window)')
plt.ylabel('H12 Value')
plt.title('H12 Values Across Genomic Windows for Susceptible Samples')
plt.show()

# %% Save samples H12 values
np.savetxt("resistant_h12_values.csv", resh12, delimiter=",")
np.savetxt("resistant_h12_values.csv", sush12, delimiter=",")

# %% Find positions with high h12 values resistant samples

max_resh12_index = np.argmax(resh12)  # Index of the maximum H12 value
max_resh12_value = resh12[max_resh12_index]  # Maximum H12 value

# Find the genomic position range for this window
window_size_bp = 1000  # Size of the window in base pairs
start_pos_of_window = pos_res_seg[max_resh12_index * window_size_bp]
end_pos_of_window = pos_res_seg[min((max_resh12_index + 1) * window_size_bp, len(pos_res_seg) - 1)]

print("Maximum H12 value:", max_resh12_value)
print("Genomic window for maximum H12:", start_pos_of_window, "-", end_pos_of_window)

# %% Find positions with high h12 values susceptible samples

max_sush12_index = np.argmax(sush12)  # Index of the maximum H12 value
max_sush12_value = sush12[max_sush12_index]  # Maximum H12 value

window_size_bp = 1000  # Size of the window in base pairs
start_pos_of_window = pos_sus_seg[max_sush12_index * window_size_bp]
end_pos_of_window = pos_sus_seg[min((max_sush12_index + 1) * window_size_bp, len(pos_res_seg) - 1)]

print("Maximum H12 value:", max_sush12_value)
print("Genomic window for maximum H12:", start_pos_of_window, "-", end_pos_of_window)

# %% Compute delta_h12

delta_h12 = resh12 - sush12

plt.plot(delta_h12)
plt.xlabel("Window index")
plt.ylabel("Delta H12 value")
plt.title("Delta H12 values across genomic windows")
plt.savefig('delta_h12.png', bbox_inches='tight')
plt.show()


# %% Plot delta_h12 values
plt.plot(delta_h12)
plt.xlabel("Window index")
plt.ylabel("Delta H12 value")
plt.title("Delta H12 values across genomic windows")
plt.savefig('delta_h12.png', bbox_inches='tight')
plt.show()

np.savetxt("delta_h12_values.csv", delta_h12, delimiter=",")











# %% PBS
# PBS uses FST to identify genomic regions showing greater evolutionary change in one group 
# (here, the resistant samples) 
# relative to a closely related group (susceptible samples) and an outgroup. While originally designed to detect
# positive selection, it has also been used to detect phenotypic association (Grau-Bové et al., 2021).
# Note, For both H12 and PBS, phenotype permutations were performed as for FST to filter out false positives
# caused by the presence of extended swept haplotypes.
# calculate in 1000 bp windows and plot against the genome. What is classed as significant?
# can I use the control samples as ac3 here?

# %% create ac3 from control samples
# select samples
con_samples = df_samples[df_samples['phenotype'] == 'control'].index.values
con_samples

# select genotypes for samples
gt_con_samples = gt.take(con_samples, axis=1)
gt_con_samples

# create allele counts array
ac_con_samples = gt_con_samples.count_alleles()
ac_con_samples

# %% compute PBS

allel.pbs(ac_sus_samples, ac_res_samples, ac_con_samples, 1000)


### To do:
# work out how to interpret H12 
# output needs to print both chromosome AND position for areas with high iHS and XP-EHH
# GWAS?

# Define variables
file="F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz"
position="2422652"
sample="NG-33833_res3_lib708245_10265_1"

# Extract the header line and find the column number for your sample
col_num=$(zgrep -m1 '^#CHROM' "$file" | tr '\t' '\n' | nl -v0 | grep -w "$sample" | cut -f1)

# Find the line with your position and extract the coverage information
zgrep -w "$position" "$file" | awk -v col="$col_num" '{print $col}' | cut -d':' -f2
