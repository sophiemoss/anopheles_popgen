
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

# %% variants could have multiple variants at the same genomic position, 
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

# %% Calculate h12 based on genomic windows of 1000bp

# Define genomic window size
window_size = 1000 

# %% Create the window boundaries
max_pos = max(pos_res_seg.max(), pos_sus_seg.max())
windows = np.arange(0, max_pos, window_size)

# %% Initialize arrays to hold H12 values for each window
h12_res = np.zeros(windows.size)
h12_sus = np.zeros(windows.size)

# %% Calculate H12 for each genomic window
for i, window_start in enumerate(windows[:-1]):
    window_end = windows[i+1]
    # Find indices of variants within the current window
    window_indices_res = (pos_res_seg >= window_start) & (pos_res_seg < window_end)
    window_indices_sus = (pos_sus_seg >= window_start) & (pos_sus_seg < window_end)
    
    # Extract haplotypes for the current window
    h_res_window = h_res_seg.compress(window_indices_res, axis=0)
    h_sus_window = h_sus_seg.compress(window_indices_sus, axis=0)
    
    # Calculate H12 for the window, if there are enough variants to calculate
    if h_res_window.shape[0] > 0:
        _, h12_res[i], _, _ = allel.moving_garud_h(h_res_window, size=h_res_window.shape[0])
    if h_sus_window.shape[0] > 0:
        _, h12_sus[i], _, _ = allel.moving_garud_h(h_sus_window, size=h_sus_window.shape[0])

# Now you have H12 values for resistant and susceptible samples for each genomic window

# %% You can then calculate delta H12
delta_h12 = h12_res - h12_sus

# %% Plot delta_h12 across the genome

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Let's assume you have these variables defined:
# delta_h12: an array of delta H12 values for each window
# window_size: the size of each window in base pairs

# Calculate the midpoints of each window
window_starts = np.arange(0, max_pos, window_size)
window_midpoints = window_starts + (window_size / 2)

# Create a DataFrame for easy plotting
plot_data = pd.DataFrame({
    'Window Midpoint': window_midpoints,
    'Delta H12': delta_h12
})

# Create the plot
plt.figure(figsize=(12, 4))
plt.plot(plot_data['Window Midpoint'], plot_data['Delta H12'], marker='o', linestyle='-', color='blue')

# Label your axes and title your plot
plt.xlabel('Genomic position (bp)')
plt.ylabel('Delta H12')
plt.title('Delta H12 Along the Genome')

# Set the limits of your x-axis if necessary
plt.xlim(0, max_pos)

# Show grid lines for better readability
plt.grid(True)

# Save the figure
plt.savefig('delta_h12.png', dpi=300)

# Show the plot
plt.show()

# I think this is not working properly to plot along the genome? It might be - plot in chromosomes and with scatter plot
# Check the logic for the windows
# %%
