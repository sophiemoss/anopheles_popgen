## scikitallel_workflow
# run this script using
# python scriptname.py /path/to/working/directory /path/to/callset.zarr chromosomename
# You make the zarr file with allel.vcf_to_zarr('phased_vcf_file_name.vcf.gz', 'output_name.zarr', fields='*', overwrite=True)
# allel.vcf_to_zarr('2019_melas_phased.vcf.gz', '2019_melas_phased.zarr', fields='*', overwrite=True)
######################## CALCULATING FST #########################
# %% adapted jupyter notebook from http://alimanfoo.github.io/2015/09/21/estimating-fst.html

import os
import argparse
import zarr
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
import pandas as pd
import allel
import sys
from datetime import datetime

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
# %% choose sample populations to work with
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
# %% filter out variants that aren't segregating in the union of the two selected populations. 
# also filter out multiallelic variants (these should already be removed during filtering of vcf but just to check)
acu = allel.AlleleCountsArray(acs[pop1][:] + acs[pop2][:])
flt = acu.is_segregating() & (acu.max_allele() == 1)
print('Filtered out variants that are not segretating, or are not biallelic. Now retaining', np.count_nonzero(flt), 'SNPs')
# %% create the new genotype array
pos = pos_all.compress(flt)
ac1 = allel.AlleleCountsArray(acs[pop1].compress(flt, axis=0)[:, :2])
ac2 = allel.AlleleCountsArray(acs[pop2].compress(flt, axis=0)[:, :2])
genotype = genotype_all.compress(flt, axis=0)
genotype
print("Created genotype array")
# check that pos and genotype are the same size
if len(pos)==len(genotype):
    print("Length of positions and genotypes in the genotype array are the same, script continuing")
else:
    print("Something is wrong with the genotype_all array as the lenght of pos_all and genotype_all are different. Stopping script.")
    sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line
# %% PLOT FST using windowed patterson fst
# allel.windowed_patterson_fst(pos, ac1, ac2, size=None, start=None, stop=None, step=None, windows=None, fill=nan)
# allel.windowed_weir_cockerham_fst(pos, g, subpops, size=None, start=None, stop=None, step=None, windows=None, fill=nan, max_allele=None)
print("Plotting Fst using allel.windowed_patterson_fst, window size 1000")
#subpoplist = [list(subpops['resistant']),
#                list(subpops['susceptible'])]
#fst, windows, counts = allel.windowed_weir_cockerham_fst(pos, genotype, subpoplist, size=1000)    # use the per-block average Fst as the Y coordinate
real_fst, real_windows, counts = allel.windowed_patterson_fst(pos, ac1, ac2, size=1000)    # use the per-block average Fst as the Y coordinate
real_y = real_fst
real_x = [np.mean(w) for w in real_windows]
# plot fst against chromosome position (bp)
fig, ax = plt.subplots(figsize=(10, 4))
sns.despine(ax=ax, offset=5)
ax.plot(real_x, real_y, 'k-', lw=.5)
ax.set_ylabel('$F_{ST}$')
ax.set_xlabel(f'Chromosome {chromosome} position (bp)')
ax.set_xlim(0, pos.max())
print("Number of Fst values:", len(real_x))
print("Number of windows:", len(real_y))
# save fst figure
timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
filename = f'fst_windowed_patterson_1000bp_{chromosome}_{pop1}_{pop2}.png'
plt.savefig(filename)
print("Saving windowed Fst plot")
# INSPECT Fst values
# %% Print the maximum FST value and its corresponding window

print("Inspecting windowed Fst plot for maximum value")
max_value = max(real_fst)
print(f"Max fst: {max_value}")
      
# %% Plot Fst values as a histogram to inspect Fst values
y = real_fst
# create histogram
fig, ax = plt.subplots(figsize=(10, 4))
ax.hist(real_y, bins=30, color='blue', edgecolor='black')  # You can adjust the number of bins as needed
ax.set_title('Histogram of Fst Values')
ax.set_xlabel('Fst Value')
ax.set_ylabel('Frequency')
# save this figure
timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
filename = f'fst_histogram_windowed_patterson_1000bp_{chromosome}_{pop1}_{pop2}.png'
plt.savefig(filename)
print("Saving histogram of Fst plot")

# %% Find the maximum and minimum fst values
# Zip together windows and fst values
zipped_windows_fst = list(zip(real_windows, real_fst))
# Find the maximum fst value and its corresponding window
max_fst_window = max(zipped_windows_fst, key=lambda real_x: real_x[1])
max_window, max_fst = max_fst_window
# Find the minimum fst value and its corresponding window
min_fst_window = min(zipped_windows_fst, key=lambda real_x: real_x[1])
min_window, min_fst = min_fst_window
# Printing or further processing
print(f"Maximum Fst value: {max_fst}, Window: {max_window}")
print(f"Minimum Fst value: {min_fst}, Window: {min_window}")
# set threshold for significant values as 3x the negative value, in the positive
hist_fst_threshold = (min_fst*3)*-1
print("Threshold for positive value being significant is:",hist_fst_threshold)

# %% save Fst values over this threshiold

hist_fst_threshold = [(window,value) for window, value in zip(real_windows, real_fst) if value >= hist_fst_threshold]
if hist_fst_threshold:
    with open(f'hist_fst_threshold_{pop1}_{pop2}_{chromosome}_window_size_1000.txt', 'w') as fst_file:
        fst_file.write("Window (Genomic bps)\FST Value\n")
        for window, value in hist_fst_threshold:
            fst_file.write(f"{window[0]}-{window[1]}\t{value}\n")
    print ("Saved FST values over histogram threshold")
else:
    print("No FST values over the histogram threshold were found")

# %% #################################################################
#################         PERMUTATIONS              ##################
######################################################################

# IMPORT METADATA
df_samples= pd.read_csv('metadata_gambiae_2022.csv',sep=',',usecols=['sample','year','country','island','phenotype'])
df_samples.head()
df_samples.groupby(by=['phenotype']).count()
print("Imported metadata")

# %%  do permutations of samples into pop1 and pop2
permuted_fst_values = []
for i in range(20):
    # Get the indices from df_samples
    df_res_sus_samples = df_samples[df_samples.phenotype.isin(['resistant', 'susceptible'])]
    indices = df_res_sus_samples.index.tolist()
    # Shuffle the indices
    np.random.shuffle(indices)
    # Split the indices into two groups
    half_length = len(indices) // 2
    pop1_indices = indices[:half_length]
    pop2_indices = indices[half_length:]
    subpops = {
        'permutation_pop1': pop1_indices,
        'permutation_pop2': pop2_indices,
    }
    permutation_pop1 = 'permutation_pop1'
    permutation_pop2 = 'permutation_pop2'
    # get allele counts
    acs = genotype_all.count_alleles_subpops(subpops)
    # filter variants out if not segregating or not biallelic     
    acu = allel.AlleleCountsArray(acs[permutation_pop1][:] + acs[permutation_pop2][:])
    flt = acu.is_segregating() & (acu.max_allele() == 1)
    # create the new genotype array
    pos = pos_all.compress(flt)
    ac1 = allel.AlleleCountsArray(acs[permutation_pop1].compress(flt, axis=0)[:, :2])
    ac2 = allel.AlleleCountsArray(acs[permutation_pop2].compress(flt, axis=0)[:, :2])
    genotype = genotype_all.compress(flt, axis=0)
    genotype
    # check that pos and genotype are the same size
    if len(pos)==len(genotype):
        print("Length of positions and genotypes in the genotype array are the same, script continuing")
    else:
        print("Something is wrong with the genotype_all array as the length of pos_all and genotype_all are different. Stopping script.")
        sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line
    
    # Calcalate and plot fst using windowed patterson fst
    fst, windows, counts = allel.windowed_patterson_fst(pos, ac1, ac2, size=1000)   # use the per-block average Fst as the Y coordinate
    print(f"Calculated fst using allel.windowed_patterson_fst, window size 1000 for permutation {i}")
    # Store Fst values for this iteration
    permuted_fst_values.append((windows,fst))

# %% Plot all fst values on the same graph from each of the permutations
fig, ax = plt.subplots(figsize=(10,4))
sns.despine(ax=ax, offset=5)
# Plot real fst values in red
ax.plot(real_x, real_fst, 'r-', lw=1.5, label='Real FST')  # 'r-' for red line
# Plot each set of permuted Fst values in grey
for windows,fst in permuted_fst_values:
    x = [np.mean(w) for w in windows]
    ax.plot(x, fst, 'k-', lw=.5, color='grey', alpha=0.5)  # grey color
ax.set_ylabel('Fst value')
ax.set_xlabel(f'Chromosome {chromosome} position (bp)')
ax.set_xlim(0, pos.max())
plt.legend()
plt.savefig('combined_fst_permutations_plot.png')
plt.show()

# %% Calculate the 99th percentile of the permuted values in each window
# Initialize an array to store the 99th percentile and window information
percentile_99th_values = []

# Number of windows
num_windows = len(permuted_fst_values[0][0])

# %% Iterate over each window
for i in range(num_windows):
    # Gather Fst values for this window from each permutation
    fst_values_for_window = [permuted[1][i] for permuted in permuted_fst_values]
    # Ignore NaN values in the percentile calculation
    fst_values_for_window = [value for value in fst_values_for_window if not np.isnan(value)]
    # Calculate the 99th percentile for this window
    percentile_99th = np.percentile(fst_values_for_window, 99) if fst_values_for_window else np.nan
    # Get the window boundaries
    window_boundaries = permuted_fst_values[0][0][i]
    # Store the window and the 99th percentile value as a tuple containing window_boundaries and 99th percentile value
    percentile_99th_values.append((window_boundaries, percentile_99th))

# percentile_99th_values contains each window, and then the value of the 99th percentile of that window

# %% Check if any of the real Fst values are above the 99th percentile value for each window
# to do this, I need to compare zipped_windows_fst with percentile_99th_values
# need to check if the value for Fst in zipped_windows_fst is greater than the percentile_99 in the percentile_99th_values list

# Initialize an empty list to store the data
significant_data = []

# Iterate through each window and its Fst value in zipped_windows_fst
for window_real_fst in zipped_windows_fst:
    window, real_fst_value = window_real_fst

    # Find the corresponding 99th percentile value for the window
    for window_percentile_values in percentile_99th_values:
        percentile_window, percentile_99th = window_percentile_values
        if np.array_equal(window, percentile_window):
            # Check if the Fst value is greater than the 99th percentile value
            if real_fst_value > percentile_99th:
                # Convert window array to a string representation for storage
                window_str = '-'.join(map(str, window))
                # Add the significant values to the list as a dictionary
                significant_data.append({'Window': window_str, 'Fst Value': real_fst_value, '99th Percentile': percentile_99th})

# Create DataFrame from the list
significant_values_df = pd.DataFrame(significant_data)

# save to a csv
significant_values_df.to_csv('significant_fst_values.csv', index=False)

# %% remove noise using minimum value * 3 using the threshold set earlier

# Assuming hist_fst_threshold is already defined and significant_values_df is your DataFrame

# bring significant values back from csv
# Read the CSV file into a DataFrame
significant_values_df = pd.read_csv('significant_fst_values_2L.csv')
# Display the first few rows of the DataFrame to verify
print(significant_values_df.head())
# Convert the 'Fst Value' column to numeric, in case it's not already
significant_values_df['Fst Value'] = pd.to_numeric(significant_values_df['Fst Value'])

# Filter the DataFrame to only include rows where 'Fst Value' is greater than hist_fst_threshold
filtered_df = significant_values_df[significant_values_df['Fst Value'] > hist_fst_threshold]

# Display the filtered DataFrame
print(filtered_df.head())

# Save the filtered DataFrame to a CSV file, if needed
filtered_df.to_csv('significant_fst_values_2L_filtered.csv', index=False)

# %% plot these on the combined fst graph somehow
