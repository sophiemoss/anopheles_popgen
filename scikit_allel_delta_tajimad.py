
# %% ################################ DELTA TAJIMA'S D #####################################

## SELECT SAMPLE POPULATION TO WORK WITH

# some selection tests only support biallelic variants, not multiallelic. 
# This filtered VCF should already be biallelic SNPs only.
# Note that the zarr file should have been made from a phased vcf
# You make the zarr file with allel.vcf_to_zarr('phased_vcf_file_name.vcf.gz', 'output_name.zarr', fields='*', overwrite=True)

# %%
import os
os.chdir('/mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf_agamP4/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022_filtering')
os.getcwd()

print("Starting script!")
# %%

import numpy as np
import allel
import zarr
import pandas as pd
import sys

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
res_samples

sus_samples = df_samples[df_samples['phenotype'] == 'susceptible'].index.values
sus_samples

# %% select genotypes for variants for those resistant samples and susceptible samples and store as a genotype dask array

gt_res_samples = gt.take(res_samples, axis=1)
gt_res_samples

gt_sus_samples = gt.take(sus_samples, axis=1)
gt_sus_samples

# %%
# create allele counts arrays ac1 and ac2
# genotype arrays were made earlier in this script: gt_sus_samples and gt_res_samples

ac_sus_samples = gt_sus_samples.count_alleles()
ac_sus_samples

ac_res_samples = gt_res_samples.count_alleles()
ac_res_samples

# %% compute delta tajima d
# Compute the difference in Tajima’s D between two populations in moving windows.
# calculate this in non-overlapping 1000 base pair windows across the genome
# for reference this took 195m 52.5s (3hrs 15)

#tajimad = allel.moving_delta_tajima_d(ac_sus_samples, ac_res_samples, 1000)

# Convert tajimad to a DataFrame
#tajimad_df = pd.DataFrame(tajimad)
# Save the DataFrame as a CSV file
#tajimad_df.to_csv('tajimad.csv', index=False)
# Later on, to load the tajimad array
tajimad_df = pd.read_csv('tajimad.csv')
# If you need it back as a numpy array:
tajimad = tajimad_df.to_numpy()

print ("Imported tajimad!")
# %% works with a small sample

#small_ac_sus_samples = ac_sus_samples[:100, :]
#small_ac_res_samples = ac_res_samples[:100, :]
#tajimad = allel.moving_delta_tajima_d(small_ac_sus_samples, small_ac_res_samples, 10)

# %% create logging function for calculating delta tajima d values

import logging
import time
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO)

def calculate_tajimas_d_with_logging(ac_sus_samples, ac_res_samples, window_size):
    logging.info("Starting Tajima's D calculation.")

    start_time = time.time()
    try:
        # Calculate Tajima's D values
        tajimad_values = allel.moving_delta_tajima_d(ac_sus_samples, ac_res_samples, window_size)
        
        # Calculate the start positions of each window
        window_starts = np.arange(0, ac_sus_samples.shape[0], window_size)
        
        # If the last window does not reach window_size, adjust the last start to the end of the array
        if window_starts[-1] + window_size > ac_sus_samples.shape[0]:
            window_starts = window_starts[:-1]  # Exclude the last start if it's out of bounds
        
        # Pair each Tajima's D value with its corresponding window start and end
        windows = [(start, min(start + window_size, ac_sus_samples.shape[0]), tajimad)
                   for start, tajimad in zip(window_starts, tajimad_values)]
        
    except Exception as e:
        logging.error(f"An error occurred during calculation: {e}")
        return None

    elapsed_time = time.time() - start_time
    logging.info(f"Tajima's D calculation completed. Time taken: {elapsed_time} seconds.")

    return windows


# Then call your function
#tajimad = calculate_tajimas_d_with_logging(ac_sus_samples, ac_res_samples, 1000)

# %%
# create x axis for plotting Tajima-D. The x axis ticks need to be the mid-points 
# of each window. I calculated tajimad in 1000bp windows, and the length of tajimad is 16425
# so there must be 16425 windows each a size of 1000 base pairs.
# Create an array of mid-points for each window:

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

window_size = 1000
mid_points = np.arange(window_size // 2, len(tajimad) * window_size, window_size)
length = len(mid_points)
print(f"The number of window midpoints is {length}")

# %% Plot Tajima D using these midpoints
# Prepare the data
y = tajimad  # Tajima's D calculated values
x = mid_points  # Mid-point of each window

# Plot
fig, ax = plt.subplots(figsize=(10, 4))
sns.despine(ax=ax, offset=5)
ax.plot(x, y, 'k-', lw=.5)
ax.set_ylabel("Delta Tajima's D")
ax.set_xlabel("Genomic position (bp)")
ax.ticklabel_format(style='plain', axis='both', useOffset=False)
plt.savefig('delta_tajimas_d_plot_wholegenome.png', bbox_inches='tight')  # Save the figure. 'bbox_inches='tight'' helps ensure that all labels are included in the saved figure.
plt.show()

# %% do this per chromosome - need to make genotype arrays for each chromosome and then transform these to allele count arrays

chrom_data = callset['variants/CHROM'][:]
print(chrom_data)

unique_chromosomes = set(chrom_data)
print('unique chromosomes:', unique_chromosomes)

# define the list of chromosomes you want to select
desired_chromosomes = ['X', '3R', '3L', '2R', '2L', 'Mt']

# select only those chromosomes that are present in unique_chromosomes
selected_chromosomes = [chrom for chrom in desired_chromosomes if chrom in unique_chromosomes]
print('selected chromosomes:', selected_chromosomes)

# %% Create dictionaries to store allele count arrays
ac_sus_by_chrom = {}
ac_res_by_chrom = {}

for chrom in selected_chromosomes:
    # Get the indices for this chromosome
    indices = np.where(chrom_data == chrom)[0]

    # Subset the genotype arrays based on the indices
    gt_sus_chrom = gt_sus_samples.take(indices, axis=0)
    gt_res_chrom = gt_res_samples.take(indices, axis=0)

    # Count alleles
    ac_sus_chrom = gt_sus_chrom.count_alleles()
    ac_res_chrom = gt_res_chrom.count_alleles()

    # Store them in the dictionaries
    #ac_sus_by_chrom[chrom] = ac_sus_chrom
    #ac_res_by_chrom[chrom] = ac_res_chrom
    #print(ac_sus_chrom.shape)
    #print(ac_res_chrom.shape)

    # Alternatively, store them in separate variables
    exec(f'ac_sus_samples_{chrom} = ac_sus_chrom')
    exec(f'ac_res_samples_{chrom} = ac_res_chrom')

 # %% Calculate delta Tajima D for each chromosome using the separate allele counts arrays
print("Calculating delta Tajima D for each chromosome!")

delta_tajimad_2L = calculate_tajimas_d_with_logging(ac_sus_samples_2L, ac_res_samples_2L, 1000)
delta_tajimad_2R = calculate_tajimas_d_with_logging(ac_sus_samples_2R, ac_res_samples_2R, 1000)
delta_tajimad_3L = calculate_tajimas_d_with_logging(ac_sus_samples_3L, ac_res_samples_3L, 1000)
delta_tajimad_3R = calculate_tajimas_d_with_logging(ac_sus_samples_3R, ac_res_samples_3R, 1000)
delta_tajimad_Mt = calculate_tajimas_d_with_logging(ac_sus_samples_Mt, ac_res_samples_Mt, 100)
delta_tajimad_X = calculate_tajimas_d_with_logging(ac_sus_samples_X, ac_res_samples_X, 1000)

print("Calculated delta_tajimad variables per chromosome!")

# %% Store delta_tajimad values for each chromosome comparison in a csv

# Convert to DataFrame
df_tajimad_2L = pd.DataFrame(delta_tajimad_2L, columns=['window_start', 'window_end', 'tajimas_d'])
df_tajimad_2R = pd.DataFrame(delta_tajimad_2R, columns=['window_start', 'window_end', 'tajimas_d'])
df_tajimad_3L = pd.DataFrame(delta_tajimad_3L, columns=['window_start', 'window_end', 'tajimas_d'])
df_tajimad_3R = pd.DataFrame(delta_tajimad_3R, columns=['window_start', 'window_end', 'tajimas_d'])
df_tajimad_Mt = pd.DataFrame(delta_tajimad_Mt, columns=['window_start', 'window_end', 'tajimas_d'])
df_tajimad_X = pd.DataFrame(delta_tajimad_X, columns=['window_start', 'window_end', 'tajimas_d'])

# Save to CSV
df_tajimad_2L.to_csv('delta_tajimad_2L.csv', index=False)
df_tajimad_2R.to_csv('delta_tajimad_2R.csv', index=False)
df_tajimad_3L.to_csv('delta_tajimad_3L.csv', index=False)
df_tajimad_3R.to_csv('delta_tajimad_3R.csv', index=False)
df_tajimad_Mt.to_csv('delta_tajimad_Mt.csv', index=False)
df_tajimad_X.to_csv('delta_tajimad_X.csv', index=False)

print("Stored delta_tajimad variables per chromosome!")

# %% Bring in delta_tajimad for each chromosome from previously saved csvs

# Read the CSV file

df_tajimad_2L = pd.read_csv('delta_tajimad_2L.csv')
df_tajimad_2R = pd.read_csv('delta_tajimad_2R.csv')
df_tajimad_3L = pd.read_csv('delta_tajimad_3L.csv')
df_tajimad_3R = pd.read_csv('delta_tajimad_3R.csv')
df_tajimad_Mt = pd.read_csv('delta_tajimad_Mt.csv')
df_tajimad_X = pd.read_csv('delta_tajimad_X.csv')

# Reconstruct the list of tuples
delta_tajimad_2L = list(df_tajimad_2L.itertuples(index=False, name=None))
delta_tajimad_2R = list(df_tajimad_2R.itertuples(index=False, name=None))
delta_tajimad_3L = list(df_tajimad_3L.itertuples(index=False, name=None))
delta_tajimad_3R = list(df_tajimad_3R.itertuples(index=False, name=None))
delta_tajimad_Mt = list(df_tajimad_Mt.itertuples(index=False, name=None))
delta_tajimad_X = list(df_tajimad_X.itertuples(index=False, name=None))

# Now delta_tajimad_2L is a list of tuples where each tuple is (window_start, window_end, tajimas_d) and can be used for plotting

# %% Make plots for each chromosome
print("Started making plots for each chromosome!")
# Chromosome 2L

window_size = 1000
mid_points = np.arange(window_size // 2, len(delta_tajimad_2L) * window_size, window_size)

y = [item[2] for item in delta_tajimad_2L]  # Delta Tajima's D calculated values
x = mid_points  # Mid-point of each window

# Plot
fig, ax = plt.subplots(figsize=(10, 4))
sns.despine(ax=ax, offset=5)
ax.plot(x, y, 'k-', lw=.5)
ax.set_ylabel("Delta Tajima's D")
ax.set_xlabel("Chromosome 2L")
ax.ticklabel_format(style='plain', axis='both', useOffset=False)
plt.savefig('delta_tajimas_d_2L.png', bbox_inches='tight')
plt.show()

# %% Chromosome 2R

window_size = 1000
mid_points = np.arange(window_size // 2, len(delta_tajimad_2R) * window_size, window_size)

y = [item[2] for item in delta_tajimad_2R]  # Delta Tajima's D calculated values
x = mid_points  # Mid-point of each window

# Plot
fig, ax = plt.subplots(figsize=(10, 4))
sns.despine(ax=ax, offset=5)
ax.plot(x, y, 'k-', lw=.5)
ax.set_ylabel("Delta Tajima's D")
ax.set_xlabel("Chromosome 2R")
ax.ticklabel_format(style='plain', axis='both', useOffset=False)
plt.savefig('delta_tajimas_d_2R.png', bbox_inches='tight')
plt.show()

# %% Chromosome 3L

window_size = 1000
mid_points = np.arange(window_size // 2, len(delta_tajimad_3L) * window_size, window_size)

y = [item[2] for item in delta_tajimad_3L]  # Delta Tajima's D calculated values
x = mid_points  # Mid-point of each window

# Plot
fig, ax = plt.subplots(figsize=(10, 4))
sns.despine(ax=ax, offset=5)
ax.plot(x, y, 'k-', lw=.5)
ax.set_ylabel("Delta Tajima's D")
ax.set_xlabel("Chromosome 3L")
ax.ticklabel_format(style='plain', axis='both', useOffset=False)
plt.savefig('delta_tajimas_d_3L.png', bbox_inches='tight')
plt.show()

# %% Chromosome 3R

window_size = 1000
mid_points = np.arange(window_size // 2, len(delta_tajimad_3R) * window_size, window_size)

y = [item[2] for item in delta_tajimad_3R]  # Delta Tajima's D calculated values
x = mid_points  # Mid-point of each window

# Plot
fig, ax = plt.subplots(figsize=(10, 4))
sns.despine(ax=ax, offset=5)
ax.plot(x, y, 'k-', lw=.5)
ax.set_ylabel("Delta Tajima's D")
ax.set_xlabel("Chromosome 3R")
ax.ticklabel_format(style='plain', axis='both', useOffset=False)
plt.savefig('delta_tajimas_d_3R.png', bbox_inches='tight')
plt.show()

# %% Chromosome Mt

window_size = 100
mid_points = np.arange(window_size // 2, len(delta_tajimad_Mt) * window_size, window_size)

y = [item[2] for item in delta_tajimad_Mt]  # Delta Tajima's D calculated values
x = mid_points  # Mid-point of each window

# Plot
fig, ax = plt.subplots(figsize=(10, 4))
sns.despine(ax=ax, offset=5)
ax.plot(x, y, 'k-', lw=.5)
ax.set_ylabel("Delta Tajima's D")
ax.set_xlabel("Chromosome Mt")
ax.ticklabel_format(style='plain', axis='both', useOffset=False)
plt.savefig('delta_tajimas_d_Mt.png', bbox_inches='tight')
plt.show()

# windows are too large for Mt and could be smaller for others 

# %% Chromosome X

window_size = 1000
mid_points = np.arange(window_size // 2, len(delta_tajimad_X) * window_size, window_size)

y = [item[2] for item in delta_tajimad_X]  # Delta Tajima's D calculated values
x = mid_points  # Mid-point of each window

# Plot
fig, ax = plt.subplots(figsize=(10, 4))
sns.despine(ax=ax, offset=5)
ax.plot(x, y, 'k-', lw=.5)
ax.set_ylabel("Delta Tajima's D")
ax.set_xlabel("Chromosome X")
ax.ticklabel_format(style='plain', axis='both', useOffset=False)
plt.savefig('delta_tajimas_d_X.png', bbox_inches='tight')
plt.show()

# %% Identify regions that are interesting in the delta tajima D values
# A negative DELTA tajima d value, means that TajimaD is higher in the sus pop (1st) than the res pop (2nd)
# This could indicate that the res pop is experiencing more recent positive selection or a selective sweep when compared to the sus pop, leading to reduced diversity around some alleles
# A positive DELTA tajima d value, means that TajimaD is higher in the res pop (2nd) than the sus pop (1st)
# This could indicate balancing selection in the res pop when compared to the sus pop.

# Identify delta tajima D that are significantly higher or lower than the average across the genome
# tajimad calculated first in this script holds the values across the whole genome, what is average? what is signiciantly higher or lower?
# tajimad has already been standardised, so we can use common thresholds for significance (scikit retuns standardized delta Tajima's D)
# so values less than -2 are significantly low and greater than 2 are significantly high
# use a loop to identify significant values from the delta tajima d for each chromosome

# Define your delta Tajima's D arrays
delta_tajimad_chromosomes = {
    '2L': delta_tajimad_2L,
    '2R': delta_tajimad_2R,
    '3L': delta_tajimad_3L,
    '3R': delta_tajimad_3R,
    'Mt': delta_tajimad_Mt,
    'X': delta_tajimad_X
}

# Set significance thresholds
lower_threshold = -2
upper_threshold = 2

# Dictionary to hold significant values for each chromosome
significant_values = {}

# Loop through each chromosome's data
for chrom, values in delta_tajimad_chromosomes.items():
    # Find significant values and include chromosome name and window information
    significant_high = [(chrom, item[0], item[1], item[2]) for item in values if item[2] > upper_threshold]
    significant_low = [(chrom, item[0], item[1], item[2]) for item in values if item[2] < lower_threshold]

    # Store in the dictionary
    significant_values[chrom] = {
        'significant_high': significant_high,
        'significant_low': significant_low
    }

# Print or process the significant values
for chrom, sig_values in significant_values.items():
    print(f"Chromosome {chrom}:")
    print("Significant High:", sig_values['significant_high'])
    print("Significant Low:", sig_values['significant_low'])
    print("\n")

# %% Store significant Delta Tajima's D values in CSV file for future use
# Flatten the data into a list of dictionaries for easier conversion to a DataFrame
data = []
for chrom, sig_values in significant_values.items():
    for value in sig_values['significant_high']:
        data.append({'Chromosome': value[0], 'Window Start': value[1], 'Window End': value[2], 'Delta Tajima\'s D': value[3], 'Significance': 'High'})

    for value in sig_values['significant_low']:
        data.append({'Chromosome': value[0], 'Window Start': value[1], 'Window End': value[2], 'Delta Tajima\'s D': value[3], 'Significance': 'Low'})

# Create DataFrame and save to CSV
df = pd.DataFrame(data)
df.to_csv('significant_tajimas_d_threshold2.csv', index=False)

# %%
