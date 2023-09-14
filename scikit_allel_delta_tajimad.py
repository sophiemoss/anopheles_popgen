
# %% ################################ DELTA TAJIMA'S D #####################################

## SELECT SAMPLE POPULATION TO WORK WITH

# some selection tests only support biallelic variants, not multiallelic. 
# This filtered VCF should already be biallelic SNPs only.
# Note that the zarr file should have been made from a phased vcf
# You make the zarr file with allel.vcf_to_zarr('phased_vcf_file_name.vcf.gz', 'output_name.zarr', fields='*', overwrite=True)

# %%
import os
os.chdir('/mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf_agamP4/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022_combinedvcf_filteringsteps')
os.getcwd()

# %%
## have now phased VCFs using beagle

import numpy as np
import allel
import zarr
import pandas as pd
import sys

## convert phased, filtered, VCF file to zarr file
# %%
# allel.vcf_to_zarr('2022gambiaevcfphased.vcf.gz', '2022gambiaevcfphased.zarr', fields='*', overwrite=True)

# %%
callset = zarr.open('2022gambiaevcfphased.zarr', mode='r')
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

# %% select resistant samples from metadata by index value and store as array

res_samples = df_samples[df_samples['phenotype'] == 'resistant'].index.values
res_samples

sus_samples = df_samples[df_samples['phenotype'] == 'susceptible'].index.values
sus_samples

# %% select genotypes for variants for those resistant samples and susceptible samples and store as a genotype dask array

gt_res_samples = gt.take(res_samples, axis=1)
gt_res_samples

gt_sus_samples = gt.take(sus_samples, axis=1)
gt_sus_samples

# allel.moving_delta_tajima_d(ac1, ac2, size, start=0, stop=None, step=None
# compute the difference in Tajima's D between two populations in moving windows

# create allele counts arrays ac1 and ac2
# genotype arrays were made earlier in this script: gt_sus_samples and gt_res_samples

ac_sus_samples = gt_sus_samples.count_alleles()
ac_sus_samples

ac_res_samples = gt_res_samples.count_alleles()
ac_res_samples

# %% compute delta tajima d
# Compute the difference in Tajima’s D between two populations in moving windows.
# calculate this in non-overlapping 100kb windows across the genome

tajimad = allel.moving_delta_tajima_d(ac_sus_samples, ac_res_samples, 100000)

# %% works with a small sample

small_ac_sus_samples = ac_sus_samples[:100, :]
small_ac_res_samples = ac_res_samples[:100, :]
tajimad = allel.moving_delta_tajima_d(small_ac_sus_samples, small_ac_res_samples, 10)

# %% adding logging
# for reference this took 127 seconds with window size of 100000 and array sizes of (17085002, 2)

import logging
import time

logging.basicConfig(level=logging.INFO)

def calculate_tajimas_d_with_logging(ac_sus_samples, ac_res_samples, window_size):
    logging.info("Starting Tajima's D calculation.")
    
    start_time = time.time()
    try:
        tajimad = allel.moving_delta_tajima_d(ac_sus_samples, ac_res_samples, window_size)
    except Exception as e:
        logging.error(f"An error occurred during calculation: {e}")
        return None

    elapsed_time = time.time() - start_time
    logging.info(f"Tajima's D calculation completed. Time taken: {elapsed_time} seconds.")
    
    return tajimad

# Then call your function
tajimad = calculate_tajimas_d_with_logging(ac_sus_samples, ac_res_samples, 100000)

# %%
# create x axis for plotting Tajima-D. The x axis ticks need to be the mid-points 
# of each window. I calcualted tajimad in 200,000bp windows, and the length of tajimad is 85
# so there must be 85 windows each a size of 200,000 base pairs.
# Create an array of mid-points for each window:

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

window_size = 100000
mid_points = np.arange(window_size // 2, len(tajimad) * window_size, window_size)

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
    ac_sus_by_chrom[chrom] = ac_sus_chrom
    ac_res_by_chrom[chrom] = ac_res_chrom
    print(ac_sus_chrom.shape)
    print(ac_res_chrom.shape)

    # Alternatively, store them in separate variables
    exec(f'ac_sus_samples_{chrom} = ac_sus_chrom')
    exec(f'ac_res_samples_{chrom} = ac_res_chrom')
    print(ac_sus_samples_3R.shape)
    print(ac_res_samples_3R.shape)

# %% Calculate delta Tajima D for each chromosome using the separate allele counts arrays

delta_tajimad_2L = calculate_tajimas_d_with_logging(ac_sus_samples_2L, ac_res_samples_2L, 100000)
delta_tajimad_2R = calculate_tajimas_d_with_logging(ac_sus_samples_2R, ac_res_samples_2R, 100000)
delta_tajimad_3L = calculate_tajimas_d_with_logging(ac_sus_samples_3L, ac_res_samples_3L, 100000)
delta_tajimad_3R = calculate_tajimas_d_with_logging(ac_sus_samples_3R, ac_res_samples_3R, 100000)
delta_tajimad_Mt = calculate_tajimas_d_with_logging(ac_sus_samples_Mt, ac_res_samples_Mt, 100000)
delta_tajimad_X = calculate_tajimas_d_with_logging(ac_sus_samples_X, ac_res_samples_X, 100000)

# %% Make plots for each chromosome

# Chromosome 2L

window_size = 100000
mid_points = np.arange(window_size // 2, len(delta_tajimad_2L) * window_size, window_size)

y = delta_tajimad_2L  # Delta Tajima's D calculated values
x = mid_points  # Mid-point of each window

# Plot
fig, ax = plt.subplots(figsize=(10, 4))
sns.despine(ax=ax, offset=5)
ax.plot(x, y, 'k-', lw=.5)
ax.set_ylabel("Delta Tajima's D")
ax.set_xlabel("Chromosome 2L")
ax.ticklabel_format(style='plain', axis='both', useOffset=False)
plt.show()

# %% Chromosome 2R

window_size = 100000
mid_points = np.arange(window_size // 2, len(delta_tajimad_2R) * window_size, window_size)

y = delta_tajimad_2R  # Delta Tajima's D calculated values
x = mid_points  # Mid-point of each window

# Plot
fig, ax = plt.subplots(figsize=(10, 4))
sns.despine(ax=ax, offset=5)
ax.plot(x, y, 'k-', lw=.5)
ax.set_ylabel("Delta Tajima's D")
ax.set_xlabel("Chromosome 2R")
ax.ticklabel_format(style='plain', axis='both', useOffset=False)
plt.show()

# %% Chromosome 3L

window_size = 100000
mid_points = np.arange(window_size // 2, len(delta_tajimad_3L) * window_size, window_size)

y = delta_tajimad_3L  # Delta Tajima's D calculated values
x = mid_points  # Mid-point of each window

# Plot
fig, ax = plt.subplots(figsize=(10, 4))
sns.despine(ax=ax, offset=5)
ax.plot(x, y, 'k-', lw=.5)
ax.set_ylabel("Delta Tajima's D")
ax.set_xlabel("Chromosome 3L")
ax.ticklabel_format(style='plain', axis='both', useOffset=False)
plt.show()

# %% Chromosome 3R

window_size = 100000
mid_points = np.arange(window_size // 2, len(delta_tajimad_3R) * window_size, window_size)

y = delta_tajimad_3R  # Delta Tajima's D calculated values
x = mid_points  # Mid-point of each window

# Plot
fig, ax = plt.subplots(figsize=(10, 4))
sns.despine(ax=ax, offset=5)
ax.plot(x, y, 'k-', lw=.5)
ax.set_ylabel("Delta Tajima's D")
ax.set_xlabel("Chromosome 3R")
ax.ticklabel_format(style='plain', axis='both', useOffset=False)
plt.show()

# %% Chromosome Mt

window_size = 100000
mid_points = np.arange(window_size // 2, len(delta_tajimad_Mt) * window_size, window_size)

y = delta_tajimad_Mt  # Delta Tajima's D calculated values
x = mid_points  # Mid-point of each window

# Plot
fig, ax = plt.subplots(figsize=(10, 4))
sns.despine(ax=ax, offset=5)
ax.plot(x, y, 'k-', lw=.5)
ax.set_ylabel("Delta Tajima's D")
ax.set_xlabel("Chromosome Mt")
ax.ticklabel_format(style='plain', axis='both', useOffset=False)
plt.show()

# windows are too large for Mt and could be smaller for others 

# %% Chromosome X

window_size = 100000
mid_points = np.arange(window_size // 2, len(delta_tajimad_X) * window_size, window_size)

y = delta_tajimad_X  # Delta Tajima's D calculated values
x = mid_points  # Mid-point of each window

# Plot
fig, ax = plt.subplots(figsize=(10, 4))
sns.despine(ax=ax, offset=5)
ax.plot(x, y, 'k-', lw=.5)
ax.set_ylabel("Delta Tajima's D")
ax.set_xlabel("Chromosome X")
ax.ticklabel_format(style='plain', axis='both', useOffset=False)
plt.show()