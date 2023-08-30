
######################## SELECTION STATISTICS #########################

## SELECT SAMPLE POPULATION TO WORK WITH

# some selection tests only support biallelic variants, not multiallelic. 
# This filtered VCF should already be biallelic SNPs only.

### Need to use phased version for conversion to haplotype array for further selection analyses.
# https://nbviewer.org/gist/alimanfoo/75b567b3d43810ef8eaef248b38b1c1c?flush_cache=true

# %%
import os
os.chdir('/mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf_agamP4/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022_combinedvcf_filteringsteps')
os.getcwd()

# %%
## have now phased VCFs using beagle

import numpy as np
np.__version__

import allel
allel.__version__

import zarr
zarr.__version__

import pandas as pd
pd.__version__

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

# %%
## working with resistant samples

res_samples = df_samples[df_samples['phenotype'] == 'resistant'].index.values
res_samples

# %%
## select genotypes for resistant samples

gt_res_samples = gt.take(res_samples, axis=1)
gt_res_samples

# %%
## select variants that are segregating within res_samples as only these will be informative

ac_res = gt_res_samples.count_alleles(max_allele=8).compute()
res_seg_variants = ac_res.is_segregating()
ac_res_seg = ac_res.compress(res_seg_variants, axis=0)
gt_res_seg = gt_res_samples.compress(res_seg_variants, axis = 0)
gt_res_seg

# %%
## this is from a phased VCF so we can convert this genotype array to haplotype array

h_res_seg = gt_res_seg.to_haplotypes().compute()
h_res_seg

# %%
## compute iHS values

# we need variant positions
pos = callset['variants/POS'][:]
pos_res_seg = pos.compress(res_seg_variants, axis=0)
pos_res_seg

# %%
# some variants in 1000 genomes have multiple variants at the same genomic position, which causes problems for some selection tests in scikit-allel. Let's check if there any of these.
np.count_nonzero(np.diff(pos_res_seg == 0))

# 0 so we are good to continue

# %%
# compute raw iHS

ihs_res_raw = allel.ihs(h_res_seg, pos_res_seg, use_threads=True, include_edges=True)
ihs_res_raw

# %%

%matplotlib inline
import matplotlib.pyplot as plt

# %%

fig, ax = plt.subplots()
ax.hist(ihs_res_raw[~np.isnan(ihs_res_raw)], bins=20)
ax.set_xlabel('Raw IHS')
ax.set_ylabel('Frequency (no. variants)');

# %% Standardize iHS

ihs_res_std = allel.standardize_by_allele_count(ihs_res_raw, ac_res_seg[:, 1])

# %% 
fig, ax = plt.subplots()
ax.hist(ihs_res_std[0][~np.isnan(ihs_res_std[0])], bins=20)
ax.set_xlabel('Standardized IHS')
ax.set_ylabel('Frequency (no. variants)');

# %% note that iHS has been calculated with unpolarized data, so only the magnitude of iHS
# is informative, not the sign.
ihs_res_std

# %% plot over the genome
fig, ax = plt.subplots(figsize=(10, 3))
ax.plot(pos_res_seg, np.abs(ihs_res_std[0]), linestyle=' ', marker='o', mfc='none', mew=.5, mec='k')
ax.set_xlabel('Genomic position (bp)')
ax.set_ylabel('$|IHS|$')
ax.set_ylim(0, 9);

# %% look for where the biggest signal is
idx_hit_max = np.nanargmax(ihs_res_std[0])
idx_hit_max

# %% genomic position of top hit
pos_res_seg[idx_hit_max]

# %% Visualise EHH decay around top hit with a Voight plot
# pull out haplotypes for region around top hit
flank_size = 2000
h_hit = h_res_seg[idx_hit_max - flank_size:idx_hit_max + flank_size + 1]
h_hit

# %%
fig = allel.fig_voight_painting(h_hit[:, h_res_seg[idx_hit_max] == 0], index=flank_size, height_factor=0.02)
fig.suptitle('Reference allele', y=1);

# %% 
fig = allel.fig_voight_painting(h_hit[:, h_res_seg[idx_hit_max] == 1], index=flank_size, height_factor=0.02)
fig.suptitle('Alternate allele', y=1);

# %%
# Display which iHS values are significant

# Where do we put the red line of significance? In Emilia's pop gen scripts, they
# are put where previous papers have had them. Default for iHS is 4. Default for XP-EHH is 5.
# Defauly for rsb is 5.
# I am using iHS significanc of 5 because otherwise it is too many SNPs at 4.

# %% Include red line in the plot showing significance

fig, ax = plt.subplots(figsize=(10, 3))
ax.plot(pos_res_seg, np.abs(ihs_res_std[0]), linestyle=' ', marker='o', mfc='none', mew=.5, mec='k')
ax.axhline(y=5, color='red', linestyle='--')
ax.set_xlabel('Genomic position (bp)')
ax.set_ylabel('$|IHS|$')
ax.set_ylim(0, 9)
ax.legend()

# %% list all positions with iHS value over certain threshold (5?)

res_ihs_positions_above_threshold_6 = pos_res_seg[ihs_res_std[0] >= 6]
res_ihs_positions_above_threshold_6

# Save positions_above_threshold to a text file
with open("res_ihs_positions_above_threshold_6.txt", "w") as file:
    for position in res_ihs_positions_above_threshold_6:
        file.write(str(position) + "\n")

res_ihs_positions_above_threshold_5 = pos_res_seg[ihs_res_std[0] >= 5]
res_ihs_positions_above_threshold_5

# Save positions_above_threshold to a text file
with open("res_ihs_positions_above_threshold_5.txt", "w") as file:
    for position in res_ihs_positions_above_threshold_5:
        file.write(str(position) + "\n")

# %% SUSCEPTIBLE POPULATION

### Now analyse iHS in susceptible population

sus_samples = df_samples[df_samples['phenotype'] == 'susceptible'].index.values
sus_samples

# %%
## select genotypes for resistant samples

gt_sus_samples = gt.take(sus_samples, axis=1)
gt_sus_samples

# %%
## select variants that are segregating within sus_samples as only these will be informative

ac_sus = gt_sus_samples.count_alleles(max_allele=8).compute()
sus_seg_variants = ac_sus.is_segregating()
ac_sus_seg = ac_sus.compress(sus_seg_variants, axis=0)
gt_sus_seg = gt_sus_samples.compress(sus_seg_variants, axis = 0)
gt_sus_seg

# %%
## this is from a phased VCF so we can convert this genotype array to haplotype array

h_sus_seg = gt_sus_seg.to_haplotypes().compute()
h_sus_seg

# %%
## compute iHS values

# we need variant positions
pos = callset['variants/POS'][:]
pos_sus_seg = pos.compress(sus_seg_variants, axis=0)
pos_sus_seg

# %%
# some variants in 1000 genomes have multiple variants at the same genomic position, which causes problems for some selection tests in scikit-allel. Let's check if there any of these.
np.count_nonzero(np.diff(pos_sus_seg == 0))

# 0 so we are good to continue

# %%
# compute raw iHS

ihs_sus_raw = allel.ihs(h_sus_seg, pos_sus_seg, use_threads=True, include_edges=True)
ihs_sus_raw

# %%

%matplotlib inline
import matplotlib.pyplot as plt

# %%

fig, ax = plt.subplots()
ax.hist(ihs_sus_raw[~np.isnan(ihs_sus_raw)], bins=20)
ax.set_xlabel('Raw IHS')
ax.set_ylabel('Frequency (no. variants)');

# %% Standardize iHS

ihs_sus_std = allel.standardize_by_allele_count(ihs_sus_raw, ac_sus_seg[:, 1])

# %% 
fig, ax = plt.subplots()
ax.hist(ihs_sus_std[0][~np.isnan(ihs_sus_std[0])], bins=20)
ax.set_xlabel('Standardized IHS')
ax.set_ylabel('Frequency (no. variants)');

# %% note that iHS has been calculated with unpolarized data, so only the magnitude of iHS
# is informative, not the sign.
ihs_sus_std

# %% plot over the genome
fig, ax = plt.subplots(figsize=(10, 3))
ax.plot(pos_sus_seg, np.abs(ihs_sus_std[0]), linestyle=' ', marker='o', mfc='none', mew=.5, mec='k')
ax.set_xlabel('Genomic position (bp)')
ax.set_ylabel('$|IHS|$')
ax.set_ylim(0, 9);

# %% look for where the biggest signal is
idx_hit_max = np.nanargmax(ihs_sus_std[0])
idx_hit_max

# %% genomic position of top hit
pos_sus_seg[idx_hit_max]

# %% Visualise EHH decay around top hit with a Voight plot
# pull out haplotypes for region around top hit
flank_size = 2000
h_hit = h_sus_seg[idx_hit_max - flank_size:idx_hit_max + flank_size + 1]
h_hit

# %%
fig = allel.fig_voight_painting(h_hit[:, h_sus_seg[idx_hit_max] == 0], index=flank_size, height_factor=0.02)
fig.suptitle('Reference allele', y=1);

# %% 
fig = allel.fig_voight_painting(h_hit[:, h_sus_seg[idx_hit_max] == 1], index=flank_size, height_factor=0.02)
fig.suptitle('Alternate allele', y=1);

# %%
# Display which iHS values are significant

# Where do we put the red line of significance? In Emilia's pop gen scripts, they
# are put where previous papers have had them. Default for iHS is 4. Default for XP-EHH is 5.
# Defauly for rsb is 5.
# I am using iHS significanc of 5 because otherwise it is too many SNPs at 4.

# %% Add red line to the plot showing significance

fig, ax = plt.subplots(figsize=(10, 3))
ax.plot(pos_sus_seg, np.abs(ihs_sus_std[0]), linestyle=' ', marker='o', mfc='none', mew=.5, mec='k')
ax.axhline(y=5, color='red', linestyle='--')
ax.set_xlabel('Genomic position (bp)')
ax.set_ylabel('$|IHS|$')
ax.set_ylim(0, 9)
ax.legend()

# %% list all positions with iHS value over certain threshold (5?)

sus_ihs_positions_above_threshold_6 = pos_sus_seg[ihs_sus_std[0] >= 6]

# Save positions_above_threshold to a text file
with open("sus_ihs_positions_above_threshold_6.txt", "w") as file:
    for position in sus_ihs_positions_above_threshold_6:
        file.write(str(position) + "\n")

sus_ihs_positions_above_threshold_5 = pos_sus_seg[ihs_sus_std[0] >= 5]

# Save positions_above_threshold to a text file
with open("sus_ihs_positions_above_threshold_5.txt", "w") as file:
    for position in sus_ihs_positions_above_threshold_5:
        file.write(str(position) + "\n")


# %%
############# Cross-population extended haplotype homozygosity (XPEHH) ###########

## Compute the unstandardized cross-population extended haplotype homozygosity score (XPEHH) for each variant.
## allel.xpehh(h1, h2, pos, map_pos=None, min_ehh=0.05, include_edges=False, gap_scale=20000, max_gap=200000, is_accessible=None, use_threads=True)
# create h1 and h2, selecting all variants instead of segregating variants only, which is what we did in iHS


# %%
## this is from a phased VCF so we can convert this genotype array to haplotype array

h_sus = gt_sus_samples.to_haplotypes().compute()
h_sus

h_res = gt_res_samples.to_haplotypes().compute()
h_res

ac_gt = gt.count_alleles(max_allele=8).compute()

# %%

# get variant positions

pos = callset['variants/POS']

# %% look at shapes of arrays
print("h_sus shape:", h_sus.shape)
print("h_res shape:", h_res.shape)

# %% compute xpehh
# xpehh_raw = allel.xpehh(h_sus, h_res, pos, map_pos=None, min_ehh=0.05, include_edges=False, gap_scale=20000, max_gap=20000, is_accessible=None, use_threads=True)

xpehh_raw = allel.xpehh(h_sus, h_res, pos, map_pos=None, include_edges=True, use_threads=True)

# %% look for where the biggest signal is
xpehh_hit_max = np.nanargmax(xpehh_raw)
xpehh_hit_max

# %% genomic position of top hit
pos[xpehh_hit_max]

# %%
%matplotlib inline
import matplotlib.pyplot as plt

# %%

fig, ax = plt.subplots()
ax.hist(xpehh_raw[~np.isnan(xpehh_raw)], bins=20)
ax.set_xlabel('Raw XP-EHH')
ax.set_ylabel('Frequency (no. variants)');

# %% Standardize XP-EHH

xpehh_std = allel.standardize_by_allele_count(xpehh_raw, ac_gt[:, 1])

# %% 
fig, ax = plt.subplots()
ax.hist(xpehh_std[0][~np.isnan(xpehh_std[0])], bins=20)
ax.set_xlabel('Standardized XP-EHH')
ax.set_ylabel('Frequency (no. variants)');

# %% note that iHS has been calculated with unpolarized data, so only the magnitude of iHS
# is informative, not the sign.
xpehh_std

# %% look at shapes

print("pos shape:", pos.shape)
print("xpehh_raw shape:", xpehh_raw.shape)

# %% plot on manhattan plot

fig, ax = plt.subplots(figsize=(10, 3))
ax.plot(pos[0], np.abs(xpehh_raw[0]), linestyle=' ', marker='o', mfc='none', mew=3, mec='k', label='$|XP-EHH|$')
ax.axhline(y=5, color='red', linestyle='--')
ax.set_xlabel('Genomic position (bp)')
ax.set_ylabel('$|XP-EHH|$')
ax.set_ylim(-1, 1)
ax.set_xlim(0,200000)
ax.legend()


# %%