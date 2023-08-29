## scikitallel_workflow
# %%

conda activate scikit

######################## CALCULATING FST #########################
## I can do this from zarr file. 

## CONVERT FILTERED UNPHASED VCF TO ZARR FILE
# Try to convert to zarr and use workflow here: https://nbviewer.org/gist/alimanfoo/75b567b3d43810ef8eaef248b38b1c1c?flush_cache=true

import numpy as np
import allel
import zarr

allel.vcf_to_zarr('example.vcf', 'example.zarr', fields='*', overwrite=True)
import zarr
import numcodecs
print('zarr', zarr.__version__, 'numcodecs', numcodecs.__version__)

callset = zarr.open('F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.zarr', mode='r')
callset.tree(expand=True)

## CONVERT ZARR FILE TO GENOTYPE ARRAY

gt = allel.GenotypeDaskArray(callset['calldata/GT'])
print(gt.shape)
#(17085002, 42, 2)
gt
#Out[25]: <GenotypeDaskArray shape=(17085002, 42, 2) dtype=int8>

## IMPORT METADATA

df_samples=pd.read_csv('metadata_gambiae_2022.csv',sep=',',usecols=['sample','year','country','island','phenotype'])
df_samples.head()
df_samples.groupby(by=['phenotype']).count()

# obtain indices of resistant samples in the dataset
resistant_samples = df_samples[df_samples['phenotype'] == 'resistant'].index.values
resistant_samples

susceptible_samples = df_samples[df_samples['phenotype'] == 'susceptible'].index.values
susceptible_samples

control_samples = df_samples[df_samples['phenotype'] == 'control'].index.values
control_samples

# Define subpopulations
subpops = [susceptible_samples, resistant_samples, control_samples]

# Calculate Fst
# I want to analyse the component of variance between populations
# blen=None is not working, I need to specify a value for blen

fst = allel.weir_cockerham_fst(gt, subpops, max_allele=1, blen=1000)

a,b,c = allel.weir_cockerham_fst(gt, subpops, max_allele=1, blen=1000)

# a is the component of variance between populations
# b is the component of variance between individuals within populations
# c is the component of variance between gametes within individuals

# in simpler terms

# a: How much genetic diversity is there within a single group of individuals from the same population?
# b: How different are the genetic traits between different groups of individuals (populations)?
# c: Do certain genetic traits tend to appear together in similar patterns across different populations?

# estimate Fst for each variant and each allele individually:

fst = a / (a+b+c)

# estimate Fst for each variant individually, averaging over alleles:

fst = (np.sum(a, axis=1) /
       (np.sum(a, axis=1) + np.sum(b, axis=1) + np.sum(c, axis=1)))

# estimate Fst averaging over all variants and alleles:

fst = np.sum(a) / (np.sum(a) + np.sum(b) + np.sum(c))

# Create genomic position array

import dask.array as da
genomic_positions = da.arange(gt.shape[0])

# create a mapping between genomic positions and Fst values

import pandas as pd

# reshape the FST arrays to be 1-dimensional
a_1d = a.flatten()
b_1d = b.flatten()
c_1d = c.flatten()

mapping = pd.DataFrame({
    'genomic_position': genomic_positions.compute(),
    'FST_a': a_1d,
    'FST_b': b_1d,
    'FST_c': c_1d
})



## Seems a crazy low number of samples so using vcftools to calculate Fst too and comparing

vcftools --gzvcf F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz \
--weir-fst-pop resistant_samples.txt \
--weir-fst-pop susceptible_samples.txt \
--out fst_resistant_vs_susceptible

CHROM   POS     WEIR_AND_COCKERHAM_FST
2L      1206    -0.0343433
2L      1402    0.0194649
2L      1444    -0.0284405
2L      1447    -0.0416782
2L      1449    -0.0245136
2L      1462    -0.0170085
2L      1463    -0.0170085
2L      1470    -0.0365088
2L      1487    -0.0265077

## this calculates for individual variants. Be arare that calculating Fst is often more common to use sliding windows to aggregate vriants over larger genomic regions.

## after filtering only kept 33 out of 42 individuals
## there are plenty of higher Fst values using vcftools, so need to understand what is going wrong in the scikit allel workflow


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
allel.vcf_to_zarr('2022gambiaevcfphased.vcf.gz', '2022gambiaevcfphased.zarr', fields='*', overwrite=True)

# %%
callset = zarr.open('2022gambiaevcfphased.zarr', mode='r')
callset.tree(expand=True)

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

# %% Add red line to the plot showing significance

fig, ax = plt.subplots(figsize=(10, 3))
ax.plot(pos_res_seg, np.abs(ihs_res_std[0]), linestyle=' ', marker='o', mfc='none', mew=.5, mec='k')
ax.axhline(y=5, color='red', linestyle='--')
ax.set_xlabel('Genomic position (bp)')
ax.set_ylabel('$|IHS|$')
ax.set_ylim(0, 9)
ax.legend()

# %% list all positions with iHS value over certain threshold (5?)

threshold = 5
positions_above_threshold = pos_res_seg[ihs_res_std[0] >= threshold]
positions_above_threshold

# %%
############# Cross-population extended haplotype homozygosity (XPEHH) ###########

## Compute the unstandardized cross-population extended haplotype homozygosity score (XPEHH) for each variant.
## allel.xpehh(h1, h2, pos, map_pos=None, min_ehh=0.05, include_edges=False, gap_scale=20000, max_gap=200000, is_accessible=None, use_threads=True)














######################## PCA FOR GENES OF INTEREST #########################

## Make PCA for just genes of interest
## Use bcftools to subset the filtered vcf just for the genes of interest and make a PCA of that

# VGSC AgamP4_2L:2358158 to 2431617
# subset VCF

bcftools view -r 2L:2358158-2431617 F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz -Oz -o VGSC_only_F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz

######################## SPECIFIC SNP DETECTION #########################

# Do samples have the L995F mutation?

bcftools query -i 'POS==2422652' -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE:%GT]\n' \
-r 2L:2422652 F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz \
> 2422652_genotypes.txt

bcftools query -i 'POS==2429745' -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE:%GT]\n' \
-r 2L:2429745 F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz \
> 2429745_genotypes.txt


######################## ADMIXTURE #########################

## Admixture
## sci-kit allel
## Need an AlleleCountsArray, from original zarr file.

callset = zarr.open('F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.zarr', mode='r')
callset.tree(expand=True)

## CONVERT ZARR FILE TO GENOTYPE ARRAY
gt = allel.GenotypeDaskArray(callset['calldata/GT'])
print(gt.shape)

## Convert genotype array to Allele Counts Array
# you need to make an allele counts array for population A (aca) and an allele counts array for population B (acb)
ac = g.count_alleles()
ac
ac.dtype
ac.shape
ac.n_variants
ac.n_alleles

## Admixture - Emilia's methods https://github.com/LSHTMPathogenSeqLab/malaria-hub/tree/master/admixture
## can alter this for diploid organism.
## additional filtering for minor allele frequency and LD requried - look in published literature for An. gambiae.
## Paper - Genome variation and population structure among 1142 mosquitoes (2023) Ag1000G - used MAF > 0.01 and LD-pruning (how was LD pruning done?)
## Paper MAF >= 1%, 10 different seeds, no LD pruning (Miles, Nature 2017)

## will need to make a database with global anopheles gambiae for admixture to be interesting, want to compare it to An. gambiae from
## The Gambia and other countries in Africa. No phenotype data (?) but can create a separate database.


## 


# %%
