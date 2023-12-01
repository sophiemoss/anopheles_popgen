
######################## SELECTION STATISTICS #########################

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

# %% select genotypes for variants for those resistant samples and store as a genotype dask array

gt_res_samples = gt.take(res_samples, axis=1)
gt_res_samples

# %% select variants that are segregating within res_samples as only these will be informative
## also some selection tests don't support multiallelic variants, so just keep biallelics
## for this pipeline the VCF is already filtered so should be no biallelic SNPs anyway

# %% compute allele counts for resistant samples and create allele counts array
ac_res = gt_res_samples.count_alleles(max_allele=8).compute()
# %% filter for those that are segregating and biallelic and store as a boolean array
# is.segregating() just finds variants where more than one allele is observed. is_non_segregating() finds non-segregating variants (where at most one allele is observed)

res_seg_variants = ac_res.is_segregating() & ac_res.is_biallelic_01()
# %% remove variants that are on Y_unplaced using a boolean mask
chrom = callset['variants/CHROM'][:]
exclude_chrom = 'Y_unplaced'
res_seg_variants = res_seg_variants & (chrom != exclude_chrom)
# %% make an allele counts array from this
ac_res_seg = ac_res.compress(res_seg_variants, axis=0)
# %% also make a genotype dask array
gt_res_seg = gt_res_samples.compress(res_seg_variants, axis = 0)
gt_res_seg

# %%
## this is from a phased VCF so we can convert this genotype array to haplotype array

h_res_seg = gt_res_seg.to_haplotypes().compute()
h_res_seg

# %% we need variant positions
pos = callset['variants/POS'][:]
pos_res_seg = pos.compress(res_seg_variants, axis=0)
pos_res_seg

# %% also store chromosome of each variant as we need this for shading the plots later
chrom = callset['variants/CHROM'][:]
chrom_res_seg = chrom.compress(res_seg_variants, axis=0)
chrom_res_seg

# %%
# some variants in 1000 genomes project have multiple variants at the same genomic position, 
# which causes problems for some selection tests in scikit-allel. 
# Let's check if there any of these.
count_multiple_variants = np.count_nonzero(np.diff(pos_res_seg == 0))

if count_multiple_variants == 0:
    print("No cases where there are multiple variants at the same genomic position, script will continue")
else:
    print("There are multiple variants at the same genomic position. This causes problems with some selection tests using sci-kit allel.")
    sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line


# %% ##################### SUSCEPTIBLE POPULATION #################

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
sus_seg_variants = ac_sus.is_segregating() & ac_sus.is_biallelic_01()
ac_sus_seg = ac_sus.compress(sus_seg_variants, axis=0)
gt_sus_seg = gt_sus_samples.compress(sus_seg_variants, axis = 0)
gt_sus_seg

# %%
## this is from a phased VCF so we can convert this genotype array to haplotype array

h_sus_seg = gt_sus_seg.to_haplotypes().compute()
h_sus_seg

# %% H12 was calculated using phased biallelic SNPs in 1000 bp windows
# SNP windows, using the garuds_h function in scikit-allel. 200 permutations?
# Calculate in 1000bp windows, look at the difference in H12
# Calculate for resistant samples

# A.miles:
# To calibrate the window sizes I ran the H12 scans with a range of different window sizes, and chose
# the smallest window size for which the mean value of H1 over all windows was below
# 0.01.
# Lucas et al (2023) to identify regions in which swept haplotypes are more frequent in resistant compared to susceptible individuals, they calculated
# the different in H12 value between groups, deltaH12.

h1, h12, h123, h2_h1 = allel.moving_garud_h(h_res_seg, 1000)

# Calculate for susceptible samples

h1, h12, h123, h2_h1 = allel.moving_garud_h(h_sus_seg, 1000)

# look for peaks of difference in H12 between resistant and susceptible samples


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