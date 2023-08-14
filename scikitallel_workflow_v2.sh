## scikitallel_workflow

conda activate scikit

## Step 1: FST. To use scikit I need to create a genotype array.
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

### CALCULATING FST

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

# estimate Fst for each variant and each allele individually:

fst = a / (a+b+c)

# estimate Fst for each variant individually, averaging over alleles:

fst = (np.sum(a, axis=1) /
       (np.sum(a, axis=1) + np.sum(b, axis=1) + np.sum(c, axis=1)))

# estimate Fst averaging over all variants and alleles:

fst = np.sum(a) / (np.sum(a) + np.sum(b) + np.sum(c))

# Using array a to identify genomic positions with fst of certain value

# Find indices of values above 0.1
above_0_1_indices = np.where(a > 0.1)

# Get the genomic positions using the identified indices
genomic_positions = above_0_1_indices[0]

if len(genomic_positions) > 0:
    print("Genomic positions with values above 0.1:", genomic_positions)
else:
    print("No values above 0.1 found in the array.")

# Find indices of values above 0.15
above_0_15_indices = np.where(a > 0.15)

# Get the genomic positions using the identified indices
genomic_positions = above_0_15_indices[0]

if len(genomic_positions) > 0:
    print("Genomic positions with values above 0.15:", genomic_positions)
else:
    print("No values above 0.15 found in the array.")

## Seems a crazy low number of samples so using vcftools to calculate Fst too and comparing

vcftools --gzvcf F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz \
--weir-fst-pop resistant_samples.txt \
--weir-fst-pop susceptible_samples.txt \
--out fst_resistant_vs_susceptible

## after filtering only kept 33 out of 42 individuals
## there are plenty of higher Fst values using vcftools, so need to understand what is going wrong in the scikit allel workflow


## CALCULATE SELECTION STATISTICS

## SELECT SAMPLE POPULATION TO WORK WITH

# some selection tests only support biallelic variants, not multiallelic. 
# This filtered VCF should already be biallelic SNPs only.

### Need to use phased version for conversion to haplotype array for further selection analyses.
https://nbviewer.org/gist/alimanfoo/75b567b3d43810ef8eaef248b38b1c1c?flush_cache=true


## Make PCA for just genes of interest
## Use bcftools to subset the filtered vcf just for the genes of interest and make a PCA of that