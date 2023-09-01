## scikitallel_workflow

conda activate scikit

######################## CALCULATING FST #########################
## Set wd
# %%
import os
os.chdir('/mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf_agamP4/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022_combinedvcf_filteringsteps')
os.getcwd()

# %% CONVERT FILTERED UNPHASED VCF TO ZARR FILE
# Try to convert to zarr and use workflow here: https://nbviewer.org/gist/alimanfoo/75b567b3d43810ef8eaef248b38b1c1c?flush_cache=true

import numpy as np
import allel
import zarr
import numcodecs
import pandas as pd

# %% Note this is using unphased VCF to convert to zarr at the moment
# allel.vcf_to_zarr('example.vcf', 'example.zarr', fields='*', overwrite=True)
# print('zarr', zarr.__version__, 'numcodecs', numcodecs.__version__)

callset = zarr.open('F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.zarr', mode='r')
callset.tree(expand=True)

# %%  CONVERT ZARR FILE TO GENOTYPE ARRAY

#gt = allel.GenotypeDaskArray(callset['calldata/GT'])
#print(gt.shape)
#gt

genotype_all = allel.GenotypeChunkedArray(callset['calldata/GT'])
genotype_all

# %%  IMPORT METADATA

df_samples=pd.read_csv('metadata_gambiae_2022.csv',sep=',',usecols=['sample','year','country','island','phenotype'])
df_samples.head()
df_samples.groupby(by=['phenotype']).count()

# %%  obtain indices of resistant samples in the dataset
resistant_samples_idx = df_samples[df_samples['phenotype'] == 'resistant'].index.values
resistant_samples_idx

susceptible_samples_idx = df_samples[df_samples['phenotype'] == 'susceptible'].index.values
susceptible_samples_idx

control_samples_idx = df_samples[df_samples['phenotype'] == 'control'].index.values
control_samples_idx

# %%  Define subpopulations
subpops = [susceptible_samples_idx, resistant_samples_idx]

# %% allele counts

# acs = genotype_all.count_alleles()
acs = genotype_all.count_alleles_subpops(subpops)

######################### CALCULATE Fst ##############################
# %% 
# I want to analyse the component of variance between populations
# blen=None is not working, I need to specify a value for blen
# a is the component of variance between populations
# b is the component of variance between individuals within populations
# c is the component of variance between gametes within individuals

# in simpler terms

# a: How much genetic diversity is there within a single group of individuals from the same population?
# b: How different are the genetic traits between different groups of individuals (populations)?
# c: Do certain genetic traits tend to appear together in similar patterns across different populations?

# %% Calcualte Fst

a,b,c = allel.weir_cockerham_fst(gt, subpops, max_allele=1, blen=1000)

# %%
snp_fst_wc = (a / (a + b + c))[:, 0]
snp_fst_wc

# %% 
# estimate Fst for each variant and each allele individually:

fst = a / (a+b+c)

# %% 
# estimate Fst for each variant individually, averaging over alleles:

fst = (np.sum(a, axis=1) /
       (np.sum(a, axis=1) + np.sum(b, axis=1) + np.sum(c, axis=1)))

# %% 
# estimate Fst averaging over all variants and alleles:

fst = np.sum(a) / (np.sum(a) + np.sum(b) + np.sum(c))

# %% 
# Create genomic position array

import dask.array as da
genomic_positions = da.arange(gt.shape[0])

# %% 
# create a mapping between genomic positions and Fst values

import pandas as pd

# %% 
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




