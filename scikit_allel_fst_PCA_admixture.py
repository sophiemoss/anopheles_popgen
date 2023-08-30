## scikitallel_workflow

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

######################### CALCULATE Fst ##############################

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

