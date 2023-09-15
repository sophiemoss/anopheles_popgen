

######################## ADMIXTURE #########################

# ## A.Miles thesis methods:
# The f3 and f4 tests were performed using a block size of 100,000 SNPs to estimate standard
# error via a block jackknife procedure (Patterson et al., 2012). Z scores reported in the
# results are computed by dividing the test statistic by the estimated standard error, and
# thus indicate the number of standard deviations from zero.


# Admixture
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
