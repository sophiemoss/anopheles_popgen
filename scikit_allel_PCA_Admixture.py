

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
