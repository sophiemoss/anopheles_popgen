# FILTERING OF COMBINED GENOTYPED VCF
# same filters used for both gambiae pipeline and melas pipeline

# FILTER 1: remove INDELS with bcftools, here -M2 indicates that bcftools should split multi-allelic sites into multiple biallelic sites, 
# keeping this information
# -m2 is used in conjunction with -M2 to apply the minor-allele-based decomposition. 
# This means that bcftools will decompose multi-allelic sites using the minor allele as the reference allele in the biallelic split.
# here -v tells bcftools to only view SNPS, so indels are excluded

bcftools view -M2 -m2 -v snps gambiae_nov2022.2023_07_05.genotyped.vcf.gz -Oz -o bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf

# you must bgzip the file before indexing if you did not use -Oz to make the output bgzip
bgzip bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf

# tabix index the compressed VCF file, creates .vcf.gz.tbi
tabix -p vcf bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz

# filter out the contigs from the VCF file, note that to produce a properly bgzipped vcf file you need the -Oz flag

bcftools view bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz --regions 2L,2R,3L,3R,Mt,X,Y_unplaced | bcftools sort -Oz -o bi_snps_chr_gambiae_nov2022.2023_07_05.genotyped.vcf.gz

# view the unique chromosome names
bcftools query -f '%CHROM\n' bi_snps_chr_gambiae_nov2022.2023_07_05.genotyped.vcf.gz | sort | uniq > unique_chromosomes_filtered.txt

tabix -p vcf bi_snps_chr_gambiae_nov2022.2023_07_05.genotyped.vcf.gz

# FILTER 2: Filter samples to keep those with 4'0% of genome with > 10x coverage, and min-ac=1 so that all variants that remain are still variants after sample removal
# I have identified samples that are above this threshold using my basic WGS stats
# create file with samples to keep: sample_40_10.txt

bcftools view -S sample_40_10.txt bi_snps_chr_gambiae_nov2022.2023_07_05.genotyped.vcf.gz --min-ac=1 -Oz -o miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz

tabix -p vcf miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz

# FILTER 3: for GATK standard filtering
# GATK filter recommendations: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
# excellent explanation of this hard filtering and additional soft filtering
# https://evomics.org/learning/population-and-speciation-genomics/2020-population-and-speciation-genomics/first-steps-in-genomic-data-analysis/#ex3.4

gatk VariantFiltration \
-R /mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.dna.toplevel.fa \
-V miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz

## note there is a repeated warning here JEXL engine undefined variable, not a problem
## because this is coming from where there are positions with no coverage
## https://gatk.broadinstitute.org/hc/en-us/community/posts/4408733963803-GATK-Variant-Filtration-undefined-variable

tabix -p vcf gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz

## ADDITIONAL SOFT FILTERING

## After hard quality filtering, we have a data set containing only variant sites that we trust with high confidence. 
## However, so far we have only removed variants that had skewed values across all samples – 
## we haven’t yet looked at the individual genotypes at each variant site. 
## For example, we may have kept a variant that has a high quality across all of the individual samples, but there 
## may be some individuals where the read depth (FMT/DP) or the genotype quality (GQ) for the individual genotype is very low, 
## so we don’t trust the calls for these individuals

## FILTER 4: Filter reads that have a read depth below 5 OR a genotype quality below 20
## The difference is, that instead of the whole variant site, we are now considering single genotypes 
## (the information in the “FORMAT” fields”) using -e 'FMT/DP<3 | FMT/GQ<20' and we will not remove the whole variant site, 
## but only set the respective genotypes to missing (./.) by using the bcftools filter -S . command.
## Here we use the logical operator “|” instead of “||”, since using “||” would mean that every genotype at a variant site 
## is set to missing, even if only one genotype doesn’t fulfill the threshold

bcftools filter -S . -e 'FMT/DP<5 | FMT/GQ<20' -O z -o DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz

## check this, the genotype here should have been set to missing ./. as we are filtering for depth below 5

bcftools query -i 'FMT/DP<5' -f '[GT=%GT\tDP=%DP\tGQ=%GQ\t]\n' DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz | less -S

## happy? index vcf
tabix -p vcf DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz

## FILTER 5 exclude -e all sites at which no alternative alleles are called for any of the samples

bcftools filter -e 'AC==0' DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz -O z -o AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz

## How many variants are left?
bcftools view -H AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz | wc -l
21024383

## compare to previous
bcftools view -H DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz | wc -l
21679588

## compare to original
bcftools view -H gambiae_nov2022.2023_07_05.genotyped.vcf.gz | wc -l
36016533

## FILTER 6
## Remove variants with a high amount of missing genotypes and filter on minor allele frequency
## The data set overall will now have lots of missing data. 
## This is because we set genotypes with low quality and low read depth to missing. 
## Therefore we will now remove all variants that have more than 20% missing genotypes or MAF < 0.01
## Filtering for MAF 0.01 means that remaining sites in the VCF have a minimum minor allele frequency of 1%
## This dataset includes 44 individuals (thus maximally 88 alleles per site). If we filter for MAF 0.01, 
## a variant is excluded if its minor allele is found in less than 0.88 genotypes (1% of 88). 
## This would remove alleles where the minor allele occurs 0.88 times or less, which for this set means a 
## variant would be removed if there are no minor alleles? So if there is no alternate allele the variant is removed.

bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.01' -O z -o F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz

# check this has worked. The minor allele count should be 1 or above in this dataset.

bcftools query F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz -f'%AC\n' | sort -g | head

# Final number of variants (SNPS) in dataset

bcftools view -H F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz | wc -l

# gambiae 17085002
# melas   11121182 

tabix -p vcf F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz

# Filtering complete!

# Phase the filtered vcf using beagle, https://faculty.washington.edu/browning/beagle/beagle_5.2_13Oct21.pdf
# could also use phasing pipeline from malariagen. Have a look at this, and understand gatk parameters above.

mamba install beagle

beagle -Xmx500g gt=F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz out=2022gambiaevcfphased

##gambie error
1569759 now,
17085002 before

##melas running with 200g - phased is now 511726 lines

# 500g of memory is an insane amount, half the entire server, for some reason it would not run without specifying this
# next time try to use less but also make sure nobody needs it if using
# it did only take 14 minutes to run but it has removed all annotations and it smaller by a factor of 10
# spoken to Joe and Emilia, neither of them have used Beagle or know why it is doing this. Speak to Jody. 