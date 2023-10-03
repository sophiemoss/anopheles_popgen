## Convert the filtered vcf to matrix
## to make biallelic matrix use .bi.GT.

python /mnt/storage11/sophie/fastq2matrix/scripts/vcf2matrix.py --vcf mt_melas_plus_global_subset_filtered.vcf.gz --no-iupacgt --threads 6

## to make multiallelic matrix use .GT.
# not making multiallelic matrix for anopheles, have already split multiallelic sites to be biallelic

## Creating a distance matrix for PCA
plink --vcf F_MISSING_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz --distance square --double-id --allow-extra-chr --out 2019_melas_dist_m --threads 6

## Create a distance matrix for PCA with melas plus global samples
plink --vcf F_MISSING_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_melas_2019_plusglobal.2023_07_05.genotyped.vcf.gz --distance square --double-id --allow-extra-chr --out 2019_melas_plus_global_dist_m --threads 6

## using just the mitochondria
## create mitochondria only vcf using bcftools, then use plink to make distance matrix
bcftools view -r Mt F_MISSING_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_melas_2019_plusglobal.2023_07_05.genotyped.vcf.gz -O z -o mito_only_melas_plusglobal.vcf.gz
# have to change the name Mt to something non-human so that plink will work and keep the variants in, usually it gets rid of all Mt variants
vim chrom_map.txt 
Mt  anop_mito
bcftools annotate --rename-chrs chrom_map.txt -O z -o mito_only_melas_plusglobal_renamed.vcf.gz mito_only_melas_plusglobal.vcf.gz
plink --vcf mito_only_melas_plusglobal_renamed.vcf.gz --distance square --double-id --out mito_only_melas_plusglobal_renamed --threads 6 --allow-extra-chr


# it will produce input files you need

# creating matrix for VGSC only VCF

python /mnt/storage11/sophie/fastq2matrix/scripts/vcf2matrix.py --vcf VGSC_only_F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz --no-iupacgt --threads 6

# Creating a distance matrix for VGSC PCA 
plink --vcf VGSC_only_F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz --distance square --double-id --allow-extra-chr --out 2019_melas_dist_m --threads 6


