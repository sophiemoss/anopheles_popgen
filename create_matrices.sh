## Convert the filtered vcf to matrix
## to make biallelic matrix use .bi.GT.

python /mnt/storage11/sophie/fastq2matrix/scripts/vcf2matrix.py --vcf F_MISSING_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz --no-iupacgt --threads 6

## to make multiallelic matrix use .GT.
# not making multiallelic matrix for anopheles, have already split multiallelic sites to be biallelic

## Creating a distance matrix for PCA
plink --vcf F_MISSING_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz --distance square --double-id --allow-extra-chr --out 2019_melas_dist_m --threads 6

# it will produce input files you need

# creating matrix for VGSC only VCF

python /mnt/storage11/sophie/fastq2matrix/scripts/vcf2matrix.py --vcf VGSC_only_F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz --no-iupacgt --threads 6

# Creating a distance matrix for VGSC PCA 
plink --vcf VGSC_only_F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz --distance square --double-id --allow-extra-chr --out 2019_melas_dist_m --threads 6


