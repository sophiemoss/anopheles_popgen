# Subset VCF to be just vgsc gene

# start: 2,358,158 stop: 2,431,617

# for just Bijagos samples

bcftools view -r 2L:2358158-2431617 F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz -Oz -o VGSC_F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz

# vcf2fasta.py

python /mnt/storage11/sophie/fastq2matrix/scripts/vcf2fasta.py --vcf VGSC_F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz --ref /mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.56.chr.gff3 --snps-no-filt --threads 10

# for samples including The Gambia and Guinea-Bissau

bcftools view -r 2L:2358158-2431617 gambiae_malariagen_GB_GM-ABC_Bijagos_merged.vcf.gz -Oz -o VGSC_gambiae_malariagen_GB_GM-ABC_Bijagos_merged.vcf.gz
tabix -p vcf VGSC_gambiae_malariagen_GB_GM-ABC_Bijagos_merged.vcf.gz

# vcf2fasta.py

python /mnt/storage11/sophie/fastq2matrix/scripts/vcf2fasta.py --vcf VGSC_gambiae_malariagen_GB_GM-ABC_Bijagos_merged.vcf.gz --ref /mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.56.chr.gff3 --snps-no-filt --threads 10

# use these fasta files in the R scripts to make haplotype networks
# these fasta files are also needed for maximum likelihood trees