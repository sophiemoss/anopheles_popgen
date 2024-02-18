## Admixture

# Identify which samples are in your vcf:
bcftools query -l yourfile.vcf.gz

# STEP 1 EXTRA FILTERING
# MAF > 0.01 filter has already been applied

# STEP 2 CONVERT CHR NAMES IN VCF TO INTEGERS

zgrep -v "^#" FMISSING_MAF_AC0_DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.ann.vcf.gz | cut -f1 | sort | uniq

zcat FMISSING_MAF_AC0_DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.ann.vcf.gz | awk 'BEGIN{OFS=FS="\t"} /^#/ {print; next} {gsub(/^2L$/, "1", $1); gsub(/^2R$/, "2", $1); gsub(/^3L$/, "3", $1); gsub(/^3R$/, "4", $1); gsub(/^Mt$/, "5", $1); gsub(/^X$/, "6", $1); gsub(/^Y_unplaced$/, "7", $1); print}' | gzip > gambiae_malariagen_GB_GM-ABC_Bijagos_merged_renamedchr.vcf.gz

# STEP 3 MAKE BED AND BIM FILES

plink --vcf gambiae_malariagen_GB_GM-ABC_Bijagos_merged_renamedchr.vcf.gz --set-missing-var-ids @:# --keep-allele-order --const-fid --allow-extra-chr --make-bed --out gambiae_global_gambiaealigned

# STEP 4 RUN ADMIXTURE

cat K_runs.txt | xargs -I {} sh -c 'admixture --cv=10 -j20 melas_global_gambiaealigned.bed {} | tee log{}.cv10.seed12345.out'
