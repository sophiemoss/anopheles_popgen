
## Filtered genotype VCF: this is where we pick the parameter for missingness

python /mnt/storage11/sophie/fastq2matrix/scripts/filter_merged_vcf.py \
--merged-file bi_snps_chr_gambiae_nov2022.2023_07_05.genotyped.vcf.gz \
--prefix gambiae_nov2022 \
--ref /mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.dna.toplevel.fa \
#--bqsr-vcf /mnt/storage9/emilia/Plasmodium_ref/Pf3D7_v3/3d7_hb3.combined.final.vcf.gz \
#--include-region /mnt/storage9/emilia/Plasmodium_ref/Pf3D7_v3/Core_genome_Pf3D7_v3_ext.bed \
--vqslod 0 \
--missing-sample-cutoff 0.2 \
--cutoff-mix-GT 0.8 \
--gff-file /mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.56.gff3 > filter_merged_vcf_log.txt 2>&1

## This output is piped to filter_merged_vcf_log.txt
## After each stage you get a checkpoint file produced which can be identified by the suffix, and is used in the analysis.

## Convert the filtered vcf to matrix
## to make biallelic matrix use .bi.GT.

python /mnt/storage9/sophie/fastq2matrix/scripts/vcf2matrix.py --vcf Pf_bijagos_Nov2021.bi.GT.miss0.4.vqslod.filt.snps.vcf.gz --no-iupacgt --threads 6

## to make multiallelic matrix use .GT.

python /mnt/storage9/sophie/fastq2matrix/scripts/vcf2matrix.py --vcf Pf_bijagos_Nov2021.GT.miss0.4.vqslod.filt.snps.vcf.gz --no-iupacgt --threads 6

## Creating a distance matrix for PCA
plink --vcf /mnt/storage9/sophie/basebijagospopgen/database/filtered/Pf_bijagos_Nov2021.bi.GT.miss0.4.vqslod.filt.snps.vcf.gz --distance square --double-id --allow-extra-chr --out Nov2021plusBijagos_dist_m --threads 6

#it will produce input files you need
#this should be from filtered vcf similarly as you wanted to do it from filtered matrix

