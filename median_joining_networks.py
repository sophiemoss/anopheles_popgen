# Subset VCF to be just vgsc gene

# start: 2,358,158 stop: 2,431,617

# for just Bijagos samples

bcftools view -r 2L:2358158-2431617 F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz -Oz -o VGSC_F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz

# subset reference to just be VGSC reference

samtools faidx Anopheles_gambiae.AgamP4.dna.toplevel.fa
samtools faidx Anopheles_gambiae.AgamP4.dna.toplevel.fa 2L:2358158-2431617 > VGSC_only_Anopheles_gambiae.AgamP4.dna.toplevel.fa
samtools faidx VGSC_only_Anopheles_gambiae.AgamP4.dna.toplevel.fa

# vcf2fasta.py for just the VGSC vcf file, using ref file.

python /mnt/storage11/sophie/fastq2matrix/scripts/vcf2fasta.py\
    --vcf VGSC_F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz \
        --ref /mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.dna.toplevel.fa \
            --threads 10 \
                --snps-no-filt

python /mnt/storage11/sophie/fastq2matrix/scripts/vcf2fasta_noiupac.py\
    --vcf VGSC_F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz \
        --ref /mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.dna.toplevel.fa \
            --threads 10 \
                --snps-no-filt


# for all samples, including Bijagos islands, The Gambia and mainland Guinea-Bissau

bcftools view -r 2L:2358158-2431617 gambiae_malariagen_GB_GM-ABC_Bijagos_merged.vcf.gz -Oz -o VGSC_gambiae_malariagen_GB_GM-ABC_Bijagos_merged.vcf.gz
tabix -p vcf VGSC_gambiae_malariagen_GB_GM-ABC_Bijagos_merged.vcf.gz

# vcf2fasta.py

python /mnt/storage11/sophie/fastq2matrix/scripts/vcf2fasta.py \
    --vcf /mnt/storage11/sophie/bijagos_mosq_wgs/malariagen_wgs/VGSC_gambiae_malariagen_GB_GM-ABC_Bijagos_merged.vcf.gz \
        --ref /mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.dna.toplevel.fa \
            --threads 10 \
                --snps-no-filt

# use these fasta files in the R scripts to make haplotype networks
# these fasta files are also needed for maximum likelihood trees (tree needed only for melas)

# for maximum likelihood tree use VCF file which is combined with malariagen

python /mnt/storage11/sophie/fastq2matrix/scripts/vcf2fasta.py \
    --vcf /mnt/storage11/sophie/bijagos_mosq_wgs/malariagen_wgs/gambiae_malariagen_GB_GM-ABC_Bijagos_merged.vcf.gz \
        --ref /mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.dna.toplevel.fa \
            --threads 10 \
                --snps-no-filt

# make fasta file from vcf for melas too

python /mnt/storage11/sophie/fastq2matrix/scripts/vcf2fasta.py \
    --vcf /mnt/storage11/sophie/bijagos_mosq_wgs/2019_melas_fq2vcf/genomics_database_melas2019plusglobal/genomics_database_melas2019plusglobal_vcf/melas2019plusglobal_vcf_filtering/F_MISSING_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_melas_2019_plusglobal.2023_07_05.genotyped.vcf.gz \
        --ref /mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.dna.toplevel.fa \
            --threads 10 \
                --snps-no-filt

###### trying with phased vcf now instead of converting to a fasta file. Otherwise I can use the --whole genome option.

# subset phased vcf to make phased vcf of just the VGSC gene


bcftools view -r 2L:2358158-2431617 2022gambiaevcfphased.vcf.gz -Oz -o VGSC_only_2022gambiae_phased.vcf.gz



# making fasta file with just VGSC reference and the --whole genome option

python /mnt/storage11/sophie/fastq2matrix/scripts/vcf2fasta_noiupac.py \
    --vcf /mnt/storage11/sophie/bijagos_mosq_wgs/malariagen_wgs/gambiae_malariagen_GB_GM-ABC_Bijagos_merged.vcf.gz \
        --ref /mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/VGSC_only_Anopheles_gambiae.AgamP4.dna.toplevel.fa \
            --threads 10 \
                --snps-no-filt \
                --whole-genome


