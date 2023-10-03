# making a maximum likelihood tree containing my melas samples to look at where they fall amongst other species
# include my melas and my anopheles gambiae sensu stricto

# 1. Split the genomic VCF to be mitochondrial only, and then convert the VCF to fasta using --whole genome option
# mitochondrial position on AgamP4_Mt vectorbase 1:15363 

# Create a list of sample names
ls *.g.vcf.gz | sed 's/.g.vcf.gz//' > samples.txt

# Loop through the sample names in samples.txt
while read -r sample; do
    # Extract only mitochondrial positions and save as a mitochondrial-only VCF
    bcftools view -r Mt "${sample}.g.vcf.gz" -O z -o "${sample}_mito_only.g.vcf.gz"
    tabix -p vcf "${sample}_mito_only.g.vcf.gz"

    # Convert to fasta file
    python /mnt/storage11/sophie/fastq2matrix/scripts/vcf2fasta.py \
        --vcf "${sample}_mito_only.g.vcf.gz" \
        --ref "/mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Mt_only_Anopheles_gambiae.AgamP4.dna.toplevel.fa" \
        --threads 10 \
        --whole-genome \
        --snps-no-filt >> "${sample}.log" 2>&1
done < samples.txt

# troubleshooting for a single sample

bcftools view -r Mt "NG-33833_5xres2_lib708244_10265_1.g.vcf.gz" -O z -o "NG-33833_5xres2_lib708244_10265_1_mito_only.g.vcf.gz"

tabix -p vcf "NG-33833_5xres2_lib708244_10265_1_mito_only.g.vcf.gz"

python /mnt/storage11/sophie/fastq2matrix/scripts/vcf2fasta.py \
        --vcf "NG-33833_5xres2_lib708244_10265_1_mito_only.g.vcf.gz" \
        --ref "/mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Mt_only_Anopheles_gambiae.AgamP4.dna.toplevel.fa" \
        --threads 10 \
        --snps-no-filt \
        --whole-genome
#

# Take fa files and combine them with fa files from ncbi to make big maximum likelihood tree
# Align using muscle
# Trim using aliview
# Make tree
# View in iTOL