
######################## SPECIFIC SNP DETECTION #########################

# Do samples have the L995F mutation?
# using bash
# 
# bcftools query -i 'POS==2358254' -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE:%GT]\n' \
# -r 2L:2422652 F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz \
# > 2358254_genotypes.txt
# 
# bcftools query -i 'POS==2429745' -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE:%GT]\n' \
# -r 2L:2429745 F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz \
# > 2429745_genotypes.txt
# 

# %% using python
import subprocess

# Define the input file with positions (CHR and POS)
positions_file = "snppositions.txt"  # Replace with your file path
vcf_filename = "F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz"  # Replace with your VCF file path

# Define the output file
output_file = "genotype_summary.txt"

# Define the genotype options
genotypes = ["0/0", "0/1", "1/1"]

# Open the output file for writing
with open(output_file, "w") as outfile:
    # Read positions from the input file
    with open(positions_file, "r") as infile:
        First = True
        for line in infile:
            if First == True:
                First = False
                continue

            chr, pos = line.strip().split()
            
            # Define the bcftools command with a separate region filter for each position

            cmd = f'bcftools query -r {chr}:{pos} -f"%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE:%GT]\n" {vcf_filename}'
            
            # Execute the bcftools command
            try:
                result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True, check=True, shell=True)
                bcftools_output = result.stdout
                
                # Count genotypes
                genotype_counts = {genotype: bcftools_output.count(genotype) for genotype in genotypes}
                
                # Write results to the output file
                outfile.write(f"CHR: {chr}\tPOS: {pos}\n")
                for genotype, count in genotype_counts.items():
                    outfile.write(f"{genotype}: {count}\n")
                outfile.write("\n")
            except subprocess.CalledProcessError as e:
                print(f"Error processing CHR: {chr}, POS: {pos}. Skipping...")

print("Done")
