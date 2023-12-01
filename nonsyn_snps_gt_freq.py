
######## Identifying non-synonymous SNPs in VCF file in specific genes of interest ########
# %% Calculate genotype counts for sample set
import subprocess
import pandas as pd

# %% Define the input file with positions (CHR and POS)
vcf_filename = "F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz"  # Replace with your VCF file path

# %% Define the input and output files
positions_file = "genes_of_interest.txt"
output_file = "non_synonymous_snps.txt"

# %% Define the genotype options
genotypes_of_interest = {"0/1", "1/0", "1/1", "0|1", "1|0"}

# %% Open the output file for writing
with open(output_file, "w") as outfile:
    # Read positions from the input file
    with open(positions_file, "r") as infile:
        next(infile) # Skip header line
        for line in infile:
            # Split the line into columns
            columns = line.strip().split("\t")
            if len(columns) >= 4:
                gene, chr, start, stop = columns
            
            # Define the bcftools command using a region filter for each gene

            cmd = f'bcftools query -r {chr}:{start}-{stop} -f"%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE:%GT]\n" {vcf_filename}'
            
            # Execute the bcftools command
            try:
                result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True, check=True, shell=True)
                bcftools_output = result.stdout.splitlines()

                for line in bcftools_output:
                    fields = line.split('\t')
                    sample_genotypes = fields[4:]
                    
                    # Check if any of the genotypes of interest are present
                    if any(gt.split(':')[1] in genotypes_of_interest for gt in sample_genotypes):
                        outfile.write(line + '\n')

            except subprocess.CalledProcessError as e:
                    print(f"Error processing gene: {gene}, CHR: {chr}, range: {start}-{stop}. Error: {e}")

print("Done")
                
# %% Next take bcftools_ouput and make a table with the number of samples with each genotype for each of these positions

