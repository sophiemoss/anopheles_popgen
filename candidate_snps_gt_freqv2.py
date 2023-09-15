
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

# %% Calculate genotype counts for sample set
import subprocess
import matplotlib

# Define the input file with positions (CHR and POS)
positions_file = "snp_positions.txt"  # Replace with your file path
vcf_filename = "/mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf_agamP4/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022_combinedvcf_filteringsteps/F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz"  # Replace with your VCF file path

# Define the output file
output_file = "genotype_summary_nosamplenames.txt"

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

            # Split the line into columns
            columns = line.strip().split("\t")

            # Assuming CHR is in the second column and POS is in the first column
            if len(columns) >= 2:
                chr = columns[1]
                pos = columns[0]
            
            # Define the bcftools command with a separate region filter for each position

            cmd = f'bcftools query -r {chr}:{pos} -f"%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE:%GT]\n" {vcf_filename}'
            
            # Execute the bcftools command
            try:
                result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True, check=True, shell=True)
                bcftools_output = result.stdout
                
                # Count genotypes
                genotype_counts = {genotype: bcftools_output.count(genotype) for genotype in genotypes}
                
                # Check if any genotype count is greater than 0
                if any(count > 0 for count in genotype_counts.values()):
                    # Write results to the output file
                    outfile.write(f"CHR: {chr}\tPOS: {pos}\n")
                    for genotype, count in genotype_counts.items():
                        outfile.write(f"{genotype}: {count}\n")
                    outfile.write("\n")
                else:
                    print(f"CHR: {chr}\tPOS: {pos} - Variant position not in VCF, no genotypes found. Skipping...")
            
            except subprocess.CalledProcessError as e:
                print(f"Error processing CHR: {chr}, POS: {pos}. Skipping...")

print("Done")

# %% Plot this on a bar chart.

import matplotlib.pyplot as plt

# Reading the file
with open("genotype_summary_nosamplenames.txt", "r") as file:
    lines = file.readlines()

# Parsing the data
positions = []
genotypes = {"0/0": [], "0/1": [], "1/1": []}
for line in lines:
    if "POS:" in line:
        pos = int(line.split()[-1])
        positions.append(pos)
    elif "0/0" in line or "0/1" in line or "1/1" in line:
        genotype, count = line.strip().split(':')
        genotypes[genotype].append(int(count))

# Plotting the data
bar_width = 0.3
index = range(len(positions))

fig, ax = plt.subplots()
bar1 = ax.bar(index, genotypes["0/0"], bar_width, label='0/0')
bar2 = ax.bar([i + bar_width for i in index], genotypes["0/1"], bar_width, label='0/1')
bar3 = ax.bar([i + bar_width*2 for i in index], genotypes["1/1"], bar_width, label='1/1')

# Labeling and formatting the chart
ax.set_xlabel('POS')
ax.set_ylabel('Count')
ax.set_title('Genotype Counts by Position')
ax.set_xticks([i + bar_width for i in index])
ax.set_xticklabels(positions, rotation=45)
ax.legend()

plt.tight_layout()
plt.show()


# %%
######### Print genotypes with sample names #####

# Define the input file with positions (CHR and POS)
positions_file = "snp_positions.txt"
vcf_filename = "/mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf_agamP4/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022_combinedvcf_filteringsteps/F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz"

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

            # Split the line into columns
            columns = line.strip().split("\t")

            # Assuming CHR is in the second column and POS is in the first column
            if len(columns) >= 2:
                chr = columns[1]
                pos = columns[0]

            # Define the bcftools command with a separate region filter for each position
            cmd = f'bcftools query -r {chr}:{pos} -f"%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE:%GT]\n" {vcf_filename}'

            # Execute the bcftools command
            try:
                result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True, check=True, shell=True)
                bcftools_output = result.stdout

                sample_data = bcftools_output.strip().split()[4:]  # Get the sample data starting from the fifth column
                samples_per_genotype = {genotype: [] for genotype in genotypes}

                for s_data in sample_data:
                    if ":" in s_data:
                        sample_name, genotype = s_data.split(":")
                        if genotype in samples_per_genotype:
                            samples_per_genotype[genotype].append(sample_name)

                outfile.write(f"CHR: {chr}\tPOS: {pos}\n")
                for genotype, samples in samples_per_genotype.items():
                    outfile.write(f"{genotype}: {', '.join(samples)}\n")
                outfile.write("\n")

            except subprocess.CalledProcessError as e:
                print(f"Error processing CHR: {chr}, POS: {pos}. Skipping...")

print("Done")

# %% Make a text file which gives you the breakdown of which samples contain 'res' in the sample name,
# which contain 'pos' and which do not contain either so are control 

import subprocess

positions_file = "snp_positions.txt"
vcf_filename = "/mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf_agamP4/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022_combinedvcf_filteringsteps/F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz"  # Replace with your VCF file path
output_file = "genotype_by_phenotype.txt"
genotypes = ["0/0", "0/1", "1/1"]

with open(output_file, "w") as outfile:
    with open(positions_file, "r") as infile:
        first = True
        for line in infile:
            if first:
                first = False
                continue

            columns = line.strip().split("\t")

            if len(columns) >= 2:
                chr = columns[1]
                pos = columns[0]

            cmd = f'bcftools query -r {chr}:{pos} -f"%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE:%GT]\n" {vcf_filename}'

            try:
                result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True, check=True, shell=True)
                bcftools_output = result.stdout

                # Create a nested dictionary to count genotypes per category
                counts = {category: {genotype: 0 for genotype in genotypes} for category in ["res", "sus", "other"]}
                
                lines = bcftools_output.strip().split("\n")
                if lines:
                    for output_line in lines:
                        fields = output_line.split("\t")
                        for field in fields[4:]:
                            sample, genotype = field.split(":")
                            if genotype in genotypes:
                                if 'res' in sample:
                                    counts['res'][genotype] += 1
                                elif 'sus' in sample:
                                    counts['sus'][genotype] += 1
                                else:
                                    counts['other'][genotype] += 1

                # Check if we have any counts to report
                if any(count > 0 for category in counts.values() for count in category.values()):
                    outfile.write(f"CHR: {chr}\tPOS: {pos}\n")
                    for category, genotype_counts in counts.items():
                        for genotype, count in genotype_counts.items():
                            outfile.write(f"{category}_{genotype}: {count}\n")
                    outfile.write("\n")
                else:
                    print(f"CHR: {chr}\tPOS: {pos} - Variant position not in VCF, no genotypes found. Skipping...")

            except subprocess.CalledProcessError as e:
                print(f"Error processing CHR: {chr}, POS: {pos}. Error: {e}")

print("Done")



# %%
