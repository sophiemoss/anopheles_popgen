
######################## SELECTION STATISTICS #########################

# some selection tests only support biallelic variants, not multiallelic. 
# This filtered VCF should already be biallelic SNPs only.
# Note that the zarr file should have been made from a phased vcf
# You make the zarr file with allel.vcf_to_zarr('phased_vcf_file_name.vcf.gz', 'output_name.zarr', fields='*', overwrite=True)

# %%
import os
os.chdir('/mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf_agamP4/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022_filtering')
os.getcwd()

# %%
import numpy as np
import allel
import zarr
import pandas as pd
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import random

## convert phased, filtered, VCF file to zarr file
# %%
# allel.vcf_to_zarr('2022gambiaevcfphased.vcf.gz', '2022gambiaevcfphased.zarr', fields='*', overwrite=True)

# %%
callset = zarr.open('2022_gambiae.zarr', mode='r')
#callset.tree(expand=True)

# %%
## convert zarr file to genotype array (GenotypeDaskArray)
genotype_array = allel.GenotypeDaskArray(callset['calldata/GT'])
print(genotype_array.shape)

# %%
## import metadata
df_samples=pd.read_csv('metadata_gambiae_2022.csv',sep=',',usecols=['sample','year','country','island','phenotype'])
df_samples.head()
df_samples.groupby(by=['phenotype']).count

# %% Filter genotype array
# 1. Make allele count array so that we can filter for biallelic and segregating variants
ac_array = genotype_array.count_alleles(max_allele=8).compute() # ac_array is an AlleleCountsArray

# %% Filter for biallelic and segregating and store as boolean array
ac_array_filter = ac_array.is_biallelic_01()

# %% remove variants that are on Y_unplaced chromosome
chrom = callset['variants/CHROM'][:]
exclude_chrom = 'Y_unplaced'
ac_array_filter = ac_array_filter & (chrom != exclude_chrom)

# %% filter the allele counts array using the filter
filtered_ac_array = ac_array.compress(ac_array_filter, axis=0)

# %% also make a genotype dask array
filtered_gt = genotype_array.compress(ac_array_filter, axis = 0)

# %% partition samples by phenotype using metadata, split by index value and store as array
res_samples = df_samples[df_samples['phenotype'] == 'resistant'].index.values
sus_samples = df_samples[df_samples['phenotype'] == 'susceptible'].index.values
notres_samples = df_samples[df_samples['phenotype'] != 'resistant'].index.values
all_samples = df_samples.index.values

# %% select genotypes for variants and store as a genotype dask array
# this array contains 23 columns (1 for each mosquito), 16 million rows for each variant, and stores the genotype for each.
gt_res_samples = filtered_gt.take(res_samples, axis=1)
gt_sus_samples = filtered_gt.take(sus_samples, axis=1)
gt_all_samples = filtered_gt.take(all_samples, axis=1)
gt_notres_samples = filtered_gt.take(notres_samples, axis=1)

# %% convert this genotype array to haplotype array (we can do this because the original data was phased)
# the haplotype array is similar to the genotype array, but there are two columns per mosquito, one for each haplotype
h_res_seg = gt_res_samples.to_haplotypes().compute()
h_sus_seg = gt_sus_samples.to_haplotypes().compute()
h_all_seg = gt_all_samples.to_haplotypes().compute()
h_notres_seg = gt_notres_samples.to_haplotypes().compute()

# %% also store chromosome of each variant as we need this for shading the plots later
chrom = callset['variants/CHROM'][:]
chrom_filtered = chrom.compress(ac_array_filter, axis=0)

# %% we need variant positions
pos = callset['variants/POS'][:]
pos_filtered = pos.compress(ac_array_filter, axis=0)

# %% some variants in 1000 genomes project have multiple variants at the same genomic position, 
# which causes problems for some selection tests in scikit-allel. 
# Let's check if there any of these.
count_multiple_variants = np.count_nonzero(np.diff(pos_filtered== 0))

if count_multiple_variants == 0:
    print("No cases where there are multiple variants at the same genomic position, script will continue")
else:
    print("There are multiple variants at the same genomic position. This causes problems with some selection tests using sci-kit allel.")
    sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line

# %% H12 was calculated using phased biallelic SNPs in 1000 bp windows along the genome
# SNP windows, using the garuds_h function in scikit-allel.
# %% Calculate h values for resistant and susceptible samples
real_res_h1, real_res_h12, real_res_123, real_res_h2_h1 = allel.moving_garud_h(h_res_seg, 1000)

real_sus_h1, real_sus_h12, real_sus_h123, real_sus_h2_h1 = allel.moving_garud_h(h_sus_seg, 1000)

real_notres_h1, real_notres_h12, real_notres_123, real_notres_h2_h1 = allel.moving_garud_h(h_notres_seg, 1000)

real_all_h1, real_all_h12, real_all_123, real_all_h2_h1 = allel.moving_garud_h(h_all_seg, 1000
                                                                                           
# %% The above variables are stored as numpy.nd arrays

max_real_res_h12 = np.max(real_res_h12)
print("Maximum res H12 value:", max_real_res_h12)

max_real_sus_h12 = np.max(real_sus_h12)
print("Maximum sus H12 value:", max_real_sus_h12)

max_real_notres_h12 = np.max(real_notres_h12)
print("Maximum not resistant H12 value:", max_real_notres_h12)

########## Plotting ##########

# %% Check the number of windows for each array 
num_windows_res = real_res_h12.shape[0]
num_windows_sus = real_sus_h12.shape[0]
num_windows_notres = real_notres_h12.shape[0]

print("Number of windows for resistant samples:", num_windows_res)
print("Number of windows for not resistant:", num_windows_notres)
print("Number of windows for susceptible:", num_windows_sus)

# %% The h12 values are calculated in windows of 1000 SNPs. Each SNP has a POS value which is in the POS array.
window_size = 1000  # Define the window size as 1000 SNPs
num_windows = len(pos_filtered) // window_size  # Calculate the number of windows in total

# %% Calculate the median genomic position for each window
# for each iteration (i) a segment of pos_res_seg is processed
median_positions = [np.median(pos_filtered[i * window_size: (i + 1) * window_size]) for i in range(num_windows)]

# %% Plot h12 values for each chromosome
# Create a dictionary to hold H12 values and positions for each chromosome
chrom_h12_data = {}
# Loop through each window
for i in range(num_windows):
    chrom = chrom_filtered[i * window_size]  # Assumes the chromosome for all SNPs in a single window is consistent 
    pos = median_positions[i]
    h12 = real_res_h12[i]
    # Add this data to the corresponding chromosome in the dictionary
    if chrom not in chrom_h12_data:
        chrom_h12_data[chrom] = {'positions': [], 'h12': []}
    chrom_h12_data[chrom]['positions'].append(pos)
    chrom_h12_data[chrom]['h12'].append(h12)
# Now plot for each chromosome
for chrom in chrom_h12_data:
    plt.figure(figsize=(10, 6))
    plt.scatter(chrom_h12_data[chrom]['positions'], chrom_h12_data[chrom]['h12'], alpha=0.6)
    plt.xlabel('Median position of SNP windows across the chromosome')
    plt.ylabel('H12 Value in Resistant Samples')
    plt.title(f'Genome Wide Selection Scan of H12 Values across Chromosome {chrom}')
    plt.show()

# %% Plot per chromosome for susceptible samples
# Create a dictionary to hold H12 values and positions for each chromosome
chrom_h12_data = {}
# Loop through each window
for i in range(num_windows):
    chrom = chrom_filtered[i * window_size]  # Assuming each window is from the same chromosome
    pos = median_positions[i]
    h12 = real_sus_h12[i]  
    # Add this data to the corresponding chromosome in the dictionary
    if chrom not in chrom_h12_data:
        chrom_h12_data[chrom] = {'positions': [], 'h12': []}
    chrom_h12_data[chrom]['positions'].append(pos)
    chrom_h12_data[chrom]['h12'].append(h12)
# Now plot for each chromosome
for chrom in chrom_h12_data:
    plt.figure(figsize=(10, 6))
    plt.scatter(chrom_h12_data[chrom]['positions'], chrom_h12_data[chrom]['h12'], alpha=0.6)
    plt.xlabel('Median position of SNP windows across the Chromosome')
    plt.ylabel('H12 Value in Susceptible Samples')
    plt.title(f'Genome Wide Selection Scan of H12 Values across Chromosome {chrom}')
    plt.show()

# %% Plot resistant and susceptible samples on the same graph
# Create a dictionary to hold H12 values and positions for each chromosome, for both resistant and susceptible
chrom_h12_data = {}
# Loop through each window
for i in range(num_windows):
    chrom = chrom_filtered[i * window_size]  # Assuming each window is from the same chromosome
    pos = median_positions[i]
    h12_res = real_res_h12[i]
    h12_sus = real_sus_h12[i]
    
    # Add this data to the corresponding chromosome in the dictionary
    if chrom not in chrom_h12_data:
        chrom_h12_data[chrom] = {'positions': [], 'h12_res': [], 'h12_sus': []}
    chrom_h12_data[chrom]['positions'].append(pos)
    chrom_h12_data[chrom]['h12_res'].append(h12_res)
    chrom_h12_data[chrom]['h12_sus'].append(h12_sus)

# Now plot for each chromosome
for chrom in chrom_h12_data:
    plt.figure(figsize=(10, 6))
    
    # Plotting resistant and susceptible H12 values
    plt.scatter(chrom_h12_data[chrom]['positions'], chrom_h12_data[chrom]['h12_res'], alpha=0.6, color='blue', label='Resistant')
    plt.scatter(chrom_h12_data[chrom]['positions'], chrom_h12_data[chrom]['h12_sus'], alpha=0.6, color='plum', label='Susceptible')
    
    plt.xlabel('Median position of SNP windows across the chromosome')
    plt.ylabel('H12 Value')
    plt.title(f'Genome Wide Selection Scan of H12 Values across Chromosome {chrom}')
    plt.legend()
    plt.show()

# %% Plot resistant and non-resistant (similar sample sizes) on the same graph
# Create a dictionary to hold H12 values and positions for each chromosome, for both resistant and susceptible
chrom_h12_data = {}
# Loop through each window
for i in range(num_windows):
    chrom = chrom_filtered[i * window_size]  # Assuming each window is from the same chromosome
    pos = median_positions[i]
    h12_res = real_res_h12[i]
    h12_notres = real_notres_h12[i]
    
    # Add this data to the corresponding chromosome in the dictionary
    if chrom not in chrom_h12_data:
        chrom_h12_data[chrom] = {'positions': [], 'h12_res': [], 'h12_notres': []}
    chrom_h12_data[chrom]['positions'].append(pos)
    chrom_h12_data[chrom]['h12_res'].append(h12_res)
    chrom_h12_data[chrom]['h12_notres'].append(h12_notres)

# Now plot for each chromosome
for chrom in chrom_h12_data:
    plt.figure(figsize=(10, 6))
    
    # Plotting resistant and susceptible H12 values
    plt.scatter(chrom_h12_data[chrom]['positions'], chrom_h12_data[chrom]['h12_res'], alpha=0.6, color='blue', label='Resistant')
    plt.scatter(chrom_h12_data[chrom]['positions'], chrom_h12_data[chrom]['h12_notres'], alpha=0.6, color='plum', label='Not resistant')
    
    plt.xlabel('Median position of SNP windows across the chromosome')
    plt.ylabel('H12 Value')
    plt.title(f'Genome Wide Selection Scan of H12 Values across Chromosome {chrom}')
    plt.legend()
    plt.show()

# %% Identify H12 values over 0.2 for the peaks that we have seen, for each chromosome

from collections import Counter
# Create a dictionary to hold median positions and H12 values for each chromosome
chrom_h12_data = {}
# Loop through each window
for i in range(num_windows):
    # Extract the subset of chromosomes and positions for this window
    window_chroms = chrom_filtered[i * window_size: (i + 1) * window_size]
    window_pos = pos_filtered[i * window_size: (i + 1) * window_size]
    # Determine the most common chromosome in this window
    chrom_counter = Counter(window_chroms)
    common_chrom, _ = chrom_counter.most_common(1)[0]
    # Calculate the median position for this window
    pos_median = np.median(window_pos)
    h12_res = real_res_h12[i]  # Replace with the relevant H12 array
    # Filter based on H12 value and store data
    if h12_res > 0.2:
        if common_chrom not in chrom_h12_data:
            chrom_h12_data[common_chrom] = {'positions': [], 'h12_values': []}
        chrom_h12_data[common_chrom]['positions'].append(pos_median)
        chrom_h12_data[common_chrom]['h12_values'].append(h12_res)

# %% Find maximum h12 values and the median positions of those windows, for each chromosome. using real_all_h12

# Initialize a dictionary to store the maximum H12 value and corresponding median position for each chromosome
max_h12_data = {}

# Loop through each window
for i in range(num_windows):
    chrom = chrom_filtered[i * window_size]  # Assumes the chromosome for all SNPs in a single window is consistent 
    pos_median = np.median(pos_filtered[i * window_size: (i + 1) * window_size])
    h12 = real_all_h12[i]

    # Check if this chromosome is already in the dictionary
    if chrom not in max_h12_data:
        max_h12_data[chrom] = {'max_h12': h12, 'median_pos': pos_median}
    else:
        # Update the maximum H12 value and corresponding position if the current H12 is greater
        if h12 > max_h12_data[chrom]['max_h12']:
            max_h12_data[chrom]['max_h12'] = h12
            max_h12_data[chrom]['median_pos'] = pos_median

# Print the results
for chrom, data in max_h12_data.items():
    print(f"Chromosome {chrom} has maximum H12 value {data['max_h12']} at median position {data['median_pos']}")

# %%  Find H12 values above 0.2 in each chromosome and find the SNP windows that these correspond to
# Get the median SNP positions for those windows
# Define window size (assuming you know the size of each window)
window_size = 1000  # Example window size
# Open a file to write outside the loop
with open('h12_peak_median_positions.csv', 'w') as file:
    # Write header line
    file.write("Chromosome,Start,End\n")
    # Loop through each chromosome in the dictionary
    for chrom, data in chrom_h12_data.items():
        for position in data['positions']:
            # Calculate start and end positions based on the median
            start = position - window_size // 2
            end = position + window_size // 2
            # Create the line to be written in CSV format
            line = f"{chrom},{start},{end}\n"
            # Write to the file
            file.write(line)
            # Optionally, print the line as well
            print(line)

# %% Make a range for each peak that you want to look at genes beneath
Write code to automate this section

# Loop through each line of the h12_peaks_locations.txt file
while IFS=',' read -r chr start end; do
    echo "Processing: $chr from $start to $end" # Debugging line
    # Use awk to filter genes within the window from the GFF file
    awk -v chr="$chr" -v start="$start" -v end="$end" '
        $1 == chr && $4 <= end && $5 >= start {print $0}' \
        /mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.56.chr.gff3 >> h12_peaks_genes.txt
done < chromosome_ranges.csv

# %% Try to do the above in python

