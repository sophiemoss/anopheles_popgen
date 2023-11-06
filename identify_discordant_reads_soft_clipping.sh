#Put this in a jupyter notebook?

# Identifying possible CNVS from Lucas et al (2023) copy number breakpoints

# We will look for soft clipped reads in the bam files at the breakpoints indicated in the supplementary information in this paper

# create list of bamfiles
# ls *_1.mkdup_sorted.bam > bam_files.txt

# 1. Sort bam files 
samtools sort AJ0126-C.bam > AJ0126-C.sorted.bam

# 2. Create coverage file from sorted bam file using samtools depth
samtools depth -aa AJ0126-C.sorted.bam > AJ0126-C.sorted.coverage

# 3. Create sam file from sorted bam,
samtools view -h AJ0126-C.sorted.bam -o AJ0126-C.sorted.sam


## Example loop for one of these commands

while read -r bamfile; do
    # Generate the output file name by appending "_filtered" before the ".bam" extension
    output="${bamfile%.*}.sam"
    
    # Apply the samtools command
    samtools view -h "$bamfile" > "$output"
    
    echo "Processed $bamfile -> $output"
done < bam_files.txt


# Make a loop to do all of this in one go:
# Sort bam files to make sorted bam file for each sample
# Create coverage file from sorted bam file using samtools depth
# Create sam file from sorted bam using samtools view

while read -r bamfile; do
    # Generate the output file name by appending "_filtered" before the ".bam" extension
    output="${bamfile%.*}_sorted.bam"
    
    # Apply the samtools command
    samtools sort "$bamfile" > "$output"
    
    echo "Processed $bamfile -> $output by sorting bamfile"

    # Next generate coverage file from the sorted BAM
    coverage_file="${output%.*}.coverage"

    # Apply the second samtool command

    samtools depth -aa "$output" > "$coverage_file"

    echo "Coverage generated for $output -> $coverage_file"

    # Apply the third samtools command
    samfile="${bamfile%.*}.sam"
    samtools view -h "$output" > "$samfile"

    echo "Sam file generated for $output > $samfile"

done < bam_files.txt

# Run the clip_identifier python script

python ClipIdentifier.py /Path/To/Your/SamFile /Path/To/Your/CoverageFile
python /mnt/storage11/sophie/gitrepos/anophelesmelas_popgen/ClipIdentifier.py AJ0126-C.sorted.primary.sam AJ0126-C.sorted.coverage

# loop through samples with clipidentifier.py script

ls *.sam | sed 's/.sam//' > samples.txt

while read -r sample; do
    # Run the ClipIdentifier.py script
    python /mnt/storage11/sophie/gitrepos/anophelesmelas_popgen/ClipIdentifier.py "$sample.sam" "$sample_sorted.coverage"   
    echo "Processed $sample with ClipIdentifier.py"
done < samples.txt

####### CREATE CSV OF SOFT-CLIPPING ##########

import pandas as pd

# Load the discordant read guide
discordant_df = pd.read_csv('discordant_read_guide.csv')

# Group the guide by Duplication_ID
grouped_discordant_df = discordant_df.groupby('Duplication_ID')

# Read the sample names from a file
with open('samples.txt', 'r') as file:
    samples = file.read().splitlines()

# Initialize an empty dictionary to hold the results
results_dict = {}

# Iterate over the grouped DataFrame
for dup_id, group in grouped_discordant_df:
    results_dict[f"{dup_id}_Pos_Range_Start"] = {}
    results_dict[f"{dup_id}_Pos_Range_End"] = {}

    for sample in samples:
        # Initialize sum for each sample for start and end positions
        results_dict[f"{dup_id}_Pos_Range_Start"][sample] = 0
        results_dict[f"{dup_id}_Pos_Range_End"][sample] = 0

        # Load the clipping data for the sample
        clipping_df = pd.read_csv(f'{sample}_sorted.Clipping.Normalised.csv')

        # Sum the NormalisedClipping for the start range
        start_group = group[group['Type'] == 'Start']
        if not start_group.empty:
            start_contig = start_group['Contig'].values[0]
            start_range = start_group[['Pos_Range_Start', 'Pos_Range_End']].values[0]
            start_sum = clipping_df[(clipping_df['Contig'] == start_contig) &
                                    (clipping_df['ClipPos'] >= start_range[0]) &
                                    (clipping_df['ClipPos'] <= start_range[1])]['NormalisedClipping'].sum()
            results_dict[f"{dup_id}_Pos_Range_Start"][sample] = start_sum

        # Sum the NormalisedClipping for the end range
        end_group = group[group['Type'] == 'End']
        if not end_group.empty:
            end_contig = end_group['Contig'].values[0]
            end_range = end_group[['Pos_Range_Start', 'Pos_Range_End']].values[0]
            end_sum = clipping_df[(clipping_df['Contig'] == end_contig) &
                                  (clipping_df['ClipPos'] >= end_range[0]) &
                                  (clipping_df['ClipPos'] <= end_range[1])]['NormalisedClipping'].sum()
            results_dict[f"{dup_id}_Pos_Range_End"][sample] = end_sum

# Convert the results dictionary to a DataFrame
results_df = pd.DataFrame.from_dict(results_dict, orient='index')
# Optional: if you want to convert the DataFrame such that Duplication_IDs are rows and Samples are columns
results_df = results_df.transpose()

# Save the results to a CSV file
results_df.to_csv('bijagos_clipping_summary.csv', index_label='Sample')


######### Testing impact of quality filtering the bam file for mapq >= 10
#1.
# AG0019-C_sorted.bam
# AG0135-C_sorted.bam
samtools view -b -q 10 AG0019-C_sorted.bam > AG0019-C_sorted_mapq10.bam
samtools view -b -q 10 AG0135-C_sorted.bam > AG0135-C_sorted_mapq10.bam

# 2. Create coverage file from sorted bam file using samtools depth
samtools depth -aa AG0019-C_sorted_mapq10.bam > AG0019-C_sorted_mapq10.coverage
samtools depth -aa AG0135-C_sorted_mapq10.bam > AG0135-C_sorted_mapq10.coverage
# 3. Create sam file from sorted bam,
samtools view -h AG0019-C_sorted_mapq10.bam -o AG0019-C_sorted_mapq10.sam
samtools view -h AG0135-C_sorted_mapq10.bam -o AG0135-C_sorted_mapq10.sam

## create normalised clipping files

while read -r sample; do
    # Construct the filenames
    sam_file="${sample}_sorted.mapq10.sam"
    coverage_file="${sample}_sorted.mapq10.coverage"

    # Run the ClipIdentifier.py script
    python /mnt/storage11/sophie/gitrepos/anophelesmelas_popgen/ClipIdentifier.py "$sam_file" "$coverage_file"
    echo "Processed $sample with ClipIdentifier.py"
done < mapq10_samples.txt

# should produce .mapq10_Normalised.csv for each sample
############# create mini-csv 

####### CREATE CSV OF SOFT-CLIPPING ##########

import pandas as pd

# Load the discordant read guide
discordant_df = pd.read_csv('discordant_read_guide.csv')

# Group the guide by Duplication_ID
grouped_discordant_df = discordant_df.groupby('Duplication_ID')

# Read the sample names from a file
with open('mapq10_samples.txt', 'r') as file:
    samples = file.read().splitlines()

# Initialize an empty dictionary to hold the results
results_dict = {}

# Iterate over the grouped DataFrame
for dup_id, group in grouped_discordant_df:
    results_dict[f"{dup_id}_Pos_Range_Start"] = {}
    results_dict[f"{dup_id}_Pos_Range_End"] = {}

    for sample in samples:
        # Initialize sum for each sample for start and end positions
        results_dict[f"{dup_id}_Pos_Range_Start"][sample] = 0
        results_dict[f"{dup_id}_Pos_Range_End"][sample] = 0

        # Load the clipping data for the sample
        clipping_df = pd.read_csv(f'{sample}_sorted.mapq10.Clipping.mapq10_Normalised.csv')

        # Sum the NormalisedClipping for the start range
        start_group = group[group['Type'] == 'Start']
        if not start_group.empty:
            start_contig = start_group['Contig'].values[0]
            start_range = start_group[['Pos_Range_Start', 'Pos_Range_End']].values[0]
            start_sum = clipping_df[(clipping_df['Contig'] == start_contig) &
                                    (clipping_df['ClipPos'] >= start_range[0]) &
                                    (clipping_df['ClipPos'] <= start_range[1])]['NormalisedClipping'].sum()
            results_dict[f"{dup_id}_Pos_Range_Start"][sample] = start_sum

        # Sum the NormalisedClipping for the end range
        end_group = group[group['Type'] == 'End']
        if not end_group.empty:
            end_contig = end_group['Contig'].values[0]
            end_range = end_group[['Pos_Range_Start', 'Pos_Range_End']].values[0]
            end_sum = clipping_df[(clipping_df['Contig'] == end_contig) &
                                  (clipping_df['ClipPos'] >= end_range[0]) &
                                  (clipping_df['ClipPos'] <= end_range[1])]['NormalisedClipping'].sum()
            results_dict[f"{dup_id}_Pos_Range_End"][sample] = end_sum

# Convert the results dictionary to a DataFrame
results_df = pd.DataFrame.from_dict(results_dict, orient='index')
# Optional: if you want to convert the DataFrame such that Duplication_IDs are rows and Samples are columns
results_df = results_df.transpose()

# Save the results to a CSV file
results_df.to_csv('test_clipping_summary.csv', index_label='Sample')

##
