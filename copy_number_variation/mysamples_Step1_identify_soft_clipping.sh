
# Identifying possible CNVS from Lucas et al (2023) by identifying soft clipping associated with different CNVs.
# This script prepares bam files by sorting them and filtering them for mapping quality >= 10
# It then uses a different script (ClipIdentifier.py, which was made by Matt Higgins) to find soft clipped reads in these bam files 
# ClipIdentifier.py normalises the % soft clipped reads against coverage of those reads


# set working directory
#cd /mnt/storage11/sophie/bijagos_mosq_wgs/malariagen_wgs/malariagen_GMABC-GW-bams
#current_dir=$(pwd)
#echo "The current working directory is: $current_dir"

# create list of bam files to work through
#ls *.bam > bam_files.txt

# EXAMPLE COMMANDS FOR SINGLE SAMPLE 
# 1. Sort all bam files using samtools
#samtools sort AJ0126-C.bam > AJ0126-C.sorted.bam

# 2. Filter bam files for mapping quality >= 10
#samtools view -b -q 10 AJ0126-C.sorted.bam > AJ0126-C.sorted_mapq10.bam

# 3. Create coverage file using samtools depth
#samtools depth -aa AJ0126-C.sorted_mapq10.bam > AJ0126-C.sorted_mapq10.coverage

# 4. Create sam file from sorted bam,
#samtools view -h AJ0126-C.sorted_mapq10.bam -o AJ0126-C.sorted_mapq10.sam

# Loop above steps for multiple samples:
# Sort bam files to make sorted bam file for each sample
# Create coverage file from sorted bam file using samtools depth
# Create sam file from sorted bam using samtools view

while read -r bamfile; do
    # Apply samtools sort command
    # Generate the output file name by appending ".sorted" before the ".bam" extension
    #output="${bamfile%.*}.sorted.bam"
    #samtools sort "$bamfile" > "$output"
    #echo "Processed $bamfile -> $output by sorting bamfile"

    # Apply samtools filter command
    filteredsortedbam="${bamfile%.*}.filtered.bam"
    samtools view -b -q 10 "$bamfile" > "$filteredsortedbam"
    echo "Processed $bamfile ->$filteredsortedbam by filtering the sorted bamfile"

    # Make the coverage file from the sorted and filtered bam file
    coverage_file="${filteredsortedbam%.*}.coverage"
    samtools depth -aa "$filteredsortedbam" > "$coverage_file"
    echo "Coverage generated for $filteredsortedbam -> $coverage_file"

    # Apply the third samtools command
    samfile="${filteredsortedbam%.*}.sam"
    samtools view -h "$filteredsortedbam" > "$samfile"

    echo "Sam file generated for $filteredsortedbam > $samfile"

done < bam_files.txt

# Run the clip_identifier python script using loop

ls *.sam | sed 's/.sam//' > samples.txt

while read -r sample; do
    # Run the ClipIdentifier.py script
    python /mnt/storage11/sophie/gitrepos/anophelesmelas_popgen/copy_number_variation/ClipIdentifier.py "$sample.sam" "$sample_sorted.coverage"   
    echo "Processed $sample with ClipIdentifier.py"
done < samples.txt