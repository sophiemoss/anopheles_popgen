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