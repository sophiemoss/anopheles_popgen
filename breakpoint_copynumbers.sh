# Identifying possible CNVS from Lucas et al (2023) copy number breakpoints

# We will look for soft clipped reads in the bam files at the breakpoints indicated in the supplementary information in this paper
samtools view -f 0x2 -F 0x800 -b input.bam > output.bam

# This command uses two flags to filter the bam file
# -f 0x2: Keeps only the reads that are flagged as "properly paired". In other words, both reads of a pair (for paired-end sequencing data) are aligned in a manner that the aligner considers "proper" (e.g., they are aligned to the same chromosome, in the expected orientation, and within a certain distance of each other).
# -F 0x800: Filters out (removes) reads that have the supplementary alignment flag set. Supplementary alignments are typically part of chimeric reads and represent alignments that are not the primary alignment for a segment of the read.
# -b: Specifies that the output should be in BAM format.
# So the command only keeps reads that are properly paired and are not supplementary alignments, excluding chimeric segments or segments that didn't align with their mate in the expected manner.

# create list of bamfiles
ls *_1.mkdup.bam > bam_files.txt

# filter bam files as above

while read -r bamfile; do
    # Generate the output file name by appending "_filtered" before the ".bam" extension
    output="${bamfile%.*}_filtered.bam"
    
    # Apply the samtools command
    samtools view -f 0x2 -F 0x800 -b "$bamfile" > "$output"
    
    echo "Processed $bamfile -> $output"
done < bam_files.txt

## Use python script (ala Matt) to identify soft clipped reads