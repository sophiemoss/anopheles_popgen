# Using samtools depth instead to be able to keep mitochondria in, think these were being filtered out by sambamba when making chromo.cov csv fles

import subprocess
import pandas as pd

# Execute the first bash command to go through bam files and use sambamba to calculate depth of coverage in windows across the genome
# windows for plasmodium were 1000, but for anopheles I could use 10x this 10,000, because the genome is 10x the size
# samtools depth ur1001_Combined.mkdup.bam > ur1001_samtools_depth.csv

# subprocess.run('for f in *mkdup.bam ; do samtools depth "$f" > "$f.samtools_depth.csv" ; done', shell=True)

# take samtools_depth.csv and convert it to be in windows instead of each position

subprocess.run('''
for f in *mkdup.bam ; do
    echo -e "# chrom\tchromStart\tchromEnd\tmeanCoverage"

    samtools depth "$f" | \
    awk 'BEGIN {OFS="\\t"; sum=0; count=0; prevChrom=""; start=0; end=0}
         {
           if ($1 != prevChrom || NR == 1) {
             if (NR > 1) {
               printf "%s\\t%d\\t%d\\t%.2f\\n", prevChrom, start, end, sum/count;
             }
             start = $2; 
             end = $2 + 999; # changed this to 999 for 1000 base pair windows
             sum = $3; 
             count = 1;
             prevChrom = $1;
           } else {
             if ($2 <= end) {
               sum += $3; 
               count++;
             } else {
               printf "%s\\t%d\\t%d\\t%.2f\\n", prevChrom, start, end, sum/count;
               start = end + 1; 
               end = start + 999; # changed this to 999 for 1000 base pair windows 
               sum = $3; 
               count = 1;
             }
           }
         }
         END {printf "%s\\t%d\\t%d\\t%.2f\\n", prevChrom, start, end, sum/count}' >> "$f.samtools_windowed_depth.csv"
done
''' , shell = True)


# Filter out contigs from the windowed depth files
# Read the file containing sample names
with open('samples.txt', 'r') as file:
    sample_names = file.read().splitlines()

# Loop through each sample name to filter out contigs
for sample_name in sample_names:
    # Construct the file paths
    input_file = sample_name + '.mkdup.bam.samtools_windowed_depth.csv'
    output_file = sample_name + '.mkdup.bam.samtools_windowed_depth_filtered.csv'

    # Read the input file
    try:
       df = pd.read_csv(input_file, sep='\t')
    except pd.errors.EmptyDataError:
      print(f"Error: {input_file} is empty. Skipping...")
      continue

    # Filter the data
    sub = df[df['# chrom'].isin(['2L', '2R', '3L', '3R', 'X', 'Mt'])]

    # Export the filtered data to a new file
    sub.to_csv(output_file, sep='\t', index=None)

# Use the other python script to create figures of coverage over the genome
subprocess.run('cat samples.txt | parallel -j 30 "python /mnt/storage11/sophie/gitrepos/anophelesmelas_popgen/chromo_coverage.py --sample {} --csvfile {}.mkdup.bam.samtools_windowed_depth_filtered.csv" > chromocoverage.txt', shell=True)
