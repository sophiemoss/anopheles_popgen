import subprocess
import pandas as pd

# use the fastq2vcfsamples.txt file which was created using:
# ls *_1.fastq.gz | sed 's/_1.fastq.gz//' > fastq2vcfsamples.txt

# Execute the first bash command to go through bam files and use sambamba to calculate depth of coverage in windows across the genome
# windows for plasmodium were 1000, but for anopheles I have used 10000 as the genome is 10x the size
subprocess.run('for f in *mkdup.bam ; do sambamba depth window -w 50000 "$f" -o "$f.chromocov" ; done', shell=True)
subprocess.run('cat fastq2vcfsamples.txt | parallel -j 1 "mv {}.mkdup.bam.chromocov {}_chromocov.csv"', shell=True)

# Read the file containing sample names
with open('fastq2vcfsamples.txt', 'r') as file:
    sample_names = file.read().splitlines()

# Loop through each sample name to filter out contigs
for sample_name in sample_names:
    # Construct the file paths
    input_file = sample_name + '_chromocov.csv'
    output_file = sample_name + '_chromocov.filtered.csv'

    # Read the input file
    df = pd.read_csv(input_file, sep='\t')

    # Filter the data
    sub = df[df['# chrom'].isin(['2L', '2R', '3L', '3R', 'X', 'Mt'])]

    # Export the filtered data to a new file
    sub.to_csv(output_file, sep='\t', index=None)

# Use the other python script to create figures of coverage over the genome
subprocess.run('cat fastq2vcfsamples.txt | parallel -j 30 "python /mnt/storage11/sophie/gitrepos/anophelesmelas_popgen/chromo_coverage.py --sample {} --csvfile {}_chromocov.filtered.csv" > chromocoverage.txt', shell=True)
    
# single line outside python
# cat fastq2vcfsamples.txt | parallel -j 30 "python chromo_coverage.py --sample {} --csvfile {}_chromocov.filtered.csv" > chromocoveragelog.txt