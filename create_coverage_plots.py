# use the fastq2vcfsamples.txt file which was created using:
# ls *_1.fastq.gz | sed 's/_1.fastq.gz//' > fastq2vcfsamples.txt

# windows for plasmodium were 1000, but for anopheles I have used 10000 as the genome is 10x the size

# OLD LINE cat fastq2vcfsamples.txt | parallel -j 30 --bar 'sambamba depth window -w 10000 {}.mkdup.bam -o {}.cov.txt' > sambamba_cov.txt_log.txt 2>&1

for f in *mkdup.bam ; do sambamba depth window -w 50000 "$f" -o "$f.chromocov" ; done #this is currently running in screen 11/7/23

cat fastq2vcfsamples.txt | parallel -j 30 "mv {}.mkdup.bam.chromocov {}_chromocov.csv

# python script to make coverage plots across the genome for WGS data

import matplotlib
import pandas as pd

# remove contigs from the chromocov files

df = pd.read_csv('bu1002_Combined.cov.txt',sep='\t')
df.head()
df.shape
print(df['# chrom'].unique())
sub = df[df['# chrom'].isin(['2L','2R','3L','3R','X','Mt'])]
sub.shape
sub['# chrom'].unique()
# export this file to a csv
sub.to_csv('bu1002_Combined.cov.filtered.txt',sep='\t',index=None)

#python chromo_coverage.py  --sample 'NG-30341_UN1011_lib602330_8122_4' --csvfile 'NG-30341_UN1011_lib602330_8122_4.cov.txt'

cat fastq2vcfsamples.txt | parallel -j 30 "python /mnt/storage11/sophie/gitrepos/anophelesmelas_popgen/chromo_coverage.py --sample {} --csvfile {}_chromocov.csv"
