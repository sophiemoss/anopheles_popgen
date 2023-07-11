# python script to make coverage plots across the genome for WGS data

## Coverage plots for WGS data:

## Before you use the code you need to use sambamba to create a cov.txt file from your bam file, so for my WGS sample it was:
sambamba depth window -w 10000 NG-31762_Bu1002_lib648117_10128_2.mkdup.bam -o NG-31762_Bu1002_lib648117_10128_2.cov.txt

## write a loop to do this

sambamba depth window -w 10000 Bu1002_Combined.mkdup.bam -o Bu1002_Combined.cov.txt

#note that the window for plasmodium was -w 1000 but the window for anopheles gambiae is -w10000 because the genome is 10x the size.

#check number of chromosomes in the cov.txt file
awk '{print $1}' NG-31762_Bu1002_lib648117_10128_2.cov.txt | sort | uniq

#create filtered.cov.txt file with only real chromsomes, not contigs

(fastq2matrix) sophie@s8:/mnt/storage8/sophie/bijagos_mosq_wgs/fastq2vcf_output$ python
Python 3.7.12 | packaged by conda-forge | (default, Oct 26 2021, 06:08:53) 
[GCC 9.4.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import pandas as pd
>>> df = pd.read_csv('NG-31762_Bu1002_lib648117_10128_2.cov.txt',sep='\t')
>>> df.head()
       # chrom  chromStart  chromEnd  readCount  meanCoverage                        sampleName
0  NC_004818.2           0     10000        597        7.1708  NG-30341_UN1011_lib602330_8122_4
1  NC_004818.2       10000     20000        590        6.6462  NG-30341_UN1011_lib602330_8122_4
2  NC_004818.2       20000     30000        525        6.3677  NG-30341_UN1011_lib602330_8122_4
3  NC_004818.2       30000     40000        506        5.8614  NG-30341_UN1011_lib602330_8122_4
4  NC_004818.2       40000     50000        663        7.3815  NG-30341_UN1011_lib602330_8122_4
>>> df.shape
(24882, 6)
>>> sub = df[df['# chrom'].isin(['NC_004818.2','NT_078266.2','NT_078268.4'])]
>>> sub.shape
(13913, 6)
>>> sub['# chrom'].unique()
array(['NC_004818.2', 'NT_078266.2', 'NT_078268.4'], dtype=object)
>>> sub.to_csv('NG-31762_Bu1002_lib648117_10128_2.cov.filtered.txt',sep='\t',index=None)
>>> 
(fastq2matrix) sophie@s8:/mnt/storage8/sophie/bijagos_mosq_wgs/fastq2vcf_output$ head NG-30341_UN1011_lib602330_8122_4.cov.filtered.txt

## Then run the python code with:
## you need to install matplotlib and pandas

conda install -c conda-forge matplotlib

conda install -c anaconda pandas

#python chromo_coverage.py  --sample 'NG-30341_UN1011_lib602330_8122_4' --csvfile 'NG-30341_UN1011_lib602330_8122_4.cov.txt'

python /mnt/storage8/sophie/bijagos_mosq_wgs/fastq2vcf_output_alignedtogambiae/coverageplots/chromo_coverage.py --sample 'NG-31762_Bu1002_lib648117_10128_2' --csvfile 'NG-31762_Bu1002_lib648117_10128_2.cov.filtered.txt'