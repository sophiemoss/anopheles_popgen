#If bamstats files are not produced automatically use the following line to produce them

#cat samples.txt | parallel -j 2 "samtools flagstat {}.bqsr.bam > {}.bqsr.bamstats"


#Eurofins needs to send over 3.75 million reads (total read pairs)
#this outputs number of lines in file, divide this by 4 because fastq file has 4 lines per paired read = gives number of paired reads

#zcat fastqfilename | wc - 1

# do this in parallel for each sample, to calculate reads from fastq files

#cat samples.txt | parallel --keep-order -j 2 "zcat {}_2.fastq.gz | wc -l" > paired_reads

#then divide these by 4 to find the total number of paired reads


## After running the fastq2matrix pipeline we can extract information from these files to gain an idea of their quality.
## We are interested in:
## Number of reads
## Percentage of reads mapping to the reference genome
## Coverage
## Number of SNPs

# To extract the number of reads
# For multiple samples, feed in a list of all the samples and use xargs

#code below is saying show the first line (head -n 1) and extract the first item (print $1)
#it creates a list (total_reads) of the number for all samples in the order of the sample list

ls *.mkdup.bamstats | xargs -i -P1 sh -c 'head -1 {} | awk -F '"'"' '"'"' '"'"'{print $1}'"'"'' > total_reads
#ls *.bqsr.bamstats | xargs -i -P1 sh -c 'head -1 {} | awk -F '"'"' '"'"' '"'"'{print $1}'"'"'' > total_reads

# Number of reads mapping to reference genome
# Same principle, but extracting the first item on the 5th line

ls *.mkdup.bamstats| xargs -i -P1 sh -c 'head -5 {} | tail -1 | awk -F '"'"' '"'"' '"'"'{print $1}'"'"'' > mapped_reads

# Look at coverage as a %
# Step 1 - Create a list of all the sample names, you just want the prefix

ls *.mkdup.bam | sed -r 's/.mkdup.bam//g' > covlist

#ls *.bqsr.bam | sed -r 's/.bqsr.bam//g' > covlist

# Step 2 - Create a text file with the bash script that you want to run in parallel

#vim coverage_script.sh
#i
#samtools depth $1.mkdup.bam -a > $1.coverage 

#samtools depth $1.bqsr.bam -a > $1.coverage 

# Step 3 - Call the command using parallel - to create coverage files
# Open the sample list, then using that list use parallel to go through the list,
# one by one (-j 1) and execute the bash script coverage_script.sh. The --bar is to show a progress bar.

 cat covlist | parallel -j 25 --bar bash coverage_script.sh {}

 # Step 4 - Extract info from coverage files
 # Now we are interested how many postions in the genome have a coverage of >5 or >10
# Use following script to go through each cov file and extract the number of positions that have a coverage of 5 or above - easily modified to 10
# For mosquitoes could use 10 but Holly is using 5, which is what Susana suggested
# This creates a list in a text file which can easily be copied into a table - divide by genome length for proportion

ls *.coverage | parallel --keep-order -j 60 --bar 'cat {} | awk '"'"'{if($3>=5)print}'"'"' | wc -l' > coverage_5
ls *.coverage | parallel --keep-order -j 1 --bar 'cat {} | awk '"'"'{if($3>=10)print}'"'"' | wc -l' > coverage_10
ls *.coverage | parallel --keep-order -j 1 --bar 'cat {} | awk '"'"'{if($3>=20)print}'"'"' | wc -l' > coverage_20

## From Matt
## Code to look at what the coverage is of wgs data by looking at the bam file and create coverage plots across genome

## Check coverage as an average read depth across all genome

samtools depth -aa bamfile.mkdup.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'

example:
#samtools depth -aa bu1002_Combined.mkdup.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'
#for multiple samples:

ls *Combined.mkdup.bam | parallel --keep-order -j 1 --bar 'samtools depth -aa {} | awk '\''{sum+=$3} END { print sum/NR}'\''' > average_read_depth


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





