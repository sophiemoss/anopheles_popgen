# anopheles_popgen
Anopheles gambiae s.s. population genomic analysis

1. This project has involved sequenced extracted DNA (post WGA)
2. Sequencing has happened multiple times for each sample to get achieve more reads
3. Sequence fastq files for each sample have been combined using SophiesFastqMerge.py 
4. Combined sequences were then checked using **check_fastq.sh** for fastq file format
5. Sequences were then run 21st June - 24th June using 60 threads on s11 with **runfastq2vcf.sh**, which uses the **fastq2vcf.py** malaria-hub script. There are 45 melas samples in the population
7. Basic statistics were then calculated for each of the run sample, see **basic_statistics_analysis_wgs.sh**
8. A genomics database was made which combines the individual sample VCFs. This was done using **merge_vcfs.py** (**see makegenomicsdb.sh**)
9. This was genotyped using **merge_vcfs.py genotype**, to create one VCF file which is genotyped and contains all variants for all samples across all positions 

_The resulting database vcf is called gambiae_nov2022.2023_07_05.genotyped.vcf.gz_

10. This resulting vcf was then filtered using steps in **filter_combinegenotyped_vcf.sh**
