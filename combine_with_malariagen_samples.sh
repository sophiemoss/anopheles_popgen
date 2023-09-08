########## MALARIAGEN DATA - phase 3 ##################

## download vcfs from malariagen for Ag1000G
wget ftp:'//ngs.sanger.ac.uk/production/ag1000g/phase2/AR1/variation/main/vcf/all/*.vcf.gz'

## concatenate them together as they are all for different chromosomes

ls *vcf.gz > to_concat.txt
bcftools concat -f to_concat.txt -O z -o concatenated_unfiltered_ag1000g_phase2.vcf.gz

## bcftools concatenate is used for when you want to merge 2 or more VCFs with the exact same samples, 
## but for example, each VCF is for a single chromosome. 
## bcftools merge is when you want to merge 2 or more VCFs with non-overlapping samples 

## UP TO HERE.
## merge my samples with malariagen samples vcf when they have been aligned to Agam_P3 and have been combined into one large VCF

bcftools merge input1.vcf.gz concatenated_unfiltered_ag1000g_phase2.vcf.gz -O z -o malariagen_phase2_bijagos_unfiltered.vcf.gz