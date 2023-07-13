##step 1, change vcf to hdf5 file for using scikit allel downstream

conda activate scikit

mamba install ipyrad -c bioconda
mamba install htslib -c bioconda
mamba install bcftools -c bioconda
mamba install vcftools -c bioconda

python
import ipyrad.analysis as ipa
import pandas as pd

# tabix index the compressed VCF file, creates .vcf.gz.tbi
tabix gambiae_nov2022.2023_07_05.genotyped.vcf.gz

# remove INDELS

vcftools --gzvcf gambiae_nov2022.2023_07_05.genotyped.vcf.gz --remove-indels --recode --recode-INFO-all --out SNPS_only_gambiae_nov2022.2023_07_05.genotyped

#convert to hdf5 and then go through filtering steps with scikit allel workflow.

# or continue filtering for:
# missingness
# normalise and then use minimum allele depth, bcftools norm 
# holly example: bcftools norm -m - AnS_0623.filter.chr.vcf.gz -Oz -o AnS_0623.namulti.vcf.gz -> remove multiallelic sites
# core genome. Holly did a bed file for the chromosomes, and also a bed file for the genes of interest for faster analysis




