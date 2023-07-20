##step 1, change vcf to hdf5 file for using scikit allel downstream

conda activate scikit

mamba install ipyrad -c bioconda
mamba install htslib -c bioconda
mamba install bcftools -c bioconda
mamba install vcftools -c bioconda

python
import ipyrad.analysis as ipa
import pandas as pd

# remove INDELS with bcftools

vcftools --gzvcf gambiae_nov2022.2023_07_05.genotyped.vcf.gz --remove-indels --recode --recode-INFO-all --out SNPS_only_gambiae_nov2022.2023_07_05.genotyped > vcftoolsfilter_log.txt 2>&1

#bcftools view -e 'TYPE="indel" | QUAL<=30' gambiae_nov2022.2023_07_05.genotyped.vcf.gz -Ou > gambiae_nov2022.2023_07_05.genotyped_SNPS_QUAL.vcf

# then zip
bgzip SNPS_only_gambiae_nov2022.2023_07_05.genotyped.recode.vcf

# tabix index the compressed VCF file, creates .vcf.gz.tbi
tabix -p vcf SNPS_only_gambiae_nov2022.2023_07_05.genotyped.recode.vcf.gz





#convert to hdf5 and then go through filtering steps with scikit allel workflow.
# instructions found here https://ipyrad.readthedocs.io/en/latest/API-analysis/cookbook-vcf2hdf5.html

import pandas as pd
import h5py
import ipyrad.analysis as ipa

# init a conversion tool
converter = ipa.vcf_to_hdf5(
    name="Gambiae2022_LD20K",
    data="gambiae_nov2022.2023_07_05.genotyped_SNPS_QUAL.vcf.gz",
    ld_block_size=20000,
)

# run the converter
converter.run()




# or continue filtering as Holly did for:
# missingness
# normalise and then use minimum allele depth, bcftools norm 
# holly example: bcftools norm -m - AnS_0623.filter.chr.vcf.gz -Oz -o AnS_0623.namulti.vcf.gz -> remove multiallelic sites
# core genome. Holly did a bed file for the chromosomes, and also a bed file for the genes of interest for faster analysis




