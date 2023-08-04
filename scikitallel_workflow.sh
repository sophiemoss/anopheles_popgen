## Step 1: Filtering vcf to remove indels and split multiallelic sites into multiple biallelic sites, and convert to hdf5 file
## Use http://alimanfoo.github.io/2017/06/14/read-vcf.html instead of ipyrad to make hdf5 file
## https://ipyrad.readthedocs.io/en/stable/API-analysis/cookbook-pca.html

conda activate scikit

mamba install ipyrad -c bioconda
mamba install htslib -c bioconda
mamba install bcftools -c bioconda
mamba install vcftools -c bioconda

## 

python

import pandas as pd
import h5py
import ipyrad.analysis as ipa
import matplotlib.pyplot as plt

## Decided not to make hdf5 files from vcf
#allel.vcf_to_hdf5('bi_snps_chr_gambiae_nov2022.2023_07_05.genotyped.vcf.gz', 'bi_snps_chr_gambiae_nov2022_v2.h5', fields='*', overwrite=True)
#http://alimanfoo.github.io/2017/06/14/read-vcf.html
# working with binary matrices instead, with scikit-allel

# sort samples using imap

imap = {
    "sus": ["NG-33833_sus10_lib708284_10265_1", "NG-33833_sus18_lib708285_10265_1", "NG-33833_sus1_lib708279_10265_1", "NG-33833_sus25_lib708286_10265_1",
         "NG-33833_sus2_lib708280_10265_1", "NG-33833_sus30_lib708290_10265_1", "NG-33833_sus33_lib708287_10265_1", "NG-33833_sus35_lib708289_10265_1",
         "NG-33833_sus3_lib708281_10265_1", "NG-33833_sus40_lib708288_10265_1", "NG-33833_sus41_lib708291_10265_1", "NG-33833_sus4_lib708282_10265_1",
         "NG-33833_sus7_lib708283_10265_1]"],
    "res": ["NG-33833_res10_lib708250_10265_1", "NG-33833_res11_lib708251_10265_1", "NG-33833_res12_lib708252_10265_1", "NG-33833_res15_lib708253_10265_1", "NG-33833_res16_lib708254_10265_1", "NG-33833_res17_lib708255_10265_1", "NG-33833_res24_lib708256_10265_1", "NG-33833_res29_lib708257_10265_1",
    "NG-33833_res35_lib708258_10265_1", "NG-33833_res36_lib708259_10265_1", "NG-33833_res3_lib708245_10265_1", "NG-33833_res43_lib708260_10265_1", "NG-33833_res45_lib708261_10265_1", "NG-33833_res46_lib708262_10265_1", "NG-33833_res49_lib708263_10265_1", "NG-33833_res50_lib708264_10265_1",
    "NG-33833_res52_lib708265_10265_1", "NG-33833_res54_lib708266_10265_1", "NG-33833_res55_lib708267_10265_1", "NG-33833_res56_lib708268_10265_1", "NG-33833_res57_lib708269_10265_1", "NG-33833_res5_lib708246_10265_1", "NG-33833_res6_lib708247_10265_1", "NG-33833_res6_lib708247_10265_1",
    "NG-33833_res8_lib708248_10265_1", "NG-33833_res9_lib708249_10265_1" ],
    "con": ["NG-33833_c12_lib708275_10265_1", "NG-33833_c14_lib708276_10265_1", "NG-33833_c16_lib708277_10265_1", "NG-33833_c17_lib708278_10265_1", "NG-33833_c2_lib708270_10265_1", "NG-33833_c4_lib708271_10265_1", "NG-33833_c5_lib708272_10265_1", "NG-33833_c7_lib708273_10265_1", "NG-33833_c9_lib708274_10265_1" ],
    "5xres": ["NG-33833_5xres1_lib708243_10265_1", "NG-33833_5xres2_lib708244_10265_1"]
}


# view the unique chromosome names
bcftools query -f '%CHROM\n' bi_snps_chr_gambiae_nov2022.2023_07_05.genotyped.vcf.gz | sort | uniq > unique_chromosomes_filtered.txt

# or continue filtering as Holly did for:
# missingness
# normalise and then use minimum allele depth, bcftools norm 
# holly example: bcftools norm -m - AnS_0623.filter.chr.vcf.gz -Oz -o AnS_0623.namulti.vcf.gz -> remove multiallelic sites
# core genome. Holly did a bed file for the chromosomes, and also a bed file for the genes of interest for faster analysis
