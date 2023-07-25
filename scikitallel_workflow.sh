## Step 1: IPYRAD filtering vcf to remove indels and split multiallelic sites into multiple biallelic sites, and convert to hdf5 file
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

# remove INDELS with bcftools

vcftools --gzvcf gambiae_nov2022.2023_07_05.genotyped.vcf.gz --remove-indels --recode --recode-INFO-all --out SNPS_only_gambiae_nov2022.2023_07_05.genotyped > vcftoolsfilter_log.txt 2>&1

# or use bcftools, here -M2 indicates that bcftools should split multi-allelic sites into multiple biallelic sites, keeping this information
# -m2 is used in conjunction with -M2 to apply the minor-allele-based decomposition. 
# This means that bcftools will decompose multi-allelic sites using the minor allele as the reference allele in the biallelic split.
# here -v tells bcftools to only view SNPS, so indels are excluded

bcftools view -M2 -m2 -v snps gambiae_nov2022.2023_07_05.genotyped.vcf.gz > bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf

# you must bgzip the file before indexing
bgzip bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf

# tabix index the compressed VCF file, creates .vcf.gz.tbi
tabix -p vcf bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz

# load the VCF as an dataframe
dfchunks = pd.read_csv(
    "bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz",
    sep="\t",
    skiprows=1000,
    chunksize=1000,
)

# show first few rows of first dataframe chunk
next(dfchunks).head()

# init a conversion tool
converter = ipa.vcf_to_hdf5(
    name="Gambiae2022_LD20K",
    data="bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz",
    ld_block_size=20000,
)

# run the converter
converter.run()


## STEP 2: Conduct PCA and plot using the hdf5 file

import pandas as pd
import h5py
import ipyrad.analysis as ipa
import matplotlib.pyplot as plt

# look at the structure of the hdf5 file using nexusformat
# pip install nexusformat

import nexusformat.nexus as nx
f = nx.nxload(‘myhdf5file.hdf5’)
print(f.tree)

# init a PCA tool and filter to allow no missing data
pca = ipa.pca(
    data="./bi_Gambiae2022_LD20K.snps.hdf5",
    mincov=1.0,
    )

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

# require that 50% of samples have data in each group
minmap = {i: 0.5 for i in imap}

pca.run()

# store the PC axes as a dataframe
df = pd.DataFrame(pca.pcaxes[0], index=pca.names)

# write the PC axes to a CSV file
df.to_csv("pca_analysis.csv")

# show the first ten samples and the first 10 PC axes
df.iloc[:10, :10].round(2)

pca.draw(0, 2); #why is this not displaying in VSC?
plt.show() # ensure the plot remains shown

#######################
# trying http://alimanfoo.github.io/2015/09/28/fast-pca.html


import h5py

# Assuming callset is an HDF5 file
with h5py.File('./bi_Gambiae2022_LD20K.snps.hdf5', 'r') as f:
    # Print the keys at the root level
    print("Keys at root level:", list(f.keys()))

    # If there is a group called '2L' in the HDF5 file, you can check its keys as well
    if '2L' in f:
        print("Keys under '2L':", list(f['2L'].keys()))
    else:
        print ("No 2L")

# problem might be because the data is not phased before it is being converted into an hdf5 file.



import gcsfs


## Malariagen snp-genotyping-vector.md pipeline

# The GCS URI to the file you want to download
gcs_uri = "gs://vo_agam_production/resources/observatory/ag.allsites.nonN.vcf.gz"

# Create a GCS filesystem object
gcs_filesystem = gcsfs.GCSFileSystem()

# Download the file to a destination
destination_file_name = "/mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022_combinedvcf_filteringsteps/"
with gcs_filesystem.open(gcs_uri, "rb") as gcs_file, open(destination_file_name, "wb") as local_file:
    local_file.write(gcs_file.read())

print("File downloaded successfully.")

## permission was denied to download this file from google cloud storage.






## some of the analyses want a zarr file


# or continue filtering as Holly did for:
# missingness
# normalise and then use minimum allele depth, bcftools norm 
# holly example: bcftools norm -m - AnS_0623.filter.chr.vcf.gz -Oz -o AnS_0623.namulti.vcf.gz -> remove multiallelic sites
# core genome. Holly did a bed file for the chromosomes, and also a bed file for the genes of interest for faster analysis
