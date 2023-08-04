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

# remove INDELS with bcftools, here -M2 indicates that bcftools should split multi-allelic sites into multiple biallelic sites, 
# keeping this information
# -m2 is used in conjunction with -M2 to apply the minor-allele-based decomposition. 
# This means that bcftools will decompose multi-allelic sites using the minor allele as the reference allele in the biallelic split.
# here -v tells bcftools to only view SNPS, so indels are excluded

bcftools view -M2 -m2 -v snps gambiae_nov2022.2023_07_05.genotyped.vcf.gz > bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf

# you must bgzip the file before indexing
bgzip bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf

# tabix index the compressed VCF file, creates .vcf.gz.tbi
tabix -p vcf bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz

## Make hdf5 file

allel.vcf_to_hdf5('bi_snps_chr_gambiae_nov2022.2023_07_05.genotyped.vcf.gz', 'bi_snps_chr_gambiae_nov2022_v2.h5', fields='*', overwrite=True)
http://alimanfoo.github.io/2017/06/14/read-vcf.html

## STEP 2: Conduct PCA and plot using the hdf5 file - IPYRAD METHOD

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
/eca-bioinf-handbook/bioinformatics-for-rad-seq-data-with-and-without-a-reference-genome.html
# show the first ten samples and the first 10 PC axes
df.iloc[:10, :10].round(2)

pca.draw(0, 2); #why is this not displaying in VSC?
plt.show() # ensure the plot remains shown

## Removing contigs and then making h5 file

# filter out the contigs from the VCF file, note that to produce a properly bgzipped vcf file you need the -Oz flag

bcftools view bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz --regions 2L,2R,3L,3R,Mt,X,Y_unplaced | bcftools sort -Oz -o bi_snps_chr_gambiae_nov2022.2023_07_05.genotyped.vcf.gz

# view the unique chromosome names
bcftools query -f '%CHROM\n' bi_snps_chr_gambiae_nov2022.2023_07_05.genotyped.vcf.gz | sort | uniq > unique_chromosomes_filtered.txt

# now the vcf file only contains chromosomes of interest, make hdf5 file 
# allel.vcf_to_hdf5('example.vcf', 'example.h5', fields='*', overwrite=True)
import h5py
allel.vcf_to_hdf5('bi_snps_chr_gambiae_nov2022.2023_07_05.genotyped.vcf.gz', 'bi_snps_chr_gamboae_nov2022_v2.h5', fields='*', overwrite=True)

## When converting a VCF file to an HDF5 file using allel.vcf_to_hdf5(), 
## the default behavior of h5py is to encode text strings as byte-strings (bytes) for efficient 
## storage and handling of different character encodings. 
## to use hf5 I need to be able to decode these strings
## Making PCA without converting to hf5, using VCF file instead.

import allel
import matplotlib.pyplot as plt
import numpy as np

# Load the VCF file
vcf_path = 'bi_snps_chr_gambiae_nov2022.2023_07_05.genotyped.vcf.gz'
gn = allel.read_vcf(vcf_path, fields=['calldata/GT'], numbers={'calldata/GT': 2})

# Reshape genotype data to 2D array
gn_2d = gn['calldata/GT'].reshape(gn['calldata/GT'].shape[0], -1)

# Find indices of NaN values
nan_indices = np.isnan(gn_2d)

# Remove rows with NaN values
gn_2d_cleaned = gn_2d[~np.any(nan_indices, axis=1)]

#Convert the data type to float64 for the allel package
gn_2d_cleaned = gn_2d_cleaned.astype(np.float64)

# Perform PCA
n_components = 2
pca_results = allel.pca(gn_2d_cleaned, n_components=n_components, copy=True, scaler='biallelic_patterson', ploidy=2)

# Extract PCA components
pca_components = pca_results['pca']

# Plot the first two principal components
plt.figure(figsize=(10, 6))
plt.scatter(pca_components[:, 0], pca_components[:, 1], alpha=0.5)
plt.title('PCA Plot')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.show()





# trying http://alimanfoo.github.io/2015/09/28/fast-pca.html





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




## creating conda enviroment that will work for h5 files

conda create -n scikit3 python=3.5
pip install h5py==2.10.0
mamba install scikit-allel==1.0.3 #need python 3.5 for this but apparently python needs upating as is depreciated
mamba install ipython


conda create -n scikit2 python=3.6
pip install h5py==2.10.0
mamba install scikit-allel
mamba install ipython

# or continue filtering as Holly did for:
# missingness
# normalise and then use minimum allele depth, bcftools norm 
# holly example: bcftools norm -m - AnS_0623.filter.chr.vcf.gz -Oz -o AnS_0623.namulti.vcf.gz -> remove multiallelic sites
# core genome. Holly did a bed file for the chromosomes, and also a bed file for the genes of interest for faster analysis
