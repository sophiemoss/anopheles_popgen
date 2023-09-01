## scikitallel_workflow

conda activate scikit

######################## CALCULATING FST #########################

# %%
import zarr
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
import pandas as pd
import allel; print('scikit-allel', allel.__version__)

# %% Set wd
os.chdir('/mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf_agamP4/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022_combinedvcf_filteringsteps')
os.getcwd()

# %% CONVERT FILTERED UNPHASED VCF TO ZARR FILE
# Note this is using unphased VCF to convert to zarr at the moment
# allel.vcf_to_zarr('example.vcf', 'example.zarr', fields='*', overwrite=True)
# print('zarr', zarr.__version__, 'numcodecs', numcodecs.__version__)

callset = zarr.open('F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.zarr', mode='r')
callset.tree(expand=True)

# %%  CONVERT ZARR FILE TO GENOTYPE ARRAY
genotype_all = allel.GenotypeChunkedArray(callset['calldata/GT'])
genotype_all

# %% load variant positions
pos_all = callset['variants/POS'][:]

# %%  IMPORT METADATA
df_samples= pd.read_csv('metadata_gambiae_2022.csv',sep=',',usecols=['sample','year','country','island','phenotype'])
df_samples.head()
df_samples.groupby(by=['phenotype']).count()

# %%  choose sample populations to work with
pop1 = 'resistant'
pop2 = 'susceptible'
n_samples_pop1 = np.count_nonzero(df_samples.phenotype == pop1)
n_samples_pop2 = np.count_nonzero(df_samples.phenotype == pop2)
print(pop1, n_samples_pop1, pop2, n_samples_pop2)

# %% dictonary mapping population names to sample indices
subpops = {
    pop1: df_samples[df_samples.phenotype == pop1].index,
    pop2: df_samples[df_samples.phenotype == pop2].index,
}

# %% get allele counts
acs = genotype_all.count_alleles_subpops(subpops)
acs

# %% filter out variants that aren't segregating in the union of the two selected populations. 
# also filter out multiallelic variants (these should already be removed during filtering of vcf but just to check)

acu = allel.AlleleCountsArray(acs[pop1][:] + acs[pop2][:])
flt = acu.is_segregating() & (acu.max_allele() == 1)
print('retaining', np.count_nonzero(flt), 'SNPs')

# %% create the new genotype array

pos = pos_all.compress(flt)
ac1 = allel.AlleleCountsArray(acs[pop1].compress(flt, axis=0)[:, :2])
ac2 = allel.AlleleCountsArray(acs[pop2].compress(flt, axis=0)[:, :2])
genotype = genotype_all.compress(flt, axis=0)
genotype

############ additional methods comparisons ##########

##### Comparing Fst estimators ####
# %% ONE: calculate Weir and Cockerham Fst without using windows
# sample indices for population 1
pop1_idx = subpops[pop1]
# sample indices for population 2
pop2_idx = subpops[pop2]
a, b, c = allel.weir_cockerham_fst(genotype, subpops=[pop1_idx, pop2_idx], max_allele=1)
snp_fst_wc_no_windows = (a / (a + b + c))[:, 0]
snp_fst_wc_no_windows


# %% TWO: calculate Weir and Cockerham Fst using windows
# sample indices for population 1
pop1_idx = subpops[pop1]
# sample indices for population 2
pop2_idx = subpops[pop2]
a, b, c = allel.weir_cockerham_fst(genotype, subpops=[pop1_idx, pop2_idx], max_allele=1, blen=20000000)
snp_fst_wc_windowed = (a / (a + b + c))[:, 0]
snp_fst_wc_windowed


# %% THREE: calculate Hudson Fst, no windows
num, den = allel.hudson_fst(ac1, ac2)
snp_fst_hudson = num / den
snp_fst_hudson

# %% Compare Fst values between Weir & Cockerham and Hudson, with no windows

fig, ax = plt.subplots(figsize=(5, 5))
sns.despine(ax=ax, offset=5)
ax.plot(snp_fst_hudson, snp_fst_wc_no_windows, color='k', marker='.', linestyle=' ')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xlabel('Hudson $F_{ST}$')
ax.set_ylabel('Weir & Cockerham $F_{ST}$')
ax.set_title('%s (%s) vs %s (%s), SNP $F_{ST}$' % (pop1, n_samples_pop1, pop2, n_samples_pop2));

# there are unequal sample sizes so the Fst estimators do appear different, not great. 
# could use Hudson instead but W&C is standard.

# %%
### Compute chromosome-wide average Fst with standard errors approximated via a block-jackknife, to compare W&C and Hudson

# Weir & Cockerham
fst_wc, se_wc, vb_wc, _ = allel.blockwise_weir_cockerham_fst(genotype, subpops=[pop1_idx, pop2_idx],                                                                    blen=2000, max_allele=1)
print('%.04f +/- %.04f (Weir & Cockerham)' % (fst_wc, se_wc))

# Hudson
fst_hudson, se_hudson, vb_hudson, _ = allel.blockwise_hudson_fst(ac1, ac2, blen=100000)
print('%.04f +/- %.04f (Hudson)' % (fst_hudson, se_hudson))

# %% SNP ascertainment - run this test to see which SNPs to include in analysis
# Bhatia et al recommend ascertaining SNPs by choosing SNPs that are segregating in a third outgroup population
# If there is no obvious outgroup, there are four choices:
# 1) choose SNPs segregating in the first population
# 2) choose SNPs segregating in the second population
# 3) choose SNPs segregating in either population
# 4) choose SNPs segregating in both populations

def compute_fst(ac1, ac2, scheme):
    
    if scheme == 'first':
        loc_asc = ac1.is_segregating()
    elif scheme == 'second':
        loc_asc = ac2.is_segregating()
    elif scheme == 'either':
        loc_asc = ac1.is_segregating() | ac2.is_segregating()
    elif scheme == 'both':
        loc_asc = ac1.is_segregating() & ac2.is_segregating()    
    n_snps = np.count_nonzero(loc_asc)
    
    ac1 = ac1.compress(loc_asc, axis=0)
    ac2 = ac2.compress(loc_asc, axis=0)
    
    fst, se, _, _ = allel.blockwise_hudson_fst(ac1, ac2, blen=100000)
    
    print('%.04f +/- %.04f (using %s SNPs segregating in %s population)' % (fst, se, n_snps, scheme))

for scheme in 'first', 'second', 'either', 'both':
    compute_fst(ac1, ac2, scheme)

# use SNPs segregating in both populations (or could change this to use SNPs segregating in the control population)

# %% Calculate Fst using W&C, windowed, with SNPs that are segregating in both populations

import matplotlib.pyplot as plt

# Assuming you have already calculated snp_fst_wc_windowed and pos

# Create a figure and axis for the plot
fig, ax = plt.subplots(figsize=(10, 4))

# Plot Fst values along the genome
ax.plot(pos, snp_fst_wc_windowed, 'k-', lw=0.5)

# Customize the plot
ax.set_ylabel('$F_{ST}$')
ax.set_xlabel('Chromosome position (bp)')
ax.set_xlim(0, pos.max())  # Set the x-axis limits based on your data

# Optionally, you can add more customization, such as titles, legends, or gridlines
# For example:
# ax.set_title('Fst Along the Genome')
# ax.grid(True)

# %% Extract W&C Fst values above 0.6 and corresponding genomic positions 
### threshold > 0.6

# Weir and Cockerham Fst values
snp_fst_wc_filtered = snp_fst_wc_windowed[snp_fst_wc_windowed > 0.6]

# Genomic positions corresponding to filtered Fst values
genomic_positions_filtered = pos.compress(snp_fst_wc_windowed > 0.6)

# Initialize a list to store the window positions
window_positions = []

# Define the window size (blen)
blen = 2000  # Adjust this based on your window size

# Iterate through filtered Fst values and calculate the window positions
for position in genomic_positions_filtered:
    # Calculate the start position of the window containing 'position'
    window_start = (position // blen) * blen
    # Calculate the end position of the window
    window_end = window_start + blen
    # Add the start and end positions of the window to the list
    window_positions.append((window_start, window_end))

# Print the start and end positions of windows with Fst > 0.6
for start, end in window_positions:
    print(f"Start of Window with Fst > 0.6: {start}, End of Window: {end}")

# %% Extract W&C Fst values above 0.6 and corresponding genomic positions 
# Weir and Cockerham Fst values
snp_fst_wc_filtered = snp_fst_wc_windowed[snp_fst_wc_windowed > 0.6]

# Genomic positions corresponding to filtered Fst values
genomic_positions_filtered = pos.compress(snp_fst_wc_windowed > 0.6)

# Initialize lists to store window properties
window_positions = []
window_max_fst = []

# Define the window size (blen)
blen = 2000  # Adjust this based on your window size

# Iterate through filtered Fst values and calculate the window positions and maximum Fst values
for position in genomic_positions_filtered:
    # Calculate the start position of the window containing 'position'
    window_start = (position // blen) * blen
    # Calculate the end position of the window
    window_end = window_start + blen
    # Add the start and end positions of the window to the list
    window_positions.append((window_start, window_end))
    
    # Find indices within the current window
    indices_in_window = (pos >= window_start) & (pos < window_end)
    # Extract Fst values within the window
    fst_values_in_window = snp_fst_wc_windowed[indices_in_window]
    # Calculate the maximum Fst within the window
    max_fst_in_window = fst_values_in_window.max()
    window_max_fst.append(max_fst_in_window)

# Print the start, end positions, and maximum Fst values of windows with Fst > 0.6
for (start, end), max_fst in zip(window_positions, window_max_fst):
    print(f"Start of Window: {start}, End of Window: {end}, Maximum Fst in Window: {max_fst}")

# %%
# Find the maximum Fst value in snp_fst_wc_windowed
max_fst_value = snp_fst_wc_windowed.max()

# note that fst does not tell you the direction of selection, but high fst values indicate that
# there is genetic differentiation at that point between the populations being compared
# to understand the direction of selection, you can do additional analyses such as Tajimas D.


# %% Also used vcftools to calculate Fst snp by snp, too.
## after filtering only kept 33 out of 42 individuals

vcftools --gzvcf F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz \
--weir-fst-pop resistant_samples.txt \
--weir-fst-pop susceptible_samples.txt \
--out fst_resistant_vs_susceptible

CHROM   POS     WEIR_AND_COCKERHAM_FST
2L      1206    -0.0343433
2L      1402    0.0194649
2L      1444    -0.0284405
2L      1447    -0.0416782
2L      1449    -0.0245136
2L      1462    -0.0170085
2L      1463    -0.0170085
2L      1470    -0.0365088
2L      1487    -0.0265077

## this calculates for individual variants. Be aware that calculating Fst is often more common to use sliding windows to aggregate vriants over larger genomic regions.


####### Notes #######
# I want to analyse the component of variance between populations
# blen=None is not working, I need to specify a value for blen
# a is the component of variance between populations
# b is the component of variance between individuals within populations
# c is the component of variance between gametes within individuals

# in simpler terms

# a: How much genetic diversity is there within a single group of individuals from the same population?
# b: How different are the genetic traits between different groups of individuals (populations)?
# c: Do certain genetic traits tend to appear together in similar patterns across different populations?


### Plot in 100kb windows across each chromosome