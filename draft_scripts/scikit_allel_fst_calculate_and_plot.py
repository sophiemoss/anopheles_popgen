## scikitallel_workflow

conda activate scikit

######################## CALCULATING FST #########################
## adapted jupyter notebook from http://alimanfoo.github.io/2015/09/21/estimating-fst.html

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
# callset.tree(expand=True)

# %% choose how much of genome to look at
# load variant positions for whole genome or section of genome
# pos_all = callset['variants/POS'][:]
# chrom_all = callset['variants/CHROM'][:]
# load variant positions for specific chromosome

chromosome_filter = callset['variants/CHROM'][:] == '3L'
pos_all = callset['variants/POS'][np.where(chromosome_filter)[0]] 

# %%  CONVERT ZARR FILE TO GENOTYPE ARRAY
# whole genome
# genotype_all = allel.GenotypeChunkedArray(callset['calldata/GT'])
# genotype_all

# create genotype array for just chrom
genotype_all = allel.GenotypeChunkedArray(callset['calldata/GT'][np.where(chromosome_filter)[0]])

# check length of pos_all and genotype_all are the same
print(len(pos_all),len(genotype_all))

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

# check that pos and genotype are the same size
print(len(pos),len(genotype))

# %% PLOT FST using windowed weird and cockerham fst

subpoplist = [list(subpops['resistant']),
                  list(subpops['susceptible'])]

fst, windows, counts = allel.windowed_weir_cockerham_fst(pos, genotype, subpoplist, size=100000)

# use the per-block average Fst as the Y coordinate
y = fst

x = [np.mean(w) for w in windows]

# plot
fig, ax = plt.subplots(figsize=(10, 4))
sns.despine(ax=ax, offset=5)
ax.plot(x, y, 'k-', lw=.5)
ax.set_ylabel('$F_{ST}$')
ax.set_xlabel('Chromosome 3L position (bp)')
ax.set_xlim(0, pos.max())

print(len(x),len(y))

# INSPECT Fst values
# %% Print the maximum FST value and its corresponding window
max_value = max(fst)
max_index = np.argmax(fst)
max_window = windows[max_index]

print("Maximum FST Value:", max_value)
print("Corresponding Window:", max_window)

# %% Calculte fst with blen = 1 and plot
# PLOT FST using windowed weird and cockerham fst - try this, not sure if x = pos

subpoplist = [list(subpops['resistant']),
                  list(subpops['susceptible'])]


# fst, se, vb, _ = allel.blockwise_hudson_fst(ac1, ac2, blen=blen)
a, b, c = allel.weir_cockerham_fst(genotype, subpoplist, blen=1)

# estimate Fst for each variant and each allele individually
# fst = a / (a + b + c)

# or estimate Fst for each variant individually, averaging over alleles

fst = (np.sum(a, axis=1) /
       (np.sum(a, axis=1) + np.sum(b, axis=1) + np.sum(c, axis=1)))

# use the per-block average Fst as the Y coordinate
y = fst

x = pos

# plot
fig, ax = plt.subplots(figsize=(10, 4))
sns.despine(ax=ax, offset=5)
ax.plot(x, y, 'k-', lw=.5)
ax.set_ylabel('$F_{ST}$')
ax.set_xlabel('Chromosome 3L position (bp)')
ax.set_xlim(0, pos.max())

print(len(x),len(y))



