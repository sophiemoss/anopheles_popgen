## scikitallel_workflow
# /mnt/storage11/sophie/gitrepos/anophelesmelas_popgen/scikit_allel_fst_calculate_and_plot_v2.py
# run this script using
# python scikit_allel_fst_calculate_and_plot.py /path/to/working/directory /path/to/callset.zarr chromosomename

######################## CALCULATING FST #########################
## adapted jupyter notebook from http://alimanfoo.github.io/2015/09/21/estimating-fst.html

import os
import argparse
import zarr
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
import pandas as pd
import allel  # Import allel at the top of the script

def main(args):
    print('scikit-allel', allel.__version__)
    
    # set working directory
    os.chdir(args.working_directory)
    print("Working directory:", os.getcwd())

    # Open the callset file
    callset = zarr.open(args.callset_file, mode='r')

    # Filter by chromosome
    chromosome_filter = callset['variants/CHROM'][:] == args.chromosome
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

    # %% Calculte fst with blen = 1 so that fst is calculated for every variant. Note this takes a bit longer to run.
    # PLOT FST using windowed weird and cockerham fst - try this, not sure if x = pos

    subpoplist = [list(subpops['resistant']),
                    list(subpops['susceptible'])]


    # fst, se, vb, _ = allel.blockwise_hudson_fst(ac1, ac2, blen=blen)
    a, b, c = allel.weir_cockerham_fst(genotype, subpoplist, blen=1)

    # estimate Fst for each variant and each allele individually
    # fst = a / (a + b + c)

    # or estimate Fst for each variant individually, averaging over alleles

    fst_pervariant = (np.sum(a, axis=1) /
        (np.sum(a, axis=1) + np.sum(b, axis=1) + np.sum(c, axis=1)))

    # use the per-block average Fst as the Y coordinate
    y = fst_pervariant

    x = pos

    # plot
    fig, ax = plt.subplots(figsize=(10, 4))
    sns.despine(ax=ax, offset=5)
    ax.plot(x, y, 'k-', lw=.5)
    ax.set_ylabel('$F_{ST}$')
    ax.set_xlabel('Chromosome 3L position (bp)')
    ax.set_xlim(0, pos.max())

    print(len(x),len(y))

    # Here you can see that fst values are much higher when calculated for each position individually

    print("Maximum FST Value:", max_value)
    print("Corresponding Window:", max_window)

# arguments

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate and plot FST using scikit-allel.")
    parser.add_argument("working_directory", type=str, help="Working directory for the script.")
    parser.add_argument("callset_file", type=str, help="Path to the callset file.")
    parser.add_argument("chromosome", type=str, help="Chromosome to analyze.")
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)