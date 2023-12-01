## scikitallel_workflow
# /mnt/storage11/sophie/gitrepos/anophelesmelas_popgen/scikit_allel_fst_calculate_and_plot_v2.py
# run this script using
# python scikit_allel_fst_calculate_and_plot.py /path/to/working/directory /path/to/callset.zarr chromosomename
# You make the zarr file with allel.vcf_to_zarr('phased_vcf_file_name.vcf.gz', 'output_name.zarr', fields='*', overwrite=True)
# allel.vcf_to_zarr('2019_melas_phased.vcf.gz', '2019_melas_phased.zarr', fields='*', overwrite=True)
######################## CALCULATING FST #########################
# %% adapted jupyter notebook from http://alimanfoo.github.io/2015/09/21/estimating-fst.html

import os
import argparse
import zarr
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
import pandas as pd
import allel
import sys
from datetime import datetime

# %%
def main(args):
    print('scikit-allel', allel.__version__)
    
    # set working directory
    os.chdir(args.working_directory)
    print("Working directory:", os.getcwd())

    # Open the callset file
    callset = zarr.open(args.callset_file, mode='r')
    print("Callset file:", args.callset_file)

    # Filter by chromosome
    chromosome_filter = callset['variants/CHROM'][:] == args.chromosome
    pos_all = callset['variants/POS'][np.where(chromosome_filter)[0]]
    print("Chromosome being analysed:", args.chromosome)
    print("Number of variants in chrom:", len(pos_all))

    # %%  CONVERT ZARR FILE TO GENOTYPE ARRAY
    # whole genome
    # genotype_all = allel.GenotypeChunkedArray(callset['calldata/GT'])
    # genotype_all

    # create genotype array for just chrom
    genotype_all = allel.GenotypeChunkedArray(callset['calldata/GT'][np.where(chromosome_filter)[0]])

    # check length of pos_all and genotype_all are the same

    print("Length of pos_all variable:",len(pos_all))
    print("Length of genotype_all variable:",len(genotype_all))
    
    if len(pos_all)==len(genotype_all):
        print("Length of positions and genotypes in the genotype array are the same, script continuing")
    else:
        print("Something is wrong with the genotype_all array as the length of pos_all and genotype_all are different. Stopping script.")
        sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line

    # %%  IMPORT METADATA
    df_samples= pd.read_csv('metadata_gambiae_2022.csv',sep=',',usecols=['sample','year','country','island','phenotype'])
    df_samples.head()
    df_samples.groupby(by=['phenotype']).count()
    print("Imported metadata")

    # %%  choose sample populations to work with
    pop1 = 'resistant'
    pop2 = 'susceptible'
    n_samples_pop1 = np.count_nonzero(df_samples.phenotype == pop1)
    n_samples_pop2 = np.count_nonzero(df_samples.phenotype == pop2)
    print("Population 1:", pop1, "Number of samples in pop1:", n_samples_pop1, "Population 2:", pop2, "Number of samples in pop2:", n_samples_pop2)

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
    print('Filtered out variants that are not segretating, or are not biallelic. Now retaining', np.count_nonzero(flt), 'SNPs')

    # %% create the new genotype array

    pos = pos_all.compress(flt)
    ac1 = allel.AlleleCountsArray(acs[pop1].compress(flt, axis=0)[:, :2])
    ac2 = allel.AlleleCountsArray(acs[pop2].compress(flt, axis=0)[:, :2])
    genotype = genotype_all.compress(flt, axis=0)
    genotype

    print("Created genotype array")

    # check that pos and genotype are the same size
    if len(pos)==len(genotype):
        print("Length of positions and genotypes in the genotype array are the same, script continuing")
    else:
        print("Something is wrong with the genotype_all array as the lenght of pos_all and genotype_all are different. Stopping script.")
        sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line

    # %% PLOT FST using windowed patterson fst
    # allel.windowed_patterson_fst(pos, ac1, ac2, size=None, start=None, stop=None, step=None, windows=None, fill=nan)
    # allel.windowed_weir_cockerham_fst(pos, g, subpops, size=None, start=None, stop=None, step=None, windows=None, fill=nan, max_allele=None)

    print("Plotting Fst using allel.windowed_patterson_fst, window size 1000")

    #subpoplist = [list(subpops['resistant']),
    #                list(subpops['susceptible'])]

    #fst, windows, counts = allel.windowed_weir_cockerham_fst(pos, genotype, subpoplist, size=1000)    # use the per-block average Fst as the Y coordinate
    fst, windows, counts = allel.windowed_patterson_fst(pos, ac1, ac2, size=1000)    # use the per-block average Fst as the Y coordinate

    y = fst

    x = [np.mean(w) for w in windows]

    # plot fst against chromosome position (bp)
    fig, ax = plt.subplots(figsize=(10, 4))
    sns.despine(ax=ax, offset=5)
    ax.plot(x, y, 'k-', lw=.5)
    ax.set_ylabel('$F_{ST}$')
    ax.set_xlabel(f'Chromosome {args.chromosome} position (bp)')
    ax.set_xlim(0, pos.max())

    print("Number of Fst values:", len(x))
    print("Number of windows:", len(y))

    # save fst figure
    timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    filename = f'fst_windowed_patterson_1000bp_{args.chromosome}_{pop1}_{pop2}.png'
    plt.savefig(filename)
    print("Saving windowed Fst plot")

    # INSPECT Fst values
    # %% Print the maximum FST value and its corresponding window
    
    print("Inspecting windowed Fst plot for maximum value")
    max_value = max(fst)


    # %% plot Fst values as a histogram to inspect Fst values
    y = fst

    # create histogram
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.hist(y, bins=30, color='blue', edgecolor='black')  # You can adjust the number of bins as needed
    ax.set_title('Histogram of Fst Values')
    ax.set_xlabel('Fst Value')
    ax.set_ylabel('Frequency')

    # save this figure
    timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    filename = f'fst_histogram_windowed_patterson_1000bp_{args.chromosome}_{pop1}_{pop2}.png'
    plt.savefig(filename)
    print("Saving histogram of Fst plot")

    # find the maximum and minimum fst values
    # Zip together windows and fst values
    zipped_windows_fst = list(zip(windows, fst))

    # Find the maximum fst value and its corresponding window
    max_fst_window = max(zipped_windows_fst, key=lambda x: x[1])
    max_window, max_fst = max_fst_window

    # Find the minimum fst value and its corresponding window
    min_fst_window = min(zipped_windows_fst, key=lambda x: x[1])
    min_window, min_fst = min_fst_window

    # Printing or further processing
    print(f"Maximum Fst value: {max_fst}, Window: {max_window}")
    print(f"Minimum Fst value: {min_fst}, Window: {min_window}")

    # set threshold for significant values as 3x the negative value, in the positive
    hist_fst_threshold = (min_fst*3)*-1
    print("Threshold for positive value being significant is:",hist_fst_threshold)

    # %% save Fst values over this threshiold
    
    hist_fst_threshold = [(window,value) for window, value in zip(windows,fst) if value >= hist_fst_threshold]

    if hist_fst_threshold:
        with open(f'hist_fst_threshold_{pop1}_{pop2}_{args.chromosome}_window_size_1000.txt', 'w') as fst_file:
            fst_file.write("Window (Genomic bps)\FST Value\n")
            for window, value in hist_fst_threshold:
                fst_file.write(f"{window[0]}-{window[1]}\t{value}\n")
        print ("Saved FST values over histogram threshold")
    else:
        print("No FST values over the histogram threshold were found")

# arguments

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate and plot FST using scikit-allel.")
    parser.add_argument("working_directory", type=str, help="Working directory for the script.")
    parser.add_argument("callset_file", type=str, help="Path to the callset file.")
    parser.add_argument("chromosome", type=str, help="Chromosome to analyze.")
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)