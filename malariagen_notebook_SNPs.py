## Example MalariaGEN analysis: SNPs https://malariagen.github.io/vector-data/ag3/examples/snps.html
pip install -qU malariagen_data seaborn
import malariagen_data

# %%
ag3 = malariagen_data.Ag3()
ag3

seq = ag3.genome_sequence("2R").compute()
seq[:100]

#%%
from collections import Counter
from functools import lru_cache
import numpy as np
import pandas as pd
import dask
import dask.array as da
# silence some dask warnings
dask.config.set(**{'array.slicing.split_large_chunks': True})
import matplotlib as mpl
import matplotlib.pyplot as plt
%matplotlib inline
import seaborn as sns
import allel
from dask.diagnostics.progress import ProgressBar
import malariagen_data

ag3 = malariagen_data.Ag3()
ag3

def plot_sequence_composition(seq_id, window_size=100_000, ax=None):
    
    # load reference sequence
    seq = ag3.genome_sequence(seq_id).compute()
    
    if ax is None:
        # make the figure size relative to largest contig
        figw = 7 * len(seq) / len(ag3.genome_sequence("2R"))
        fig, ax = plt.subplots(figsize=(figw, 1))
 
    # convert to upper-case
    seq = np.char.upper(seq)

    # locate nucleotides
    is_a = seq == b'A'
    is_c = seq == b'C'
    is_g = seq == b'G'
    is_t = seq == b'T'
    is_n = seq == b'N'
    # check there's nothing unexpected
    is_other = ~is_a & ~is_c & ~is_g & ~is_t & ~is_n
    assert np.sum(is_other) == 0

    # construct windows
    bins = np.arange(0, len(seq), window_size)

    # count nucleotides in each window
    h_a, _ = np.histogram(np.nonzero(is_a)[0] + 1, bins=bins)
    h_c, _ = np.histogram(np.nonzero(is_c)[0] + 1, bins=bins)
    h_g, _ = np.histogram(np.nonzero(is_g)[0] + 1, bins=bins)
    h_t, _ = np.histogram(np.nonzero(is_t)[0] + 1, bins=bins)
    h_n, _ = np.histogram(np.nonzero(is_n)[0] + 1, bins=bins)

    # plot frequence of nucleotides within each bin
    left = bins[:-1]
    bottom = 0
    width = np.diff(bins)
    palette = sns.color_palette('colorblind')
    colors = [palette[i] for i in [2, 0, 3, 8]] + ['k']
    for h, c, l in zip([h_a, h_t, h_g, h_c, h_n], colors, 'ATGCN'):
        ax.bar(left, h, width=width, bottom=bottom, color=c, align='edge', label=l)
        bottom += h
        
    # tidy up plot
    ax.set_xlim(0, len(seq))
    ax.set_yticks(ax.get_ylim())
    ax.set_yticklabels(['0%', '100%'])
    ax.set_title(seq_id)
    # convert X axis to Mbp
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: x//1e6))
    ax.set_xlabel("Position (Mbp)")
    
    # add centromere/telomere annotation
    if seq_id in {'2L', '3L'}:
        ltxt = "centromere"
        rtxt = "telomere"
    else:
        ltxt = "telomere"
        rtxt = "centromere"
    ax.annotate(ltxt, xy=(0, 1), xycoords="axes fraction", 
                xytext=(0, 2), textcoords='offset points', va='bottom', ha='left')
    ax.annotate(rtxt, xy=(1, 1), xycoords="axes fraction", 
                xytext=(0, 2), textcoords='offset points', va='bottom', ha='right')
    
    # add legend - reverse order so matches the plot
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(reversed(handles), reversed(labels), 
              loc='center left', bbox_to_anchor=(1, .5), 
              prop=dict(family='monospace'), ncol=1,
              frameon=False)
        
plot_sequence_composition('2R')

# %%
plot_sequence_composition('2L')
# %%
plot_sequence_composition('X')

# %%
ds_snps = ag3.snp_calls(region="2R", sample_sets="3.0")
ds_snps


# %%
# variant positions here holds the 1 based positions (coordinates) where the SNPs were called
# eg. first ten positions

pos = ds_snps["variant_position"].data
pos[:10].compute()

# %%
# these positions do not increase by 1 the whole way through the genome, because
# there are gaps in the reference genome. Visualise that by plotting distance between adjacent sites:

def plot_gaps(contig):
    ds_snps = ag3.snp_calls(region=contig)
    pos = ds_snps["variant_position"].data.compute()
    d = np.diff(pos)
    loc_gap = d > 1
    x = pos[1:][loc_gap]
    y = d[loc_gap]
    fig, ax = plt.subplots(figsize=(7, 2))
    ax.plot(x, y, color='black', linewidth=.5, marker='o', linestyle=' ', markersize=2, mfc='none')
    ax.set_title("2R")
    # convert X axis to Mbp
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: x//1e6))
    ax.set_xlabel("Position (Mbp)")
    ax.set_ylabel('Gap size (bp)')

# %%
plot_gaps(contig="2R")

#%%
len(seq) - len(pos)

# %%
alleles = ds_snps["variant_allele"].data
alleles[:10].compute()

# %%
seq[pos[:10].compute() - 1]


# access SNP genotypes for chromosome arm 2R for all wild samples in Ag3

# %%
ds_snps = ag3.snp_calls(region="2R", sample_sets="3.0")
gt = ds_snps["call_genotype"].data
gt.shape
# %%
gt[9, 0].compute()
gt[9, 2392].compute()
gt[9, 248].compute()