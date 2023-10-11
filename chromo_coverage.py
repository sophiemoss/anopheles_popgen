# example command line for this script

# python /mnt/storage8/sophie/bijagos_mosq_wgs/fastq2vcf_output_alignedtogambiae/coverageplots/chromo_coverage.py \
#       --sample 'NG-31762_Bu1002_lib648117_10128_2' \
#              --csvfile 'NG-31762_Bu1002_lib648117_10128_2.cov.filtered.txt'


import matplotlib.pyplot as plt
import pandas as pd
import glob
import argparse
from matplotlib.offsetbox import AnchoredText
import numpy as np
'''
# How to generate coverage file for different window sizes using sambamba tool.
sambamba depth window -w 20000 input.bam -o output.txt
'''
def main(args):
  Title = args.sample #What you want to call the file
  CoverageDF = pd.read_csv(args.csvfile,sep='\t') # Path to input file
  UniqueChromosome = CoverageDF['# chrom'].unique().tolist()
  fig, axs = plt.subplots(len(UniqueChromosome))
  fig.set_figheight(50)
  fig.set_figwidth(15)
  for ChromoIndex in range(len(UniqueChromosome)):
          x = CoverageDF[CoverageDF['# chrom']==UniqueChromosome[ChromoIndex]].loc[:,'chromEnd']
          y = CoverageDF[CoverageDF['# chrom']==UniqueChromosome[ChromoIndex]].loc[:,'meanCoverage']
          axs[ChromoIndex].set_title(UniqueChromosome[ChromoIndex])
          axs[ChromoIndex].plot(x, y)
          axs[ChromoIndex].axhline(y=5, color='red', linestyle='dotted')
          axs[ChromoIndex].set_ylim(0,100) #changed y axis limit
  plt.subplots_adjust(hspace = 0.8)
  plt.tight_layout()
  plt.savefig('{}.png'.format(Title))
  plt.clf()
#############################################################
parser = argparse.ArgumentParser()
parser.add_argument('--sample', type=str, help = 'title of image', required=True)
parser.add_argument('--csvfile', type=str, help = 'input csv file', required=True)
parser.set_defaults(func=main)
args = parser.parse_args()
main(args)