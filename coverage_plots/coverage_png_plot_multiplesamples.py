# usage
# python coverage_png_plot_multisamples_all.py --sample '<sample_name>' --csvfile 'csv1.csv,csv2.csv,csvN.csv'
# set the y limit if the plots are cut off at the top 

import matplotlib.pyplot as plt
import pandas as pd
import argparse

def main(args):
    Title = args.sample  # Title for the output plot
    csv_files = args.csvfiles.split(',')  # Split input CSV files by comma
    
    UniqueChromosomes = set()  # Store unique chromosome names
    
    for csv_file in csv_files:
        CoverageDF = pd.read_csv(csv_file, sep='\t')  # Read CSV file
        UniqueChromosomes.update(CoverageDF['# chrom'].unique())  # Add unique chromosome names
        
    for chromosome in UniqueChromosomes:
        fig, ax = plt.subplots()  # Create a single plot for each chromosome
        
        for csv_file in csv_files:
            CoverageDF = pd.read_csv(csv_file, sep='\t')  # Read CSV file
            chromosome_data = CoverageDF[CoverageDF['# chrom'] == chromosome]
            x = chromosome_data['chromEnd']
            y = chromosome_data['meanCoverage']
            
            ax.plot(x, y, label=f'{csv_file}')
        
        ax.axhline(y=5, color='red', linestyle='dotted')
        ax.set_ylim(0, 9000) # made this very high for mitochondrial genome, make this smaller for other chromosomes
        ax.set_title(f"Coverage Plot - {chromosome}")
        ax.set_xlabel("Position")
        ax.set_ylabel("Mean Coverage")
        ax.legend()  # Add a legend to distinguish CSV files
        
        plt.tight_layout()
        plt.savefig(f'{Title}_{chromosome}.png')
        plt.clf()

parser = argparse.ArgumentParser()
parser.add_argument('--sample', type=str, help='title of images', required=True)
parser.add_argument('--csvfiles', type=str, help='input csv files separated by commas', required=True)
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)