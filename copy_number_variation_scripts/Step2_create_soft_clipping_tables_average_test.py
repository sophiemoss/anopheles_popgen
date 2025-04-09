####### CREATE CSV OF SOFT-CLIPPING USING IPYTHON ##########


##### directory: /mnt/storage11/sophie/bijagos_mosq_wgs/2019_melas_fq2vcf_gambiae_aligned/soft_clipping #####


# %% create samples.txt file
# Create a list of sample names
# ls *.Clipping._Normalised.csv | sed 's/.Clipping._Normalised.csv//' > samples.txt
# %%

import pandas as pd

# %% Load the discordant read guide
discordant_df = pd.read_csv('discordant_read_guide_v2.csv')

# %% Group the guide by Duplication_ID
grouped_discordant_df = discordant_df.groupby('Duplication_ID')

# %%  Read the sample names from a file
with open('samples.txt', 'r') as file:
    samples = file.read().splitlines()

# %%  Initialize an empty dictionary to hold the results
results_dict = {}

print("Iterating over dataframe")

# %%  Iterate over the grouped DataFrame
for dup_id, group in grouped_discordant_df:
    results_dict[f"{dup_id}_Pos_Range_Start"] = {}
    results_dict[f"{dup_id}_Pos_Range_End"] = {}

    for sample in samples:
        # Initialize sum for each sample for start and end positions
        results_dict[f"{dup_id}_Pos_Range_Start"][sample] = 0
        results_dict[f"{dup_id}_Pos_Range_End"][sample] = 0

        # Load the clipping data for the sample
        clipping_df = pd.read_csv(f'{sample}.Clipping._Normalised.csv')

        # Average the NormalisedClipping for the start range
        start_group = group[group['Type'] == 'Start']
        if not start_group.empty:
            start_contig = start_group['Contig'].values[0]
            start_range = start_group[['Pos_Range_Start', 'Pos_Range_End']].values[0]
            start_avg = clipping_df[(clipping_df['Contig'] == start_contig) &
                                    (clipping_df['ClipPos'] >= start_range[0]) &
                                    (clipping_df['ClipPos'] <= start_range[1])]['NormalisedClipping'].mean()
            results_dict[f"{dup_id}_Pos_Range_Start"][sample] = start_avg

        # Average the NormalisedClipping for the end range
        end_group = group[group['Type'] == 'End']
        if not end_group.empty:
            end_contig = end_group['Contig'].values[0]
            end_range = end_group[['Pos_Range_Start', 'Pos_Range_End']].values[0]
            end_avg = clipping_df[(clipping_df['Contig'] == end_contig) &
                                  (clipping_df['ClipPos'] >= end_range[0]) &
                                  (clipping_df['ClipPos'] <= end_range[1])]['NormalisedClipping'].mean()
            results_dict[f"{dup_id}_Pos_Range_End"][sample] = end_avg

print("Collected all results, now converting to dataframe")

# %%  Convert the results dictionary to a DataFrame
results_df = pd.DataFrame.from_dict(results_dict, orient='index')
# %%  Optional: if you want to convert the DataFrame such that Duplication_IDs are rows and Samples are columns
results_df = results_df.transpose()

# %%  Save the results to a CSV file
results_df.to_csv('all_genes_clipping_summary.csv', index_label='Sample')