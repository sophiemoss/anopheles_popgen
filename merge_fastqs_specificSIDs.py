# This script needs to be run where the folders are, or put in the full pathname to find the folders
# This script looks for fastq files containing the ID that you input, so be careful if there are duplicated ids from the same sample
# forward fastq files are combined with other forward fastq files
# reverse fastq files are combined with other reverse fastq files

import glob
import subprocess
import sys
import os

# search and return fastq filepaths that match the provided criteria
# criteria include folders, subset_id and the read_suffix

def get_fastq_files(folders, subset_id, read_suffix):
    fastq_files = []
    for folder in folders:
        pattern = os.path.join(folder, '**', f'*{subset_id}*{read_suffix}.fastq.gz')
        fastq_files.extend(glob.glob(pattern, recursive=True))
    return fastq_files


# This code checks if the number of command line arguments is less than 4. If so, it means the required arguments are not provided, 
# and it prints an example command for the user. 
# Otherwise, it creates variables folders_flag and subset_id_flag which are empty, ready to store the command line flag values. 
# Then, it iterates over the command line arguments and assigns the folder names to the folders variable and the subset IDs to the subset_ids variable

if len(sys.argv) < 4:
    print('''Example Command:
python merge_fastqs_specificSIDs.py --folders Folder1 Folder2 Folder3 --subset-id eg.so1016
''')
else:
    folders_flag = "--folders"
    subset_id_flag = "--subset-id"

folders = []
subset_ids = []

for i in range(1, len(sys.argv)):
    arg = sys.argv[i]
    if arg == folders_flag:
        folders = sys.argv[i+1:]
    elif arg == subset_id_flag:
        subset_ids = sys.argv[i+1:]

# This code block checks if there are subset IDs provided. 
# If there are, it iterates over each subset ID. For each subset ID, it iterates over the read suffixes "_1" and "_2" (forward and reverse reads). 
# It calls the get_fastq_files function to obtain the list of fastq files matching the subset ID and read suffix. 
# It then prints the list of files to be merged. 
# Next, it creates the output file name based on the subset ID and read suffix. 
# It constructs the shell command to concatenate the files using cat and redirects the output to the output file. 
# Finally, it executes the command using subprocess.run and prints a message indicating the successful merge or the absence of files for that specific subset ID and read suffix.

If no subset IDs are provided, it prints a message stating that subset IDs were not provided

if subset_ids:
    for subset_id in subset_ids:
        print(f"Subset ID: {subset_id}")
        for read_suffix in ["_1", "_2"]:
            individual_files = get_fastq_files(folders, subset_id, read_suffix)
            precom = ' '.join(individual_files)

            if len(individual_files) > 0:
                print(f"Fastq files to be merged for Subset ID {subset_id}{read_suffix}:")
                for file in individual_files:
                    print(file)
                output_file = f"{subset_id}_Combined{read_suffix}.fastq.gz"
                command = f"cat {precom} > {output_file}"
                subprocess.run(command, shell=True)
                print(f"Combined Fastq files for Subset ID: {subset_id}{read_suffix}")
            else:
                print(f"No Fastq files found for Subset ID: {subset_id}{read_suffix}")
else:
    print("Subset IDs not provided.")