# This script needs to be run where the folders are, or put in the full pathname to find the folders
# This script looks for fastq files containing the ID that you input, so be careful if there are duplicated ids from the same sample
# forward fastq files are combined with other forward fastq files
# reverse fastq files are combined with other reverse fastq files
import glob
import subprocess
import sys
import os

def get_fastq_files(folders, subset_id, read_suffix):
    fastq_files = []
    for folder in folders:
        pattern = os.path.join(folder, '**', f'*{subset_id}*{read_suffix}.fastq.gz')
        fastq_files.extend(glob.glob(pattern, recursive=True))
    return fastq_files

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