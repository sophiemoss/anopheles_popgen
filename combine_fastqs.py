## This script is designed for combining fastq files from the same sample, which have been generated at different times/on different runs
## The script will search --dir and any subdirectories for samples with the sample sample ID
## It will then match the fastq files and combine them, making sure that the correct fastqs combine with eachother so that all read names in the forward and reverse files are the same
## Jody Phelan created this script
## To run this script, use:
## python ~/combine_fastqs.py --dir /mnt/storage11/sophie/bijagos_mosq_wgs/mosq_raw_data/melas_2019/raw_fastq_files --r1 "(bu1035).+_1.fastq.gz" --r2 "(bu1035).+_2.fastq.gz"
## replace bu1035 with the sample of interest that needs combining


import re
import sys
import os
from collections import defaultdict
import argparse
import subprocess as sp

parser = argparse.ArgumentParser(description='Find fastq files in a directory and return a list of tuples with the sample name and the path to the fastq files from both pairs.')
parser.add_argument('--dir',nargs="+", type=str, help='Directory to search for fastq files')
parser.add_argument('--r1', metavar='r1_pattern', type=str, help='Pattern to match R1 files')
parser.add_argument('--r2', metavar='r2_pattern', type=str, help='Pattern to match R2 files')
parser.add_argument('--dry-run',action='store_true')
args = parser.parse_args()

class Sample:
    def __init__(self,prefix,r1,r2):
        self.prefix = prefix
        self.r1 = sorted(r1)
        self.r2 = sorted(r2)
        if len(r1)!=len(r2):
            raise ValueError("Number of R1 and R2 files do not match")
        if len(r1)==1:
            self.multiple = False
        else:
            self.multiple = True

    def __repr__(self):
        return f"Sample(prefix={self.prefix},r1={self.r1},r2={self.r2})"


def sort_out_paried_files(files,r1_suffix,r2_suffix):
    prefixes = defaultdict(lambda:{"r1":[],"r2":[]})

    for f in files:
        tmp1 = re.search("%s" % r1_suffix,f)
        tmp2 = re.search("%s" % r2_suffix,f)
        p = None
        if tmp1:
            p = tmp1.group(1).split("/")[-1]
            prefixes[p]['r1'].append(f)
        elif tmp2:
            p = tmp2.group(1).split("/")[-1]
            prefixes[p]['r2'].append(f)

    runs = []
    for p,vals in prefixes.items():
        if len(vals['r1'])!=1 or len(vals['r2'])!=1:
            if len(vals['r1'])!=len(vals['r2']):
                raise ValueError(f"Number of R1 and R2 files for sample {p} do not match")
            vals['r1'].sort()
            vals['r2'].sort()
            runs.append(
                Sample(p,vals['r1'],vals['r2'])
            )
        else:
            runs.append(
                Sample(p,vals['r1'],vals['r2'])
            )
    return runs

def find_fastq_files(directories,r1_pattern,r2_pattern):
    """
    Find fastq files in a directory and return a
    list of tuples with the sample name and the
    path to the fastq files from both pairs.
    """
    files = []
    for d in directories:
        for a,b,c in os.walk(d):
            for f in c:
                files.append(f"{os.path.abspath(a)}/{f}")
#    print(files)
    fastq_files = sort_out_paried_files(files,r1_pattern,r2_pattern)

    return fastq_files


with open("combine.log","w") as O:
    for s in find_fastq_files(args.dir,args.r1,args.r2):
        combined_file_1 = s.prefix+"_1.fastq.gz"
        combined_file_2 = s.prefix+"_2.fastq.gz"
        O.write("\t".join([s.prefix,combined_file_1,combined_file_2,",".join(s.r1),",".join(s.r2)])+"\n")
        cmd = f"cat {' '.join(s.r1)} > {combined_file_1}"
        sys.stdout.write(f"Running: {cmd}\n")
        if not args.dry_run:
            sp.call(cmd,shell=True)
        cmd = f"cat {' '.join(s.r2)} > {combined_file_2}"
        sys.stdout.write(f"Running: {cmd}\n")
        if not args.dry_run:
            sp.call(cmd,shell=True)
