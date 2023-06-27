# H


import glob
import subprocess
import sys

if len(sys.argv)==1:
    print('''Example Command:
    
    python SophiesFastqMerge.py Folder1 Folder2 Folder3 etc
    
    Goal: Merge Fastq Files from Eurofins runs according to sample and pairing.
    ''')

else:
    FoldersOfInterest = sys.argv[1:]

    print('Folders of Interest:')

    # Identify Sample IDs in the original Path
    SampleIDs = []
    for Dir in FoldersOfInterest:
        print(Dir)
        for F in glob.glob(f'{Dir}/*_1.fastq.gz'):
            SampleIDs.append(F.split('/')[-1].split('_')[1])

    # Core Sample IDs 
    UniqueSampleIDs = list(set(SampleIDs))

    for SID in UniqueSampleIDs:
        
        print(SID)

        for P in ['_1.fastq.gz','_2.fastq.gz']:

            IndividualFiles = []
            for Dir in FoldersOfInterest:
                IndividualFiles += glob.glob(f'{Dir}/*{SID}*{P}')

            PreCom = ' '.join(IndividualFiles)

            Command = f'cat {PreCom} > {SID}_Combined{P}'
            subprocess.run(Command,shell=True)



    

