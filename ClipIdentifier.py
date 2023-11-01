'''

--------------------
Example Command
--------------------

python ClipIdentifier.py /Path/To/SamFile /Path/To/CoverageFile

Note: The CoverageFile has been generated using 'samtools depth -aa' command

----------------------------------------
Generating SAM File with Long Reads
----------------------------------------

Note - Exclude Not Primary Alignments with F 256

samtools view -h -F 256 Po221_SpadesHybrid_Scrubbed_Krakened_No_Human.sorted.bam -o Po221_SpadesHybrid_Scrubbed_Krakened_No_Human.sorted.primary.sam

---------------
Logic Overview
---------------

1) Parse SAM File

2) Loop over sam file and extract cigar string for each read alignment

	3) Split Cigar String into Alphabetic Letters and Numbers using re.find all
	e.g. 223M1D335M1D479M3D64M1S becomes ['223', 'M', '1', 'D', '335', 'M', '1', 'D', '479', 'M', '3', 'D', '64', 'M', '1', 'S']

	4) See if soft-clipping at start and end of the read.
	If present write the chromo and position to temporary output file.

5) Read temporary output file in as dataframe 

6) Use pandas groupby function to determine the number of reads clipped 
   at each unique position. This creates condensed dataframe which now 
   has 3 columns (Contig, ClipPos and size) whereby Size equals count of reads
   in same files which are clipped at that given position.

7) Read Coverage File (Generated via samtools depth -aa)

8) Normalised number of clipped reads at given position by coverage.
   Save output as csv

-------------------------------------------

'''




import sys
import regex as re
import pandas as pd
from collections import defaultdict 


SamFile = open(sys.argv[1],'r')
CoverageFile = sys.argv[2]


OutF = sys.argv[1].replace('.sam','.Clipping.temp')
Out = open(OutF,'w')

Out.write('Contig,ClipPos\n')
for L in SamFile.readlines():
	if L[0] != '@':
		Split = L.split('\t')
		Cigar = Split[5]
		if 'S' in Cigar:


			Values = re.findall(r"[^\W\d_]+|\d+",Cigar)

			# Soft Clipping at start of Read
			if Values[1]=='S':
				# Find Soft Clipping Intercept at Start of Read in Reference
				# Contig aka First Mapped Position
				Out.write(f'{Split[2].strip()},{Split[3].strip()}\n')

			# Soft Clipping at End of Read
			if Values[-1]=='S':

				AC = int(Split[3].strip())

				# Find Soft Clipping Intercept at End of Read in Reference Contig
				# aka last mapped position

				for a in range(1,len(Values),2):

					# Match of Deletion are relative to reference
					if Values[a] in ['M','D']:

						AC+=int(Values[a-1])
				
				# Correction
				AC += -1

				Out.write(f'{Split[2].strip()},{str(AC)}\n')


df = pd.read_csv(OutF)

Condensed = df.groupby(df.columns.tolist(), as_index=False).size()
Condensed.to_csv(OutF.replace('.temp','.csv'),index=None)


Cov = pd.read_csv(CoverageFile,sep='\t',header=None)
Combined = pd.merge(Condensed,Cov,how='left',left_on=['Contig','ClipPos'],right_on=[0,1])
Combined = Combined.rename(columns={2:'Coverage','size':'ClippedCount'})

Combined = Combined[['Contig',
					 'ClipPos',
					 'ClippedCount',
					 'Coverage']]

Combined['NormalisedClipping'] = (Combined['ClippedCount'] * 100)/Combined['Coverage'] 
Combined.to_csv(OutF.replace('.temp','.Normalised.csv'),index=None)



