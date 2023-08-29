## download additional samples from ncbi using fasterq-dump
conda activate fasterq-dump
fasterq-dump SRR000001

## to make a genomics database of sample VCFs, use the following
ls *.g.vcf.gz | sed 's/.g.vcf.gz//' > fastq2vcfsamples.txt

/mnt/storage11/sophie/fastq2matrix/scripts/merge_vcfs.py import --sample-file fastq2vcfsamples.txt --ref GCF_000005575.2_AgamP3_genomic.fasta --prefix gambiae_bijagos_2022 --vcf-dir .

## UP TO HERE, THEN NEED TO MERGE VCF AS BELOW AND THEN GO TO COMINE WITH MALARIAGEN SH TO COMBINE WITH THE MALARIAGEN VCF

## now merge VCF files

merge_vcfs.py genotype --ref GCF_000005575.2_AgamP3_genomic.fasta --prefix gambiae_bijagos_2022 > mergevcf_log.txt 2>&1

# resulting vcf is called gambiae_bijagos_2022.2023_07_25.genotyped.vcf.gz

# concerned that this genotyping does not use a set of validated variants like I had for plasmodium
# this has a genotyping pipeline:
# https://github.com/malariagen/pipelines/blob/v0.0.4/docs/specs/snp-genotyping-vector.md
# pipeline being sent from MalariaGEN.